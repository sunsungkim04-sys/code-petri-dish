"""코드 배양접시 — 되감기 2 판정. 사전등록 `되감기2-사전등록-2026-09-15.md` §2~§3 을 옮겼다.
시드 101~200 결과 0건 시점(2026-09-15)에 썼다. 판정선을 바꾸면 노트에 개정으로 적는다.

  Q2' 도착점   정체 층(안정성 먼저 → 가장 흔한 끝 우세 코드 공유 비율) · 이력 층(자기 과거 순위 r̄ 대 균등 귀무)
  Q3' 길       길 층(가장 흔한 길 비율) · 순서 층(함께 거친 교대 코드 쌍의 앞뒤 일치율 C 대 뒤섞기 귀무 · 파생 쌍 제외판)
  Q4' 조건 의존 p · q · r̄ 의 조건 차이 — 접시 재추출 부트스트랩
  Q1 재현은 이 스크립트 밖(main_analyze.py 를 같은 정의로).

계획: 조건 3(space · energy · mat 밀도 8) × 시드 101~200 × 50,000틱 · sim v0.3.1 · μ 0.01
  - 계획한 300접시가 다 모이지 않거나 설정이 다르거나 재료 보존 위반이 있으면 판정하지 않는다
  - 값이 없으면(생존 부족 · 쌍 0) '판정 불가' — NaN 을 무효과로 읽지 않는다(본실험 🔻 결함)
  - 백분위는 선형 보간 · 난수 random.Random(20260916) 고정 · 순수 파이썬 3.8
실행: python3 rewind2_analyze.py [petri 폴더, 기본 ~/petri]  →  표 인쇄 + _RESULT_rewind2.json
드라이런(판정 아님 · 코드 경로 확인용, 출력은 버린다):
  PETRI_DRYRUN_DIR=pilot2 PETRI_DRYRUN_SEEDS=9001-9020 PETRI_DRYRUN_MINSURV=10 python3 rewind2_analyze.py
"""
import json
import math
import os
import random
import sys
from collections import Counter
from itertools import combinations

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
_dry_dir = os.environ.get("PETRI_DRYRUN_DIR")
_dry_seeds = os.environ.get("PETRI_DRYRUN_SEEDS")
DRY = bool(_dry_dir)
FOLDER = _dry_dir if DRY else "rewind2"
SEEDS = list(range(int(_dry_seeds.split("-")[0]), int(_dry_seeds.split("-")[1]) + 1)) if DRY else list(range(101, 201))
VERSIONS = {"0.3.0", "0.3.1"} if DRY else {"0.3.1"}
TICKS, MU, MAT_DENSITY = 50000, 0.01, 8
CONDS = ["space", "energy", "mat"]
MIN_SURV = int(os.environ.get("PETRI_DRYRUN_MINSURV", "30")) if DRY else 30
LAG_FROM, END = 40000, 50000
STABLE, ONE_PLACE, SPLIT = 0.90, 0.90, 0.50
SAME_PATH, RARE_PATH = 0.90, 0.10
N_RANK_NULL, N_SHUF, N_BOOT = 10000, 1000, 2000
RNG = random.Random(20260916)
SEED_CODE = "rascld"


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def bray(a, b):
    keys = set(a) | set(b)
    den = sum(a.get(k, 0) + b.get(k, 0) for k in keys)
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys) / den if den else float("nan")


def one_edit(a, b):
    if a == b:
        return True
    if abs(len(a) - len(b)) > 1:
        return False
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) == 1
    short, long_ = (a, b) if len(a) < len(b) else (b, a)
    return any(short == long_[:i] + long_[i + 1:] for i in range(len(long_)))


# ---------------- 자료 · 완결성 가드 ----------------
dishes = {c: [] for c in CONDS}
problems = []
for c in CONDS:
    for s in SEEDS:
        f = os.path.join(BASE, FOLDER, f"{c}_d8_mu0p01_s{s:05d}.json")
        if not os.path.exists(f):
            problems.append(f"없음 {os.path.basename(f)}")
            continue
        z = json.load(open(f))
        if (z["sim_version"] not in VERSIONS or z["ticks_planned"] != TICKS or z["opts"]["mu"] != MU
                or z["opts"]["density"] != MAT_DENSITY or z["cond"] != c or z["opts"]["seed"] != s):
            problems.append(f"설정 다름 {os.path.basename(f)}")
        if z["conservation_violations"]:
            problems.append(f"재료 보존 위반 {os.path.basename(f)}")
        dishes[c].append(z)
print(f"{'🟡 드라이런 — 판정 아님 · ' if DRY else ''}폴더 {FOLDER} · 시드 {SEEDS[0]}~{SEEDS[-1]}")
for c in CONDS:
    print(f"  {c:>6}: {len(dishes[c])}/{len(SEEDS)}")
if problems:
    print(f"🚨 완결성 가드 {len(problems)}건 — 판정하지 않는다")
    for p in problems[:20]:
        print("   ", p)
    sys.exit(1)


def tops_at(z, t):
    for tt, top in z["tops"]:
        if tt == t:
            return top
    return None


def pop_at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r[1]
    return None


def comp(z, t):
    top, pop = tops_at(z, t), pop_at(z, t)
    return {k: cnt / pop for k, cnt in top} if top and pop else None


def dominant(z):
    return z["final_counts"][0][0]


def stability(z):
    fin = dominant(z)
    recs = [top[0][0] for t, top in z["tops"] if t >= LAG_FROM and top]
    return sum(k == fin for k in recs) / len(recs) if recs else float("nan")


def path(z):
    ev, order, prev, prev_t, run = {}, [], None, None, 0
    for t, top in z["tops"]:
        k = top[0][0] if top else None
        if k is None:
            prev, run = None, 0
            continue
        run = run + 1 if k == prev else 1
        if run == 2 and k != SEED_CODE and k not in ev:
            ev[k] = prev_t
            order.append(k)
        prev, prev_t = k, t
    return tuple(order), ev


def concord_items(evs, exclude_derived):
    items = []
    for i, j in combinations(range(len(evs)), 2):
        common = sorted(set(evs[i]) & set(evs[j]))
        for u, v in combinations(common, 2):
            if exclude_derived and one_edit(u, v):
                continue
            items.append((i, j, u, v))
    return items


def concordance(evs, items):
    conc = tot = 0
    for i, j, u, v in items:
        a = evs[i][u] - evs[i][v]
        b = evs[j][u] - evs[j][v]
        if a == 0 or b == 0:
            continue
        tot += 1
        conc += (a > 0) == (b > 0)
    return conc, tot


def shuffled(evs):
    out = []
    for e in evs:
        codes, times = list(e.keys()), list(e.values())
        RNG.shuffle(times)
        out.append(dict(zip(codes, times)))
    return out


def order_test(evs, exclude_derived):
    items = concord_items(evs, exclude_derived)
    conc, tot = concordance(evs, items)
    if tot == 0:
        return {"items": len(items), "tot": 0, "C": float("nan"), "null975": float("nan"), "verdict": "판정 불가 — 잴 코드 쌍이 없다"}
    null = []
    for _ in range(N_SHUF):
        cc, tt = concordance(shuffled(evs), items)
        null.append(cc / tt if tt else float("nan"))
    C, n975 = conc / tot, pct(null, .975)
    # 드라이런에서 찾아 결과 전에 고침(09-15): 쌍이 너무 적으면 귀무 97.5% 가 천장 1.0 에 닿아 C 가 넘을 수 없다 → '반복 없음' 이 아니라 판정 불가
    if n975 >= 1.0:
        verdict = "판정 불가 — 코드 쌍이 너무 적어 귀무가 천장(1.0)에 닿는다"
    else:
        verdict = "순서가 반복된다" if C > n975 else "순서 반복 없음"
    return {"items": len(items), "tot": tot, "C": C, "null975": n975, "null50": pct(null, .5), "verdict": verdict}


def testable(t):
    return not t["verdict"].startswith("판정 불가")


# ---------------- 조건마다 ----------------
res = {}
print("\n== 생존")
for c in CONDS:
    surv = [z for z in dishes[c] if z["extinct_at"] < 0]
    print(f"  {c:>6}: {len(surv)}/{len(dishes[c])}")
    res[c] = {"n": len(surv), "surv": surv}

print("\n== Q2' 도착점 — 정체 층 (안정성 중앙 < 0.90 이면 떠돈다 · 아니면 공유 비율 Wilson 하한 ≥ 0.90 한 곳 · 상한 ≤ 0.50 갈라짐 · 그 밖 대체로 한 곳)")
for c in CONDS:
    r, surv, n = res[c], res[c]["surv"], res[c]["n"]
    if n < MIN_SURV:
        r["q2_id"] = {"verdict": f"판정 불가 — 생존 {n} < {MIN_SURV}"}
        print(f"  {c:>6}: {r['q2_id']['verdict']}")
        continue
    stab = [stability(z) for z in surv]
    doms = [dominant(z) for z in surv]
    modal, k = Counter(doms).most_common(1)[0]
    lo, hi = wilson(k, n)
    stab_med = pct(stab, .5)
    if stab_med < STABLE:
        verdict = "끝 우세 코드가 자리 잡지 않는다(떠돈다) — 공유 비율로 판정하지 않음"
    elif lo >= ONE_PLACE:
        verdict = "한 곳에 도착한다"
    elif hi <= SPLIT:
        verdict = "여러 곳으로 갈라진다"
    else:
        verdict = "대체로 한 곳 · 갈래가 있다"
    r["doms"] = doms
    r["q2_id"] = {"stability_median": stab_med, "stability_lt1": sum(x < 1 for x in stab), "modal": modal, "k": k,
                  "p": k / n, "lo": lo, "hi": hi, "top5": Counter(doms).most_common(5), "verdict": verdict}
    print(f"  {c:>6}: 안정성 중앙 {stab_med:.3f} (1 미만 {r['q2_id']['stability_lt1']}/{n}) · 가장 흔한 끝 우세 {modal} {k}/{n} = {k / n:.0%} [{lo:.0%}, {hi:.0%}] → {verdict}")
    print(f"          끝 우세 상위: " + " · ".join(f"{a} {b}" for a, b in Counter(doms).most_common(5)))

print("\n== Q2' 도착점 — 이력 층 (끝 조성 대 40,000틱 조성 Bray–Curtis · 자기 과거 순위 r̄ < 균등 귀무 2.5% 면 이력이 남는다)")
for c in CONDS:
    r, surv, n = res[c], res[c]["surv"], res[c]["n"]
    if n < MIN_SURV:
        r["q2_hist"] = {"verdict": f"판정 불가 — 생존 {n} < {MIN_SURV}"}
        print(f"  {c:>6}: {r['q2_hist']['verdict']}")
        continue
    ends = [comp(z, END) for z in surv]
    pasts = [comp(z, LAG_FROM) for z in surv]
    if any(e is None for e in ends) or any(p is None for p in pasts):
        r["q2_hist"] = {"verdict": "판정 불가 — 조성 기록 없음"}
        print(f"  {c:>6}: {r['q2_hist']['verdict']}")
        continue
    ranks, own_d, other_d = [], [], []
    for i in range(n):
        own = bray(ends[i], pasts[i])
        others = [bray(ends[i], pasts[j]) for j in range(n) if j != i]
        ranks.append((sum(d < own for d in others) + 0.5 * sum(d == own for d in others)) / (n - 1))
        own_d.append(own)
        other_d.append(pct(others, .5))
    rbar = sum(ranks) / n
    null = [sum(RNG.randrange(n) for _ in range(n)) / (n * (n - 1)) for _ in range(N_RANK_NULL)]
    n025 = pct(null, .025)
    verdict = "접시마다 제 이력이 1만 틱 넘게 남는다" if rbar < n025 else "1만 틱이면 이력이 지워진다(접시끼리 바꿔 끼울 수 있다)"
    r["ranks"] = ranks
    r["q2_hist"] = {"rbar": rbar, "null025": n025, "nearest_own": sum(x == 0 for x in ranks), "own_median": pct(own_d, .5),
                    "other_median": pct(other_d, .5), "verdict": verdict}
    print(f"  {c:>6}: r̄ {rbar:.3f} · 귀무 2.5% {n025:.3f} · 자기 과거가 가장 가까운 접시 {r['q2_hist']['nearest_own']}/{n} · "
          f"거리 중앙 자기 {pct(own_d, .5):.3f} 대 남 {pct(other_d, .5):.3f} → {verdict}")

print("\n== Q3' 길 — 길 층 (가장 흔한 길을 그대로 걸은 비율 Wilson 하한 ≥ 0.90 매번 · 상한 ≤ 0.10 드묾 · 그 밖 우연)")
for c in CONDS:
    r, surv, n = res[c], res[c]["surv"], res[c]["n"]
    if n < MIN_SURV:
        r["q3_path"] = {"verdict": f"판정 불가 — 생존 {n} < {MIN_SURV}"}
        print(f"  {c:>6}: {r['q3_path']['verdict']}")
        continue
    paths = [path(z) for z in surv]
    r["paths"] = paths
    pc = Counter(p for p, _ in paths)
    mpath, kq = pc.most_common(1)[0]
    lo, hi = wilson(kq, n)
    verdict = "매번 같은 길" if lo >= SAME_PATH else ("같은 길은 드물다" if hi <= RARE_PATH else "같은 길은 우연에 달렸다")
    lens = Counter(len(p) for p, _ in paths)
    r["q3_path"] = {"modal": list(mpath), "k": kq, "q": kq / n, "lo": lo, "hi": hi, "distinct": len(pc),
                    "len_dist": dict(sorted(lens.items())), "verdict": verdict}
    shown = " → ".join(mpath) if mpath else "(교대 없음)"
    print(f"  {c:>6}: 가장 흔한 길 [{shown}] {kq}/{n} = {kq / n:.0%} [{lo:.0%}, {hi:.0%}] · 서로 다른 길 {len(pc)} · 길이 분포 {dict(sorted(lens.items()))} → {verdict}")

print("\n== Q3' 길 — 순서 층 (함께 거친 교대 코드 쌍의 앞뒤 일치율 C > 뒤섞기 귀무 97.5% 면 반복 · 둘이 갈리면 파생 제외판)")
for c in CONDS:
    r, n = res[c], res[c]["n"]
    if n < MIN_SURV:
        r["q3_order"] = {"verdict": f"판정 불가 — 생존 {n} < {MIN_SURV}"}
        print(f"  {c:>6}: {r['q3_order']['verdict']}")
        continue
    evs = [ev for _, ev in r["paths"]]
    full = order_test(evs, False)
    excl = order_test(evs, True)
    if not testable(excl):
        final = full["verdict"] + (" ⚠️ 파생 제외판은 판정 불가 — 파생 결합을 걷어내지 못했다" if testable(full) else "")
    elif not testable(full) or full["verdict"] == excl["verdict"]:
        final = excl["verdict"]
    else:
        final = excl["verdict"] + " (파생 제외판을 따름)"
    r["q3_order"] = {"full": full, "excl": excl, "verdict": final}
    for name, t in (("전체", full), ("파생 제외", excl)):
        print(f"  {c:>6} {name:>5}: 코드 쌍 {t['tot']:,} · C {t['C']:.3f} · 귀무 97.5% {t['null975']:.3f} → {t['verdict']}")
    print(f"  {c:>6} → {final}")

print("\n== Q4' 조건 의존 — 접시 재추출 부트스트랩 2,000회 · 차이 95% 구간이 0 을 벗어나면 다르다")


def modal_share(xs):
    return Counter(xs).most_common(1)[0][1] / len(xs)


ok = [c for c in CONDS if res[c]["n"] >= MIN_SURV]
stats = {
    "Q2' 끝 우세 공유 몫 p": lambda c, idx: modal_share([res[c]["doms"][i] for i in idx]),
    "Q3' 가장 흔한 길 몫 q": lambda c, idx: modal_share([res[c]["paths"][i][0] for i in idx]),
    "Q2' 자기 과거 순위 r̄": lambda c, idx: sum(res[c]["ranks"][i] for i in idx) / len(idx),
}
boots = {c: [[RNG.randrange(res[c]["n"]) for _ in range(res[c]["n"])] for _ in range(N_BOOT)] for c in ok}
q4 = {}
for a, b in combinations(CONDS, 2):
    for name, fn in stats.items():
        label = f"{a} − {b}"
        if a not in ok or b not in ok or any(key not in res[x] for x in (a, b) for key in ("doms", "paths", "ranks")):
            q4[f"{label} {name}"] = {"verdict": "판정 불가 — 비교할 값 없음"}
            print(f"  {label:>16} {name}: 판정 불가 — 비교할 값 없음")
            continue
        full_a = fn(a, range(res[a]["n"]))
        full_b = fn(b, range(res[b]["n"]))
        diffs = [fn(a, ia) - fn(b, ib) for ia, ib in zip(boots[a], boots[b])]
        lo, hi = pct(diffs, .025), pct(diffs, .975)
        verdict = "다르다" if (lo > 0 or hi < 0) else "차이 검출 안 됨"
        q4[f"{label} {name}"] = {"a": full_a, "b": full_b, "lo": lo, "hi": hi, "verdict": verdict}
        print(f"  {label:>16} {name}: {full_a:.3f} 대 {full_b:.3f} · 차이 [{lo:+.3f}, {hi:+.3f}] → {verdict}")

out = {"dry": DRY, "folder": FOLDER, "seeds": [SEEDS[0], SEEDS[-1]], "q4": q4, "conds": {}}
for c in CONDS:
    out["conds"][c] = {k: v for k, v in res[c].items() if k in ("n", "q2_id", "q2_hist", "q3_path", "q3_order")}
name = "_RESULT_rewind2_dryrun.json" if DRY else "_RESULT_rewind2.json"
json.dump(out, open(os.path.join(BASE, name), "w"), ensure_ascii=False, indent=1)
print(f"\n저장: {name}")
