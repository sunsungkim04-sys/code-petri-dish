"""코드 배양접시 — 테이프 되감기 본실험 판정.

사전등록 `되감기-사전등록-2026-09-15.md` 의 §2 개정판(조건별 중립 문턱) · §3(Q1~Q4) 을 그대로 옮겼다.
본실험 자료 0건 시점(2026-09-15)에 썼다. 판정선을 바꾸면 노트에 개정으로 적는다.

🔻 09-15 결과 후 교정(판정선은 그대로, 스크립트 결함만): Q3 에서 공통 사건 ≥ 3 인 쌍이 0 이면 τ 가 NaN 인데
   `NaN > 문턱` 이 거짓이 되어 '순서 반복 없음' 으로 찍혔다 — 데이터 없음을 무효과로 읽은 것.
   → 쌍이 없으면 '판정 불가', Q4 에서 비교할 값이 없으면 '판정 불가' 로 낸다.

계획: 조건 3(space · energy · mat 밀도 8) × 시드 1~100 × 50,000틱 · sim v0.3.0 · μ 0.01
  - 계획한 300접시가 다 모이지 않거나, 계획 밖 설정 파일이 섞였거나, 재료 보존 위반이 있으면 **판정하지 않고 끝낸다**
  - 백분위는 선형 보간(numpy 기본과 같음) · 난수는 random.Random(20260915) 고정 · 순수 파이썬 3.8
실행: python3 main_analyze.py [결과 폴더, 기본 ~/petri/main]  →  표 인쇄 + ../_RESULT_main.json
"""
import glob
import json
import math
import os
import random
import statistics as st
import sys
from itertools import combinations

DIR = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri/main")
CONDS = ["space", "energy", "mat"]
# 드라이런(판정 아님): PETRI_DRYRUN_SEEDS="9001-9020" PETRI_DRYRUN_MINSURV=10 — 파일럿 자료로 코드 경로만 확인하고 출력은 버린다
_dry = os.environ.get("PETRI_DRYRUN_SEEDS")
# 되감기 2 의 Q1 재현(09-15 추가 · 판정선 불변 · 기본값이면 본실험과 같은 출력):
#   PETRI_SEEDS="101-200" PETRI_VERSION=0.3.1 PETRI_OUT=_RESULT_rewind2_q1.json python3 main_analyze.py ~/petri/rewind2
_seeds = _dry or os.environ.get("PETRI_SEEDS")
SEEDS = list(range(int(_seeds.split("-")[0]), int(_seeds.split("-")[1]) + 1)) if _seeds else list(range(1, 101))
TICKS, VERSION, MU, MAT_DENSITY = 50000, os.environ.get("PETRI_VERSION", "0.3.0"), 0.01, 8
OUT = os.environ.get("PETRI_OUT", "_RESULT_main.json")
N_NULL, N_BOOT = 1000, 2000
MIN_SURV = int(os.environ.get("PETRI_DRYRUN_MINSURV", "30")) if _dry else 30
RNG = random.Random(20260915)
EVENTS = ["E_j", "E_noc", "E_e", "E_h", "E_short", "E_seed"]
TAG = {"E_j": "j", "E_noc": "noc", "E_e": "e", "E_h": "h"}
NEUTRAL_H = {"space", "mat"}          # 이 조건에서 h 는 쉬기와 같다 → E_h 는 사건이 아니라 두 번째 중립 대조
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


def cosine(a, b):
    dot = sum(v * b.get(k, 0) for k, v in a.items())
    na = math.sqrt(sum(v * v for v in a.values()))
    nb = math.sqrt(sum(v * v for v in b.values()))
    return dot / (na * nb) if na and nb else float("nan")


def kmer_vec(entries, k=3):
    d = {}
    for key, cnt in entries:
        for i in range(len(key) - k + 1):
            s = key[i:i + k]
            d[s] = d.get(s, 0) + cnt
    return d


def letter_vec(entries):
    d = {}
    for key, cnt in entries:
        for ch in key:
            d[ch] = d.get(ch, 0) + cnt
    return d


def tau_b(x, y):
    n0 = n1 = n2 = s = 0
    m = len(x)
    for i in range(m):
        for j in range(i + 1, m):
            n0 += 1
            dx, dy = x[i] - x[j], y[i] - y[j]
            if dx == 0:
                n1 += 1
            if dy == 0:
                n2 += 1
            if dx * dy > 0:
                s += 1
            elif dx * dy < 0:
                s -= 1
    den = math.sqrt((n0 - n1) * (n0 - n2))
    return s / den if den > 0 else None


def fmt_p(k, n):
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {k / n:.0%} [{lo:.0%}, {hi:.0%}]" if n else "—"


# ---------------------------------------------------------------- 불러오기 · 완결성
data = {c: {} for c in CONDS}
problems = []
for f in glob.glob(os.path.join(DIR, "*.json")):
    z = json.load(open(f))
    c, s = z["cond"], z["opts"]["seed"]
    name = os.path.basename(f)
    if c not in CONDS or s not in SEEDS:
        problems.append(f"계획 밖 조건/시드: {name}")
        continue
    if z["sim_version"] != VERSION or z["ticks_planned"] != TICKS or abs(z["opts"]["mu"] - MU) > 1e-12:
        problems.append(f"계획 밖 설정(sim {z['sim_version']} · 틱 {z['ticks_planned']} · μ {z['opts']['mu']}): {name}")
        continue
    if c == "mat" and z["opts"]["density"] != MAT_DENSITY:
        problems.append(f"계획 밖 재료 밀도 {z['opts']['density']}: {name}")
        continue
    if s in data[c]:
        problems.append(f"중복: {name}")
        continue
    data[c][s] = z

print(f"폴더 {DIR}")
missing = {c: [s for s in SEEDS if s not in data[c]] for c in CONDS}
for c in CONDS:
    print(f"  {c:>6}: {len(data[c])}/{len(SEEDS)}" + (f"  ← 빠진 시드 {missing[c][:10]}{'…' if len(missing[c]) > 10 else ''}" if missing[c] else ""))
cons = sum(z["conservation_violations"] for z in data["mat"].values())
checks = sum(z["conservation_checks"] for z in data["mat"].values())
print(f"  재료 보존 위반 {cons}/{checks}")
if problems:
    print("🚨 계획 밖 파일:", *problems[:10], sep="\n  ")
if any(missing.values()) or problems or cons:
    print("\n판정하지 않는다 — 계획한 300접시가 다 모이지 않았거나, 계획 밖 파일이 섞였거나, 재료 보존 위반이 있다.")
    sys.exit(1)


# ---------------------------------------------------------------- 접시별 양
def specials_series(z):
    pop_at = {r[0]: r[1] for r in z["samples"]}
    hdr = z["special_header"]
    out = {}
    for tag in hdr:
        col = hdr.index(tag) + 1
        out[tag] = [(row[0], (row[col][1] / pop_at[row[0]]) if (row[col] and pop_at.get(row[0])) else 0.0) for row in z["specials"]]
    return out


def l2(series):
    return max((min(series[i][1], series[i + 1][1]) for i in range(len(series) - 1)), default=0.0)


def first_exceed(series, theta):
    for i in range(len(series) - 1):
        if min(series[i][1], series[i + 1][1]) > theta:
            return series[i][0]
    return None


def e_short(z):
    ix = {h: i for i, h in enumerate(z["header"])}
    run, start = 0, None
    for r in z["samples"]:
        if r[ix["pop"]] > 0 and r[ix["mean_len"]] <= 5.0:
            if run == 0:
                start = r[0]
            run += 1
            if run >= 10:
                return start
        else:
            run = 0
    return None


def e_seed(z):
    for tick, top in z["tops"]:
        if SEED_CODE not in {k for k, _ in top}:
            return tick
    return None


result = {"plan": dict(conds=CONDS, seeds=[SEEDS[0], SEEDS[-1]], ticks=TICKS, sim=VERSION, mu=MU, mat_density=MAT_DENSITY), "conds": {}}
surv, ev, theta = {}, {}, {}
for c in CONDS:
    zs = [data[c][s] for s in SEEDS]
    surv[c] = [z for z in zs if z["extinct_at"] < 0]
    ser = [specials_series(z) for z in surv[c]]
    L2x = [l2(s["x"]) for s in ser]
    theta[c] = pct(L2x, 0.95)
    ev[c] = []
    for z, s in zip(surv[c], ser):
        d = {e: first_exceed(s[TAG[e]], theta[c]) for e in TAG}
        d["E_short"] = e_short(z)
        d["E_seed"] = e_seed(z)
        ev[c].append(d)
    result["conds"][c] = dict(planned=len(zs), survived=len(surv[c]), extinct=len(zs) - len(surv[c]),
                              theta=theta[c], L2_x_median=pct(L2x, 0.5),
                              L2_h_median=pct([l2(s["h"]) for s in ser], 0.5), L2_n_median=pct([l2(s["n"]) for s in ser], 0.5))

print("\n== 생존 · 중립 문턱 θ (= 생존 접시 L2_x 95 백분위)")
for c in CONDS:
    r = result["conds"][c]
    print(f"  {c:>6}: 생존 {r['survived']}/{r['planned']} · θ {r['theta']:.2%} · L2 중앙 x {r['L2_x_median']:.2%} · h {r['L2_h_median']:.2%} · n {r['L2_n_median']:.2%}")

# ---------------------------------------------------------------- 계측기 점검 (§2 개정판)
print("\n== 계측기 점검 — 자리만 · 재료에서 중립인 h 의 사건 발생률 (구간 하한 > 5% 면 Q1 보류)")
hold = {}
for c in CONDS:
    n = len(surv[c])
    if c in NEUTRAL_H:
        k = sum(d["E_h"] is not None for d in ev[c])
        lo, hi = wilson(k, n)
        hold[c] = lo > 0.05
        print(f"  {c:>6}: E_h {fmt_p(k, n)} → {'🚨 Q1 판정 보류' if hold[c] else '정상'}")
    else:
        hold[c] = False
result["q1_hold"] = hold

# ---------------------------------------------------------------- Q1
print("\n== Q1 사건 반복성 — 생존 접시 기준 발생 비율 · Wilson 95% (하한 ≥ 90% 매번 · 상한 ≤ 10% 거의 안 · 그 밖 우연)")
q1 = {}
for c in CONDS:
    n = len(surv[c])
    q1[c] = {}
    for e in EVENTS:
        times = [d[e] for d in ev[c] if d[e] is not None]
        k = len(times)
        lo, hi = wilson(k, n)
        control = (e == "E_h" and c in NEUTRAL_H)
        if control:
            verdict = "중립 대조(판정 없음)"
        elif hold[c]:
            verdict = "보류(계측기)"
        else:
            verdict = "매번 일어난다" if lo >= 0.90 else ("거의 안 일어난다" if hi <= 0.10 else "우연에 달렸다")
        q1[c][e] = dict(k=k, n=n, lo=lo, hi=hi, verdict=verdict,
                        t_median=pct(times, 0.5) if times else None, t_q1=pct(times, 0.25) if times else None, t_q3=pct(times, 0.75) if times else None)
        tt = f" · 시각 중앙 {q1[c][e]['t_median']:,.0f} [{q1[c][e]['t_q1']:,.0f}, {q1[c][e]['t_q3']:,.0f}]" if times else ""
        print(f"  {c:>6} {e:>7}: {fmt_p(k, n):>24} → {verdict}{tt}")
result["q1"] = q1

# ---------------------------------------------------------------- Q2
print("\n== Q2 끝 상태 수렴 — 조각(3-mer) 층 주 판정")
q2, simmat = {}, {}
for c in CONDS:
    zs = surv[c]
    n = len(zs)
    fin = [kmer_vec(z["final_counts"]) for z in zs]
    snap = [kmer_vec(z["snap"][1]) if z["snap"] else {} for z in zs]
    sm = {}
    for i, j in combinations(range(n), 2):
        v = cosine(fin[i], fin[j])
        if v == v:
            sm[(i, j)] = v
    simmat[c] = sm
    within = [v for v in (cosine(fin[i], snap[i]) for i in range(n)) if v == v]
    med_b, med_w, p5_w = pct(list(sm.values()), 0.5), pct(within, 0.5), pct(within, 0.05)
    if med_b != med_b or med_w != med_w:
        verdict = "판정 불가 — 유사도 값 없음"
    else:
        verdict = "같은 곳에 도착한다" if med_b >= med_w else ("갈라진다" if med_b < p5_w else "중간")
    doms = [z["final_counts"][0][0] for z in zs if z["final_counts"]]
    top_dom = max(set(doms), key=doms.count) if doms else None
    k_dom = doms.count(top_dom) if doms else 0
    lv = [letter_vec(z["final_counts"]) for z in zs]
    let_b = pct([v for v in (cosine(lv[i], lv[j]) for i, j in combinations(range(n), 2)) if v == v], 0.5)
    q2[c] = dict(n=n, between_median=med_b, within_median=med_w, within_p5=p5_w, verdict=verdict,
                 top_dominant=top_dom, top_dominant_share=fmt_p(k_dom, len(doms)), letter_between_median=let_b)
    print(f"  {c:>6}: 접시 사이 중앙 {med_b:.3f} · 접시 안(끝 대 2,000틱 전) 중앙 {med_w:.3f} · 5% {p5_w:.3f} → {verdict}")
    print(f"          코드 층: 가장 흔한 우세 코드 {top_dom} 공유 {q2[c]['top_dominant_share']} · 글자 층 접시 사이 중앙 {let_b:.3f}")
result["q2"] = q2


# ---------------------------------------------------------------- Q3
def q3_stats(evlist, use):
    sets = [{e: d[e] for e in use if d[e] is not None} for d in evlist]
    n = len(sets)
    jac, pairs = [], []
    for i, j in combinations(range(n), 2):
        A, B = set(sets[i]), set(sets[j])
        if A | B:
            jac.append(len(A & B) / len(A | B))
        common = sorted(A & B)
        if len(common) >= 3:
            pairs.append((i, j, common))

    def pair_taus(ss):
        out = {}
        for i, j, common in pairs:
            t = tau_b([ss[i][e] for e in common], [ss[j][e] for e in common])
            if t is not None:
                out[(i, j)] = t
        return out

    obs = pair_taus(sets)
    obs_mean = st.mean(obs.values()) if obs else float("nan")
    null = []
    if obs:
        for _ in range(N_NULL):
            sh = []
            for d in sets:
                keys, vals = list(d.keys()), list(d.values())
                RNG.shuffle(vals)
                sh.append(dict(zip(keys, vals)))
            t = pair_taus(sh)
            if t:
                null.append(st.mean(t.values()))
    thr = pct(null, 0.975)
    testable = bool(obs) and thr == thr
    return dict(n_pairs_tau=len(obs), jaccard_mean=st.mean(jac) if jac else float("nan"), tau_mean=obs_mean,
                null_p975=thr, testable=testable, repeats=(testable and obs_mean > thr)), obs


print("\n== Q3 순서 반복성 — 공통 사건 ≥ 3 인 쌍의 평균 Kendall τ-b 대 뒤섞기 귀무 97.5% (E_seed 뺀 판정이 갈리면 그쪽)")
q3, taumat = {}, {}
for c in CONDS:
    base = [e for e in EVENTS if not (e == "E_h" and c in NEUTRAL_H)]
    full, obs_full = q3_stats(ev[c], base)
    nosd, obs_nosd = q3_stats(ev[c], [e for e in base if e != "E_seed"])
    if not full["testable"] and not nosd["testable"]:
        final, src, taumat[c] = None, "—", {}
    elif full["testable"] and nosd["testable"] and full["repeats"] != nosd["repeats"]:
        final, src, taumat[c] = nosd, "E_seed 제외", obs_nosd
    elif full["testable"]:
        final, src, taumat[c] = full, "전체", obs_full
    else:
        final, src, taumat[c] = nosd, "E_seed 제외", obs_nosd
    verdict = "판정 불가 — 공통 사건 ≥ 3 인 접시 쌍 없음" if final is None else ("순서가 반복된다" if final["repeats"] else "순서 반복 없음")
    q3[c] = dict(with_seed=full, without_seed=nosd, verdict_from=src, verdict=verdict)
    for label, r in (("전체", full), ("E_seed 제외", nosd)):
        tau = f"{r['tau_mean']:+.3f}" if r["testable"] else "—"
        thr = f"{r['null_p975']:+.3f}" if r["testable"] else "—"
        print(f"  {c:>6} {label:>9}: τ 평균 {tau} (쌍 {r['n_pairs_tau']}) · 귀무 97.5% {thr} · Jaccard 평균 {r['jaccard_mean']:.3f}")
    print(f"  {c:>6} → {verdict} ({src} 기준)")
result["q3"] = q3

# ---------------------------------------------------------------- Q4
print(f"\n== Q4 조건 의존 — 접시 재추출 부트스트랩 {N_BOOT}회 · 차이 95% 구간이 0 을 벗어나면 '모자라는 것이 반복성을 바꾼다'")
eligible = [c for c in CONDS if len(surv[c]) >= MIN_SURV]
for c in CONDS:
    if c not in eligible:
        print(f"  {c}: 생존 {len(surv[c])} < {MIN_SURV} → 비교에서 뺀다")


def boot_dist(c):
    n = len(surv[c])
    dist = {"q2": [], "q3": []}
    for e in EVENTS:
        dist[e] = []
    for _ in range(N_BOOT):
        idx = [RNG.randrange(n) for _ in range(n)]
        for e in EVENTS:
            dist[e].append(sum(ev[c][i][e] is not None for i in idx) / n)
        sims, taus = [], []
        for p, q in combinations(range(n), 2):
            i, j = idx[p], idx[q]
            if i == j:
                continue
            key = (i, j) if i < j else (j, i)
            if key in simmat[c]:
                sims.append(simmat[c][key])
            if key in taumat[c]:
                taus.append(taumat[c][key])
        dist["q2"].append(pct(sims, 0.5))
        dist["q3"].append(st.mean(taus) if taus else float("nan"))
    return dist


boots = {c: boot_dist(c) for c in eligible}
q4 = {}
for a, b in combinations(eligible, 2):
    q4[f"{a}-{b}"] = {}
    stats = [e for e in EVENTS if not (e == "E_h" and (a in NEUTRAL_H or b in NEUTRAL_H))] + ["q2", "q3"]
    for sname in stats:
        label = {"q2": "Q2 접시 사이 유사도 중앙", "q3": "Q3 평균 τ"}.get(sname, f"Q1 {sname} 발생 비율")
        diffs = [x - y for x, y in zip(boots[a][sname], boots[b][sname]) if x == x and y == y]
        if not diffs:
            q4[f"{a}-{b}"][sname] = dict(lo=None, hi=None, differs=None)
            print(f"  {a:>6} − {b:<6} {label:>24}: 판정 불가 — 비교할 값 없음")
            continue
        lo, hi = pct(diffs, 0.025), pct(diffs, 0.975)
        differs = (lo > 0 or hi < 0)
        q4[f"{a}-{b}"][sname] = dict(lo=lo, hi=hi, differs=differs)
        print(f"  {a:>6} − {b:<6} {label:>24}: [{lo:+.3f}, {hi:+.3f}] → {'다르다' if differs else '차이 검출 안 됨'}")
result["q4"] = q4

json.dump(result, open(os.path.join(DIR, "..", OUT), "w"), ensure_ascii=False, indent=1, default=str)
print(f"\n저장: {OUT}")
