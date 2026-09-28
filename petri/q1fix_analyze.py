"""사건 문턱 고정 판정 — 사전등록 사건문턱-고정-사전등록-2026-09-15.md. 보정 접시 결과를 보기 전에 썼다.

  θ_c = 보정 전용 접시(시드 301~500) 생존 접시 L2_x 의 95 백분위 — 판정 접시(1~100 · 101~200)는 θ 에 쓰지 않는다
  사건 = 본실험 §2 개정판의 L2 사건 넷(E_j · E_noc · E_e · E_h). E_short · E_seed 는 θ 를 안 써서 다시 판정하지 않는다
  주 판정 = 판정 접시 1~200 을 합친 발생 비율 · Wilson 95% (하한 ≥ 0.90 매번 · 상한 ≤ 0.10 거의 안 · 그 밖 우연)
    자리만 · 재료의 E_h 는 중립 대조 — 합친 구간 하한 > 0.05 면 그 조건 Q1 보류(본실험 규칙 그대로)
  재현 = 1~100 과 101~200 을 **같은 θ_c** 로 따로 판정해 범주가 같나
  안정도 = 보정 접시 생존 목록 재추출 2,000회 θ* 로 합친 판정을 다시 내 주 판정과 같은 범주인 비율 — 0.90 미만이면 '문턱 민감'
  자기검산 = 이 스크립트의 L2 · 백분위로 본실험 · 되감기 2 의 옛 θ 를 다시 재서 결과 JSON 과 같아야 한다 — 다르면 중단
  값이 없으면 판정 불가(분모 0)
실행: python3 q1fix_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_q1fix.json
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
CONDS = ["space", "energy", "mat"]
EVENTS = [("E_j", "j"), ("E_noc", "noc"), ("E_e", "e"), ("E_h", "h")]
NEUTRAL_H = {"space", "mat"}
SETS = {
    "main": ("main", range(1, 101), {"0.3.0"}, "_RESULT_main.json"),
    "rewind2": ("rewind2", range(101, 201), {"0.3.1"}, "_RESULT_rewind2_q1.json"),
    "calib": ("calib", range(301, 501), {"0.3.1"}, None),
}
N_BOOT = 2000
RNG = random.Random(20260919)


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


def category(k, n):
    if n == 0:
        return "판정 불가"
    lo, hi = wilson(k, n)
    return "매번 일어난다" if lo >= 0.90 else ("거의 안 일어난다" if hi <= 0.10 else "우연에 달렸다")


def l2(z, tag):
    pop_at = {r[0]: r[1] for r in z["samples"]}
    col = z["special_header"].index(tag) + 1
    s = []
    for row in z["specials"]:
        pop = pop_at.get(row[0], 0)
        s.append((row[col][1] / pop) if (row[col] and pop) else 0.0)
    return max((min(s[i], s[i + 1]) for i in range(len(s) - 1)), default=0.0)


L2 = {}
problems = []
for name, (folder, seeds, versions, _) in SETS.items():
    L2[name] = {c: [] for c in CONDS}
    for c in CONDS:
        for s in seeds:
            f = os.path.join(BASE, folder, f"{c}_d8_mu0p01_s{s:05d}.json")
            if not os.path.exists(f):
                problems.append(f"없음 {folder}/{os.path.basename(f)}")
                continue
            z = json.load(open(f))
            if z["sim_version"] not in versions or z["ticks_planned"] != 50000 or z["opts"]["mu"] != 0.01 or z["opts"]["density"] != 8:
                problems.append(f"설정 다름 {folder}/{os.path.basename(f)}")
            if z["conservation_violations"]:
                problems.append(f"재료 보존 위반 {folder}/{os.path.basename(f)}")
            if z["extinct_at"] < 0:
                L2[name][c].append({tag: l2(z, tag) for tag in ("x", "j", "noc", "e", "h")})
print("자료: " + " · ".join(f"{n} " + "/".join(str(len(L2[n][c])) for c in CONDS) + " 생존" for n in SETS))
if problems:
    print(f"🚨 완결성 가드 {len(problems)}건 — 판정하지 않는다")
    for p in problems[:20]:
        print("   ", p)
    sys.exit(1)

print("\n== 자기검산 — 옛 θ 를 이 스크립트로 다시 재기")
for name in ("main", "rewind2"):
    R = json.load(open(os.path.join(BASE, SETS[name][3])))
    for c in CONDS:
        mine = pct([d["x"] for d in L2[name][c]], .95)
        theirs = R["conds"][c]["theta"]
        ok = abs(mine - theirs) < 1e-12
        print(f"  {name:>7} {c:>6}: 다시 잰 θ {mine:.6%} · 결과 JSON {theirs:.6%} → {'일치' if ok else '🚨 불일치'}")
        if not ok:
            print("🚨 L2 · θ 구현이 판정 스크립트와 다르다 — 판정하지 않는다")
            sys.exit(1)

out = {"theta_c": {}, "conds": {}}
print("\n== 보정 접시 θ_c (시드 301~500 생존 접시 L2_x 95 백분위) · 옛 θ 와 나란히")
for c in CONDS:
    th = pct([d["x"] for d in L2["calib"][c]], .95)
    old = [pct([d["x"] for d in L2[n][c]], .95) for n in ("main", "rewind2")]
    out["theta_c"][c] = th
    print(f"  {c:>6}: θ_c {th:.2%} (보정 접시 {len(L2['calib'][c])}) · 옛 θ 본실험 {old[0]:.2%} · 되감기 2 {old[1]:.2%}")

print("\n== 사건 판정 — 같은 θ_c · 합친 1~200(주) · 1~100 · 101~200 · 안정도")
for c in CONDS:
    th = out["theta_c"][c]
    pooled = L2["main"][c] + L2["rewind2"][c]
    calib_x = [d["x"] for d in L2["calib"][c]]
    boot_th = [pct([calib_x[RNG.randrange(len(calib_x))] for _ in calib_x], .95) for _ in range(N_BOOT)]
    out["conds"][c] = {}
    hold = False
    for ev, tag in EVENTS:
        k = sum(d[tag] > th for d in pooled)
        n = len(pooled)
        lo, hi = wilson(k, n)
        k1 = sum(d[tag] > th for d in L2["main"][c])
        k2 = sum(d[tag] > th for d in L2["rewind2"][c])
        n1, n2 = len(L2["main"][c]), len(L2["rewind2"][c])
        if ev == "E_h" and c in NEUTRAL_H:
            status = "정상" if not (lo > 0.05) else "🚨 발동 — 이 조건 Q1 보류"
            hold = hold or lo > 0.05
            print(f"  {c:>6} {ev:>5}: 중립 대조 {k}/{n} = {k / n:.0%} [{lo:.0%}, {hi:.0%}] → 계측기 점검 {status} · 1~100 {k1}/{n1} · 101~200 {k2}/{n2}")
            out["conds"][c][ev] = {"k": k, "n": n, "lo": lo, "hi": hi, "check": status, "k1": k1, "k2": k2}
            continue
        cat = category(k, n)
        cat1, cat2 = category(k1, n1), category(k2, n2)
        same = sum(category(sum(d[tag] > t for d in pooled), n) == cat for t in boot_th) / N_BOOT
        out["conds"][c][ev] = {"k": k, "n": n, "lo": lo, "hi": hi, "cat": cat, "k1": k1, "n1": n1, "cat1": cat1, "k2": k2, "n2": n2,
                               "cat2": cat2, "replicated": cat1 == cat2, "stability": same}
        print(f"  {c:>6} {ev:>5}: 합친 {k}/{n} = {k / n:.0%} [{lo:.0%}, {hi:.0%}] → {cat} · 안정도 {same:.0%}{' ⚠️ 문턱 민감' if same < 0.90 else ''}"
              f" · 1~100 {k1}/{n1} {cat1} · 101~200 {k2}/{n2} {cat2} → {'재현' if cat1 == cat2 else '재현 안 됨'}")
    if hold:
        print(f"  {c:>6} → 🚨 계측기 점검 발동 — 이 조건의 위 판정은 보류")
        out["conds"][c]["hold"] = True

json.dump(out, open(os.path.join(BASE, "_RESULT_q1fix.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_q1fix.json")
