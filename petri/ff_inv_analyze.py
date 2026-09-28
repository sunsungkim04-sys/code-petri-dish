# -*- coding: utf-8 -*-
"""해부 7 ② 인과 판정 — 재료 먼저 규칙에서 씨앗의 오류 불이익이 사라지나. 사전등록 해부7-사전등록-2026-09-16.md §2 · §4.
결과 전에 썼다.

  원래 규칙 = 해부 6C 기록(inv6/) · 재료 먼저 = inv7ff/ — **같은 20 배경 · 같은 μ(0 · 0.1 · 0.3 · 0.5 · 1%)** 로 짝짓는다.
  Δ = s(팔) − WT 팔 11개 평균(dms_analyze.py 와 같은 s) · 기울기 = Δ 를 μ 에 최소제곱 · 배경 재추출 2,000회 95%

  F4 핵심 : 씨앗 기울기 차 = 재료먼저 − 원래. 줄어든 몫 R = 1 − 기울기_재료먼저/기울기_원래
            차 하한 > 0 이고 R ≥ 1/2 → '씨앗의 오류 불이익이 대부분 사라진다 — 가설이 맞다'
            차 하한 > 0 이고 R < 1/2 → '일부만 준다'
            차 하한 ≤ 0 → '줄지 않는다 — 가설 기각'
  F5 음성 : μ 0 의 씨앗 Δ 차(재료먼저 − 원래) 구간이 0 을 품는다 → 통과 · 못 품으면 🚩 '오류가 없어도 규칙 바꿈이 차이를 만든다'
  F6 음성 : 재료먼저 racldx 기울기 구간이 0 을 품는다 → 통과 · 못 품으면 🚩
  서술    : rsacld 기울기 두 규칙 · Δ 표 · 계통 순도
  점검    : 재료먼저 파일의 find_first 표시 · 배경 체크섬 · 복제 충실도 · 재료 보존 · 팔 14 · μ 일치 · 완결성
실행: python3 ff_inv_analyze.py [petri 폴더]  →  _RESULT_ff_inv.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOTS = {"orig": os.environ.get("PETRI_ROOT_ORIG", "inv6"), "ff": os.environ.get("PETRI_ROOT", "inv7ff")}
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
RNG = random.Random(20261007)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


F, bad = {}, []
for rule, root in ROOTS.items():
    for mu in MUS:
        for f in sorted(glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f))
            F[(rule, mu, z["seed"])] = z
            if abs(z["mu_assay"] - mu) > 1e-12: bad.append("μ %s %s %d" % (rule, mu, z["seed"]))
            if bool(z.get("find_first")) != (rule == "ff"): bad.append("규칙 표시 %s %s %d" % (rule, mu, z["seed"]))
            if not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (rule, mu, z["seed"]))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (rule, mu, z["seed"]))
            if len(z["arms"]) != 14: bad.append("팔 수 %s %s %d" % (rule, mu, z["seed"]))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (rule, mu, z["seed"]))
bgs = sorted({k[2] for k in F if k[0] == "ff"})
complete = len(bgs) == NBG and all((r, mu, s) in F for r in ROOTS for mu in MUS for s in bgs)
print("== 해부 7 ② 인과 — 재료 먼저 규칙의 침입 사다리 — 판정")
n_used = sum(1 for k in F if k[2] in bgs)
print("  파일 %d/%d(두 규칙 × μ %d × 배경 %d) · 점검 문제 %d" % (n_used, 2 * len(MUS) * NBG, len(MUS), len(bgs), len(bad)))
if not complete or bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)

IDX = [[RNG.randrange(len(bgs)) for _ in bgs] for _ in range(N_BOOT)]
delta = {r: {c: {mu: [] for mu in MUS} for c in CODES} for r in ROOTS}
pure = {r: {c: {mu: [] for mu in MUS} for c in CODES + ["ref"]} for r in ROOTS}
for r in ROOTS:
    for mu in MUS:
        for s in bgs:
            z = F[(r, mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                last = a["series"][-1]
                sh = last[3] / last[2] if last[2] else float("nan")
                if a["kind"] == "mut":
                    delta[r][a["key"]][mu].append(s_of(a) - base)
                    pure[r][a["key"]][mu].append(sh)
                elif a["kind"] == "ref":
                    pure[r]["ref"][mu].append(sh)


def slope(ys):
    xm = sum(MUS) / len(MUS)
    ym = sum(ys) / len(ys)
    return sum((x - xm) * (y - ym) for x, y in zip(MUS, ys)) / sum((x - xm) ** 2 for x in MUS)


def slope_boot(r, c):
    pt = slope([sum(delta[r][c][mu]) / len(bgs) for mu in MUS])
    bs = [slope([sum(delta[r][c][mu][i] for i in ix) / len(ix) for mu in MUS]) for ix in IDX]
    return pt, bs


def ci(bs):
    return pct(bs, 0.025), pct(bs, 0.975)


out = {"delta": {}, "slopes": {}, "F4": {}, "F5": {}, "F6": {}, "purity": {}}
print("\n   %6s | %-30s | %-30s" % ("μ", "씨앗 Δ 원래 → 재료먼저", "rsacld Δ 원래 → 재료먼저"))
for mu in MUS:
    row = []
    for c in ("rascld", "rsacld"):
        a = sum(delta["orig"][c][mu]) / len(bgs); b = sum(delta["ff"][c][mu]) / len(bgs)
        out["delta"]["%s|%s" % (c, mu)] = {"orig": a, "ff": b}
        row.append("%+.3f → %+.3f" % (a, b))
    x = sum(delta["ff"]["racldx"][mu]) / len(bgs)
    out["delta"]["racldx|%s" % mu] = {"orig": sum(delta["orig"]["racldx"][mu]) / len(bgs), "ff": x}
    print("   %5.2f%% | %-30s | %-30s | racldx 재료먼저 %+.3f" % (mu * 100, row[0], row[1], x))

sl = {}
for r in ROOTS:
    for c in CODES:
        sl[(r, c)] = slope_boot(r, c)
        out["slopes"]["%s|%s" % (r, c)] = [sl[(r, c)][0], *ci(sl[(r, c)][1])]
print("\n   기울기(Δ / μ):")
for c in CODES:
    print("     %-7s 원래 %+.1f [%+.1f, %+.1f] · 재료먼저 %+.1f [%+.1f, %+.1f]" % (c, *out["slopes"]["orig|%s" % c], *out["slopes"]["ff|%s" % c]))

po, pf = sl[("orig", "rascld")], sl[("ff", "rascld")]
dbs = [b - a for a, b in zip(po[1], pf[1])]
dm = pf[0] - po[0]
dlo, dhi = ci(dbs)
R = 1 - pf[0] / po[0] if po[0] != 0 else float("nan")
Rbs = [1 - b / a for a, b in zip(po[1], pf[1]) if a != 0]
if dlo > 0 and R >= 0.5:
    v4 = "씨앗의 오류 불이익이 대부분 사라진다 — 가설이 맞다"
elif dlo > 0:
    v4 = "일부만 준다"
else:
    v4 = "줄지 않는다 — 가설 기각"
out["F4"] = {"diff": [dm, dlo, dhi], "R": R, "R_ci": list(ci(Rbs)), "verdict": v4}
print("\n-- F4 씨앗 기울기 차(재료먼저 − 원래) %+.1f [%+.1f, %+.1f] · 줄어든 몫 R %.2f [%.2f, %.2f]" % (dm, dlo, dhi, R, *ci(Rbs)))
print("   → **%s**" % v4)

d0 = [b - a for a, b in zip(delta["orig"]["rascld"][0.0], delta["ff"]["rascld"][0.0])]
m0 = sum(d0) / len(d0)
bs0 = [sum(d0[i] for i in ix) / len(ix) for ix in IDX]
f5 = ci(bs0)[0] <= 0 <= ci(bs0)[1]
out["F5"] = {"diff": [m0, *ci(bs0)], "passes": f5}
print("-- F5 μ 0 씨앗 Δ 차 %+.3f [%+.3f, %+.3f] → %s" % (m0, *ci(bs0), "통과 — 오류가 없으면 규칙 바꿈이 차이를 안 만든다" if f5 else "🚩 오류가 없어도 규칙 바꿈이 차이를 만든다"))
xs = out["slopes"]["ff|racldx"]
f6 = xs[1] <= 0 <= xs[2]
out["F6"] = {"slope": xs, "passes": f6}
print("-- F6 재료먼저 racldx 기울기 %+.1f [%+.1f, %+.1f] → %s" % (*xs, "통과" if f6 else "🚩 끼어든 x 도 깎인다"))

print("\n-- 서술: 끝 기록의 계통 순도 (원래 / 재료먼저)")
for mu in MUS:
    cells = []
    for c in ["ref"] + CODES:
        a = [x for x in pure["orig"][c][mu] if not math.isnan(x)]; b = [x for x in pure["ff"][c][mu] if not math.isnan(x)]
        cells.append("%s %.0f/%.0f%%" % (c, 100 * sum(a) / len(a), 100 * sum(b) / len(b)))
        out["purity"]["%s|%s" % (c, mu)] = [sum(a) / len(a), sum(b) / len(b)]
    print("   %5.2f%% | %s" % (mu * 100, " · ".join(cells)))

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_ff_inv.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_ff_inv.json"))
