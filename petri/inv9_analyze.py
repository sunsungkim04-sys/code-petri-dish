# -*- coding: utf-8 -*-
"""해부 9 판정 — 기억하는 주사위 세계(A) · 재료 규칙을 끈 세계(B)의 침입 μ 사다리. 사전등록 해부9-사전등록-2026-09-17.md §2~§4.
결과 전에 썼다(파일럿은 노트 §5 — B 는 파일럿에서 설계를 바꿨다: space 공동체에 끼운 racld 계통이 전부 죽어, 같은 mat 공동체에서 팔 동안만 재료 규칙을 끄는 것으로).

  같은 mat 배경 20(시드 3~22) · 같은 μ(0 · .001 · .003 · .005 · .01) · 같은 팔(ref + neutral 10 + rascld · rsacld · racldx) — 네 세계가 **짝**이다.
    orig  inv6/    원래 규칙(해부 6C 기록)
    ff    inv7ff/  재료 먼저(해부 7 ② 기록) — 주사위 몫 R → 1 · 거름 F 는 남음
    mem   inv9mem/ 기억하는 주사위(sim v0.3.4 rememberDie · 팔 동안만) — 주사위 몫 R → 1 · 기억한 오류는 글자를 기다려 실현되므로 거름도 거의 없음(파일럿 F 0.94 · 부풀림 0.98)
    nm    inv9nm/  재료 규칙 끔(dms.js --no-mat 1 · 팔 동안만: w.mat = null) — 헛손질 0 · 거름 0 · 재료 생태 0
  Δ = s(팔) − WT 팔 11개 평균 · 기울기 = μ 다섯 점 최소제곱 · 배경 재추출 2,000회 95%(네 세계 같은 인덱스 → 짝)

  M1 (A 핵심): 씨앗 기울기 차(기억 − 원래) 하한 > 0 이고 R_mem = 1 − b_mem/b_orig ≥ 1/2
              → '주사위를 글자당 한 번으로 만드는 것만으로 불이익 대부분이 준다(재료 생태는 그대로)' · 하한 > 0 · R < 1/2 '일부' · 하한 ≤ 0 '안 준다'
  M2 (A 대 재료 먼저): 씨앗 기울기 차(기억 − 재료 먼저) 구간이 0 을 품는다 → '두 조작이 구별 안 됨' · 하한 > 0 '기억 쪽이 덜 깎임' · 상한 < 0 '재료 먼저 쪽이 덜 깎임'
  M3 (버그 없음 점검): μ 0 씨앗 Δ 차(기억 − 원래) 구간이 0 을 품고, 기억 세계의 racldx 기울기 구간이 0 을 품는다
  S1 (B 핵심): 재료 끔 세계의 씨앗 기울기 구간이 0 을 품는다 → '헛손질 없는 세계에서는 오류 불이익이 없다'
              못 품되 |b_nm| 의 구간 끝 < |b_orig|/3 → '작다(원래의 1/3 미만)' · 그 밖 '있다 — 헛손질 밖의 이유가 있다'
  S2 (서술): 재료 끔의 μ 0 Δ(씨앗 · rsacld · racldx — 재료 생태가 없을 때의 고리 이득) · rsacld 기울기 · 끼운 계통이 사라진 팔 수 · 특이 몫 대조: (b_orig − b_nm) 대 (b_orig − b_mem)
  점검: 파일 · μ 일치 · 규칙 표시(remember_die · no_mat · find_first) · 배경 체크섬 · 충실도 · 팔 14 · 재료 보존(nm 은 참으로 둔다)
실행: python3 inv9_analyze.py [petri 폴더]  →  _RESULT_inv9.json
마른 실행: PETRI_ROOT_MEM=stage9/pilot9/mem PETRI_ROOT_NM=stage9/pilot9/nm PETRI_NBG=1 PETRI_MUS=0,0.01 python3 inv9_analyze.py
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOTS = {"orig": os.environ.get("PETRI_ROOT_ORIG", "inv6"), "ff": os.environ.get("PETRI_ROOT_FF", "inv7ff"),
         "mem": os.environ.get("PETRI_ROOT_MEM", "inv9mem"), "nm": os.environ.get("PETRI_ROOT_NM", "inv9nm")}
FLAG = {"orig": None, "ff": "find_first", "mem": "remember_die", "nm": "no_mat"}
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
RNG = random.Random(20260921)


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
            for r2, fl in FLAG.items():
                if fl and bool(z.get(fl)) != (rule == r2): bad.append("규칙 표시 %s(%s) %s %d" % (rule, fl, mu, z["seed"]))
            if not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (rule, mu, z["seed"]))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (rule, mu, z["seed"]))
            if len(z["arms"]) != 14: bad.append("팔 수 %s %s %d" % (rule, mu, z["seed"]))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (rule, mu, z["seed"]))
bgs = sorted({k[2] for k in F if k[0] == "mem"})
complete = len(bgs) == NBG and all((r, mu, s) in F for r in ROOTS for mu in MUS for s in bgs)
print("== 해부 9 — 기억하는 주사위 세계 · 재료 끔 세계의 침입 사다리 — 판정")
print("  배경 %d · μ %d · 파일 %d/%d · 점검 문제 %d" % (len(bgs), len(MUS), sum(1 for k in F if k[2] in bgs), 4 * len(MUS) * NBG, len(bad)))
if not complete or bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)

delta = {r: {c: {mu: [] for mu in MUS} for c in CODES} for r in ROOTS}
lost = {r: 0 for r in ROOTS}
for r in ROOTS:
    for mu in MUS:
        for s in bgs:
            z = F[(r, mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                if a["kind"] == "mut":
                    delta[r][a["key"]][mu].append(s_of(a) - base)
                if a["series"][-1][2] == 0:
                    lost[r] += 1
IDX = [[RNG.randrange(NBG) for _ in range(NBG)] for _ in range(N_BOOT)]


def mean(r, c, mu, ix=None):
    v = delta[r][c][mu]
    return sum(v[i] for i in ix) / len(ix) if ix is not None else sum(v) / len(v)


def slope(ys):
    mx = sum(MUS) / len(MUS)
    my = sum(ys) / len(ys)
    return sum((m - mx) * (y - my) for m, y in zip(MUS, ys)) / sum((m - mx) ** 2 for m in MUS)


def sl(r, c, ix=None):
    return slope([mean(r, c, mu, ix) for mu in MUS])


def ci(v):
    return pct(v, .025), pct(v, .975)


out = {"delta": {}, "slopes": {}, "lost": lost, "bgs": bgs}
print("\n-- 씨앗 Δ 표 (배경 평균 [95%])")
print("   %6s | %-24s | %-24s | %-24s | %-24s" % ("μ", "원래", "재료 먼저", "**기억 주사위**", "**재료 끔**"))
for mu in MUS:
    cells = []
    for r in ROOTS:
        m = mean(r, "rascld", mu)
        lo, hi = ci([mean(r, "rascld", mu, ix) for ix in IDX])
        cells.append("%+.3f [%+.3f, %+.3f]" % (m, lo, hi))
        out["delta"]["%s|rascld|%g" % (r, mu)] = [m, lo, hi]
    print("   %6g | %s" % (mu, " | ".join(cells)))
for r in ROOTS:
    for c in CODES:
        pt = sl(r, c)
        lo, hi = ci([sl(r, c, ix) for ix in IDX])
        out["slopes"]["%s|%s" % (r, c)] = [pt, lo, hi]
print("\n-- 기울기 (Δ / μ · 다섯 점)")
for c in CODES:
    print("   %-7s " % c + " · ".join("%s %+.1f [%+.1f, %+.1f]" % (r, *out["slopes"]["%s|%s" % (r, c)]) for r in ROOTS))

# M1
d_mo = [sl("mem", "rascld", ix) - sl("orig", "rascld", ix) for ix in IDX]
R_mem = [1 - sl("mem", "rascld", ix) / sl("orig", "rascld", ix) for ix in IDX]
pt_d = sl("mem", "rascld") - sl("orig", "rascld")
pt_R = 1 - sl("mem", "rascld") / sl("orig", "rascld")
lo, hi = ci(d_mo)
rlo, rhi = ci(R_mem)
if lo > 0 and pt_R >= 0.5:
    v1 = "주사위를 글자당 한 번으로 만드는 것만으로 불이익 대부분이 준다(재료 생태는 그대로)"
elif lo > 0:
    v1 = "일부만 준다"
else:
    v1 = "줄지 않는다"
print("\n-- M1 기억 세계: 씨앗 기울기 차(기억 − 원래) %+.1f [%+.1f, %+.1f] · R_mem %.2f [%.2f, %.2f] → **%s**" % (pt_d, lo, hi, pt_R, rlo, rhi, v1))
out["M1"] = {"diff": [pt_d, lo, hi], "R": [pt_R, rlo, rhi], "verdict": v1}
# M2
d_mf = [sl("mem", "rascld", ix) - sl("ff", "rascld", ix) for ix in IDX]
pt2 = sl("mem", "rascld") - sl("ff", "rascld")
lo2, hi2 = ci(d_mf)
v2 = "두 조작이 구별 안 됨" if lo2 <= 0 <= hi2 else ("기억 쪽이 덜 깎임" if lo2 > 0 else "재료 먼저 쪽이 덜 깎임")
print("-- M2 기억 대 재료 먼저: 씨앗 기울기 차 %+.1f [%+.1f, %+.1f] → **%s**" % (pt2, lo2, hi2, v2))
out["M2"] = {"diff": [pt2, lo2, hi2], "verdict": v2}
# M3
d0 = [mean("mem", "rascld", 0.0, ix) - mean("orig", "rascld", 0.0, ix) for ix in IDX]
lo3, hi3 = ci(d0)
xs = out["slopes"]["mem|racldx"]
ok3 = (lo3 <= 0 <= hi3) and (xs[1] <= 0 <= xs[2])
print("-- M3 버그 없음 점검: μ 0 씨앗 Δ 차 %+.3f [%+.3f, %+.3f] · 기억 racldx 기울기 %+.1f [%+.1f, %+.1f] → %s"
      % (mean("mem", "rascld", 0.0) - mean("orig", "rascld", 0.0), lo3, hi3, *xs, "✅ 통과" if ok3 else "🚩 실패"))
out["M3"] = {"mu0_diff": [mean("mem", "rascld", 0.0) - mean("orig", "rascld", 0.0), lo3, hi3], "racldx": xs, "passes": ok3}
# S1
nm = out["slopes"]["nm|rascld"]
bo = abs(out["slopes"]["orig|rascld"][0])
if nm[1] <= 0 <= nm[2]:
    v4 = "헛손질 없는 세계에서는 오류 불이익이 없다"
elif max(abs(nm[1]), abs(nm[2])) < bo / 3:
    v4 = "작다 — 원래의 1/3 미만"
else:
    v4 = "있다 — 헛손질 밖의 이유가 있다"
print("-- S1 재료 끔 세계: 씨앗 기울기 %+.1f [%+.1f, %+.1f] (원래 %+.1f) → **%s**" % (*nm, out["slopes"]["orig|rascld"][0], v4))
out["S1"] = {"slope": nm, "verdict": v4}
d_on = [sl("nm", "rascld", ix) - sl("orig", "rascld", ix) for ix in IDX]
d_om = [sl("mem", "rascld", ix) - sl("orig", "rascld", ix) for ix in IDX]
print("-- S2 서술: 재료 끔 μ 0 Δ 씨앗 %+.3f · rsacld %+.3f · racldx %+.3f · rsacld 기울기 %+.1f [%+.1f, %+.1f] · 끼운 계통 사라진 팔 %s"
      % (mean("nm", "rascld", 0.0), mean("nm", "rsacld", 0.0), mean("nm", "racldx", 0.0), *out["slopes"]["nm|rsacld"], lost))
print("          기울기 차(세계 − 원래 · 양수 = 덜 깎임): 재료 끔 %+.1f [%+.1f, %+.1f] 대 기억 %+.1f [%+.1f, %+.1f]"
      % (sl("nm", "rascld") - sl("orig", "rascld"), *ci(d_on), sl("mem", "rascld") - sl("orig", "rascld"), *ci(d_om)))
out["S2"] = {"nm_minus_orig": [sl("nm", "rascld") - sl("orig", "rascld"), *ci(d_on)], "mem_minus_orig": [sl("mem", "rascld") - sl("orig", "rascld"), *ci(d_om)]}
json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_inv9.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_inv9.json")
