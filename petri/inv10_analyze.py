# -*- coding: utf-8 -*-
"""해부 10 판정 C · G — 우주선 돌연변이 끔(C) · 둘째 시조 세계(G)의 침입 μ 사다리. 사전등록 해부10-사전등록-2026-09-18.md §2 · §5.
결과 전에 썼다(파일럿은 노트 §6).

  C  같은 mat 배경 20(시드 3~22) · μ 5 · 팔 14 — 짝: orig(inv6/) · ff(inv7ff/) · c0(inv10c0/ 원래 규칙 + 우주선 끔) · ffc0(inv10ffc0/ 재료 먼저 + 우주선 끔)
  G  --ancestor racld 로 50,000틱 키운 배경(시드 1101~1120) · 원래 규칙 · μ 5 · 팔 14 (inv10anc/) — 본실험 기록이 없어 배경 체크섬은 'no-main'
     배경 규칙: 키운 끝 우세가 racld 가 아니거나 멸종한 배경은 뺀다(수 보고)
  Δ = s(팔) − WT 팔 11개 평균 · 기울기 = μ 다섯 점 최소제곱 · 배경 재추출 2,000회 95%(C 의 네 세계는 같은 인덱스 → 짝)

  C1 (재료 먼저 + 우주선 끔): 씨앗 기울기 구간이 0 을 품는다 → '재료 먼저 세계에 남은 불이익(−10)은 우주선 몫 — 전부 사라진다'
      못 품되 차(ffc0 − ff) 하한 > 0 → '일부는 우주선 몫' · 그 밖 '우주선 몫 아님'
  C2 (원래 + 우주선 끔): 씨앗 기울기 차(c0 − orig) 구간이 0 을 품는다 → '원래 세계의 불이익에 우주선 몫은 검출되지 않는다' · 하한 > 0 '일부 우주선 몫' · 상한 < 0 '우주선을 끄면 더 깎인다'
  C3 (버그 없음 점검): μ 0 씨앗 Δ 차(c0 − orig · ffc0 − ff) 구간이 0 을 품는다 · racldx 기울기는 서술(해부 9 M3 에서 세계에 따라 다름을 봤다)
  G1 (둘째 시조 재현): 씨앗 μ 0 Δ 하한 > 0 이고 씨앗 기울기 상한 < 0 → 'F4 의 두 전제(μ 0 씨앗 우위 · μ 불이익)가 다른 시조 세계에서 재현'
      하나만 → '절반 재현' · 둘 다 아님 → '재현 안 됨'
  G2 (서술): 기울기 비 b_anc/b_orig · 뺀 배경 수 · 끝 우세 코드
  점검: 파일 · μ · 규칙 표시(cosmic_off · find_first · ancestor) · C 는 배경 체크섬 · 충실도 · 팔 14 · 재료 보존 · G 는 충실도 · 팔 14 · 재료 보존 · main 없음
실행: python3 inv10_analyze.py [petri 폴더]  →  _RESULT_inv10.json
마른 실행: PETRI_ROOT_C0=stage10/pilot10/c0 PETRI_ROOT_FFC0=stage10/pilot10/ffc0 PETRI_ROOT_ANC=stage10/pilot10/anc PETRI_NBG=1 PETRI_MUS=0.01 python3 inv10_analyze.py
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOTS = {"orig": os.environ.get("PETRI_ROOT_ORIG", "inv6"), "ff": os.environ.get("PETRI_ROOT_FF", "inv7ff"),
         "c0": os.environ.get("PETRI_ROOT_C0", "inv10c0"), "ffc0": os.environ.get("PETRI_ROOT_FFC0", "inv10ffc0"),
         "anc": os.environ.get("PETRI_ROOT_ANC", "inv10anc")}
EXPECT = {"orig": (False, False, None), "ff": (True, False, None), "c0": (False, True, None), "ffc0": (True, True, None), "anc": (False, False, "racld")}
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
RNG = random.Random(20260922)


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


F, bad, dropped_anc = {}, [], {}
for rule, root in ROOTS.items():
    ff_e, c0_e, anc_e = EXPECT[rule]
    for mu in MUS:
        for f in sorted(glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f))
            s = z["seed"]
            if abs(z["mu_assay"] - mu) > 1e-12: bad.append("μ %s %s %d" % (rule, mu, s))
            if bool(z.get("find_first")) != ff_e or bool(z.get("cosmic_off")) != c0_e or z.get("ancestor") != anc_e:
                bad.append("규칙 표시 %s %s %d" % (rule, mu, s))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (rule, mu, s))
            if rule == "anc":
                if z["main_checksum"] is not None: bad.append("anc 에 main 이 있다 %s %d" % (mu, s))
                if z["grow_extinct"] >= 0 or not z["grow_top"] or z["grow_top"][0][0] != "racld" or len(z["arms"]) != 14:
                    dropped_anc[s] = "멸종" if z["grow_extinct"] >= 0 else ("끝 우세 %s" % (z["grow_top"][0][0] if z["grow_top"] else "-"))
                    continue
            else:
                if not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (rule, mu, s))
                if len(z["arms"]) != 14: bad.append("팔 수 %s %s %d" % (rule, mu, s))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (rule, mu, s))
            F[(rule, mu, s)] = z
bg_c = sorted({k[2] for k in F if k[0] == "c0"})
bg_a = sorted({k[2] for k in F if k[0] == "anc"})
bg_a = [s for s in bg_a if all(("anc", mu, s) in F for mu in MUS)]
complete = len(bg_c) == NBG and all((r, mu, s) in F for r in ("orig", "ff", "c0", "ffc0") for mu in MUS for s in bg_c)
print("== 해부 10 C · G — 우주선 끔 · 둘째 시조 침입 사다리 — 판정")
print("  C 배경 %d · G 배경 %d(뺀 것 %d: %s) · μ %d · 점검 문제 %d" % (len(bg_c), len(bg_a), len(dropped_anc), ", ".join("%d %s" % kv for kv in sorted(dropped_anc.items())) or "없음", len(MUS), len(bad)))
if not complete or bad or len(bg_a) < max(1, NBG // 2):   # 사전등록: 남은 배경이 NBG/2(=10) 미만이면 판정하지 않는다
    print("🚨 점검 실패 — 판정하지 않는다 (C 완결 %s · G 배경 %d)" % (complete, len(bg_a)))
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)

BG = {"orig": bg_c, "ff": bg_c, "c0": bg_c, "ffc0": bg_c, "anc": bg_a}
delta = {r: {c: {mu: [] for mu in MUS} for c in CODES} for r in ROOTS}
for r in ROOTS:
    for mu in MUS:
        for s in BG[r]:
            z = F[(r, mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                if a["kind"] == "mut":
                    delta[r][a["key"]][mu].append(s_of(a) - base)
IDX = [[RNG.randrange(len(bg_c)) for _ in bg_c] for _ in range(N_BOOT)]
IDXA = [[RNG.randrange(len(bg_a)) for _ in bg_a] for _ in range(N_BOOT)]


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


IX = {r: (IDXA if r == "anc" else IDX) for r in ROOTS}
out = {"delta": {}, "slopes": {}, "bg_c": bg_c, "bg_a": bg_a, "dropped_anc": dropped_anc}
print("\n-- 씨앗 Δ 표 (배경 평균)")
print("   %6s | %-9s | %-9s | %-11s | %-11s | %-11s" % ("μ", "원래", "재료 먼저", "**우주선 끔**", "**재료먼저+끔**", "**둘째 시조**"))
for mu in MUS:
    row = []
    for r in ROOTS:
        m = mean(r, "rascld", mu)
        lo, hi = ci([mean(r, "rascld", mu, ix) for ix in IX[r]])
        out["delta"]["%s|rascld|%g" % (r, mu)] = [m, lo, hi]
        row.append("%+.3f" % m)
    print("   %6g | %s" % (mu, " | ".join("%-11s" % x for x in row)))
for r in ROOTS:
    for c in CODES:
        pt = sl(r, c)
        lo, hi = ci([sl(r, c, ix) for ix in IX[r]])
        out["slopes"]["%s|%s" % (r, c)] = [pt, lo, hi]
print("\n-- 기울기 (Δ / μ · 다섯 점)")
for c in CODES:
    print("   %-7s " % c + " · ".join("%s %+.1f [%+.1f, %+.1f]" % (r, *out["slopes"]["%s|%s" % (r, c)]) for r in ROOTS))

# C1
s_ffc0 = out["slopes"]["ffc0|rascld"]
d_fc = [sl("ffc0", "rascld", ix) - sl("ff", "rascld", ix) for ix in IDX]
lo, hi = ci(d_fc)
if s_ffc0[1] <= 0 <= s_ffc0[2]:
    v1 = "재료 먼저 세계에 남은 불이익은 우주선 몫 — 전부 사라진다"
elif lo > 0:
    v1 = "일부는 우주선 몫"
else:
    v1 = "우주선 몫 아님"
print("\n-- C1 재료 먼저 + 우주선 끔: 씨앗 기울기 %+.1f [%+.1f, %+.1f] · 차(ffc0 − ff) %+.1f [%+.1f, %+.1f] → **%s**"
      % (*s_ffc0, sl("ffc0", "rascld") - sl("ff", "rascld"), lo, hi, v1))
out["C1"] = {"slope": s_ffc0, "diff": [sl("ffc0", "rascld") - sl("ff", "rascld"), lo, hi], "verdict": v1}
# C2
d_co = [sl("c0", "rascld", ix) - sl("orig", "rascld", ix) for ix in IDX]
lo2, hi2 = ci(d_co)
v2 = "원래 세계의 불이익에 우주선 몫은 검출되지 않는다" if lo2 <= 0 <= hi2 else ("일부 우주선 몫" if lo2 > 0 else "우주선을 끄면 더 깎인다")
print("-- C2 원래 + 우주선 끔: 씨앗 기울기 차(c0 − orig) %+.1f [%+.1f, %+.1f] → **%s**" % (sl("c0", "rascld") - sl("orig", "rascld"), lo2, hi2, v2))
out["C2"] = {"diff": [sl("c0", "rascld") - sl("orig", "rascld"), lo2, hi2], "verdict": v2}
# C3
d0a = ci([mean("c0", "rascld", 0.0, ix) - mean("orig", "rascld", 0.0, ix) for ix in IDX])
d0b = ci([mean("ffc0", "rascld", 0.0, ix) - mean("ff", "rascld", 0.0, ix) for ix in IDX])
ok3 = d0a[0] <= 0 <= d0a[1] and d0b[0] <= 0 <= d0b[1]
print("-- C3 버그 없음 점검: μ 0 씨앗 Δ 차 c0−orig [%+.3f, %+.3f] · ffc0−ff [%+.3f, %+.3f] → %s · (서술) racldx 기울기 c0 %+.1f · ffc0 %+.1f"
      % (*d0a, *d0b, "✅ 통과" if ok3 else "🚩 실패", out["slopes"]["c0|racldx"][0], out["slopes"]["ffc0|racldx"][0]))
out["C3"] = {"mu0_c0": list(d0a), "mu0_ffc0": list(d0b), "passes": ok3}
# G1
m0 = out["delta"]["anc|rascld|0"]
sa = out["slopes"]["anc|rascld"]
okA, okB = m0[1] > 0, sa[2] < 0
v4 = "F4 의 두 전제(μ 0 씨앗 우위 · μ 불이익)가 다른 시조 세계에서 재현" if okA and okB else ("절반 재현(%s)" % ("μ 0 우위만" if okA else "μ 불이익만") if okA or okB else "재현 안 됨")
tops = {}
for s in bg_a:
    t = F[("anc", 0.0, s)]["grow_top"][0][0]
    tops[t] = tops.get(t, 0) + 1
print("-- G1 둘째 시조(racld 시조 세계 · 배경 %d): 씨앗 μ 0 Δ %+.3f [%+.3f, %+.3f] · 기울기 %+.1f [%+.1f, %+.1f] → **%s**" % (len(bg_a), *m0, *sa, v4))
print("-- G2 서술: 기울기 비 b_anc/b_orig %.2f · 뺀 배경 %d · 끝 우세 %s · rsacld 기울기 %+.1f · racldx %+.1f"
      % (sa[0] / out["slopes"]["orig|rascld"][0], len(dropped_anc), tops, out["slopes"]["anc|rsacld"][0], out["slopes"]["anc|racldx"][0]))
out["G1"] = {"mu0": m0, "slope": sa, "verdict": v4, "ratio": sa[0] / out["slopes"]["orig|rascld"][0], "tops": tops}
json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_inv10.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_inv10.json")
