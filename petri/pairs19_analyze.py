#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 19 판정 — 꽉 찬 접시가 속도 읽기를 삼키나 (수명 사다리 × 밀도 × 짝 · μ 0). 사전등록 해부19-사전등록-2026-09-28.md §4.
  자료: 수명 150 · 600 = {ROOT}D_<밑>/d<D>/a<A>/dms_mat_<밑>_s*_list[_d32]_gmu0.json (이 실험)
        수명 300     = {REF}D_<밑>/d<D>/mu0/dms_mat_<밑>_s*_list[_d32]_gmu0.json (해부 18 기록 · 재사용 · 같은 시드 3~22)
  G0  설정(age0 = age_var = A · grow 50000·f · ticks 3000·f · every 250·f · 시조 racld · 성장 μ 0 · 밀도 · μ 0) · 복제 충실도 · 배경 멸종 없음
      · 정상 상태: **밑 코드 racld 의 ref 팔**(WT = racld 끼움 = 배경 그대로)의 창 안 개체 수 드리프트(끝 − 처음)/처음 — 파일 중앙값 |·| ≤ DRIFT_MAX
         (🔀 파일럿 뒤 정의 수정 2회: mut 팔 → ref 팔 → racld 의 ref 팔. `e` 가족(rascled)은 WT 자신이 먹는 코드라 ref 팔도 개체 수를 30% 줄인다 — 배경 검사가 아니다)
      · 같은 시드 · 밀도 · 수명의 여섯 밑 코드 파일은 grow_checksum 이 전부 같다(--wt 는 배경에 안 닿는다)
  P1  조작 확인: 밀도 32 빈 칸 몫 v = 1 − pop/칸 — 시드 짝 비 v(150)/v(300) 하한 ≥ V150_MIN · v(600)/v(300) 상한 ≤ V600_MAX
  V1  밀도 32 합친 |ΔΔ(0)|(시드마다 짝 60 평균): 150 − 300 시드 짝 차 하한 > 0   ·   V1b 300 − 600 하한 > 0
  V2  대비 [(150 − 300)₃₂ − (150 − 300)₈] 하한 > 0   ·   V2b [(300 − 600)₃₂ − (300 − 600)₈] 하한 > 0
  V3  서술: 밀도 8 줄 150 − 300 · 600 − 300 (|평균| ≤ FLAT_TOL 이면 '평평')
  V4  서술: 회복 몫 ρ = (150 − 300)₃₂ / (300₈ − 300₃₂) (시드마다 · 부트스트랩)
  V5  서술: 밑 코드별 (150 − 300)₃₂ 하한 > 0 수 · 밀도 32 수명 150 짝의 ΔΔ(0) 상한 < 0 수 · 각 칸의 상위 코드 · 개체 수
  V6  서술(🟡 파일럿 뒤 추가 · 판정 아님): 두 밀도 각각 (600 − 150) 의 시드 짝 차 — 파일럿(2 배경)이 예측과 반대(수명이 길수록 |ΔΔ(0)| 이 큼)를 시사했다
실행: python3 pairs19_analyze.py [petri]  (환경 PETRI_ROOT PETRI_REF PETRI_NBASES PETRI_NBG PETRI_PILOT)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv19")
REF = os.environ.get("PETRI_REF", "inv18")
NBASES = os.environ.get("PETRI_NBASES", "racld rascld rsacld racldx acld rascled").split()
NBG = int(os.environ.get("PETRI_NBG", "20"))
PILOT = os.environ.get("PETRI_PILOT") == "1"
AGES = [150, 300, 600]
DENS = ["8", "32"]
CELLS = 10732
GT = "_gmu0"
def suffix(d): return ("" if d == "8" else "_d" + d) + GT
DRIFT_MAX, V150_MIN, V600_MAX, FLAT_TOL = 0.05, 1.5, 0.67, 0.10   # 🟡 후보 — 파일럿 뒤 §4b 에서 고정
N_BOOT = 2000
RNG = random.Random(20261003)
bad, notes, EXTINCT = [], [], []
OUT = {"root": ROOT, "ref": REF, "pilot": PILOT, "thresholds": {"DRIFT_MAX": DRIFT_MAX, "V150_MIN": V150_MIN, "V600_MAX": V600_MAX, "FLAT_TOL": FLAT_TOL}}


def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None: out.setdefault(k, L)
    return out


def s_of(arm):
    r = arm["series"]; n0, m0 = r[0][1], r[0][2]; nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot_mean(v):
    if not v: return (float("nan"),) * 3
    m = sum(v) / len(v)
    bs = [sum(v[RNG.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fmt(m, lo, hi): return "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)


def pairs_of(codes):
    out = []; Ls = sorted(set(codes.values()))
    for L1, L2 in zip(Ls, Ls[1:]):
        if L2 - L1 != 1: continue
        for a in codes:
            for b in codes:
                if codes[a] == L1 and codes[b] == L2: out.append((a, b))
    return out


def read_dms(dirpath, b, d, a):
    """seed -> {code: s} · seed -> (top, pop) · seed -> 빈 칸 몫 · seed -> 드리프트 중앙값"""
    S, TOP, VAC, DRIFT, CK = {}, {}, {}, {}, {}
    f = a / 300.0
    for fn in sorted(glob.glob(os.path.join(dirpath, "dms_mat_%s_s*_list%s.json" % (b, suffix(d))))):
        z = json.load(open(fn)); s = z["seed"]; tag = "%s d%s a%d s%d" % (b, d, a, s)
        if z["wt"] != b or not z["fidelity"]["same"]: bad.append("dms 설정/충실도 " + tag); continue
        if z["grow_extinct"] >= 0: EXTINCT.append((b, d, a, s, z["grow_extinct"])); continue   # 🟡 09-28 15:50 판정 전 수정: 멸종 배경은 개수와 함께 제외(해부 17 과 같은 처리) — 그 시드의 짝 값이 없으니 판정에 못 들어간다
        exp_age = None if a == 300 else a
        meta_ok = (z.get("age0") == exp_age and z.get("age_var") == exp_age and z["grow"] == round(50000 * f) and z["ticks"] == round(3000 * f)
                   and z["every"] == round(250 * f) and z.get("ancestor") == "racld" and z.get("grow_mu") == 0
                   and z.get("density") == (None if d == "8" else float(d)) and z["mu_assay"] == 0 and z["arms_mode"] == "list")
        if not meta_ok: bad.append("설정 불일치 " + tag + " " + str({k: z.get(k) for k in ("age0", "age_var", "grow", "ticks", "every", "ancestor", "grow_mu", "density", "mu_assay")})); continue
        TOP[s] = (z["grow_top"][0][0], z["grow_pop"]); VAC[s] = 1.0 - z["grow_pop"] / CELLS
        S[s] = {arm["key"]: s_of(arm) for arm in z["arms"] if arm["kind"] == "mut"}
        CK[s] = z["grow_checksum"]
        ref = [arm for arm in z["arms"] if arm["kind"] == "ref"]
        if ref and b == "racld": DRIFT[s] = (ref[0]["series"][-1][1] - ref[0]["series"][0][1]) / ref[0]["series"][0][1]
    return S, TOP, VAC, DRIFT, CK


# ---- 읽기 ----
S, TOPS, VACS, DRIFTS, CKS = {}, {}, {}, {}, {}
for b in NBASES:
    for d in DENS:
        for a in AGES:
            dirpath = os.path.join(BASE, "%sD_%s" % (REF, b), "d" + d, "mu0") if a == 300 else os.path.join(BASE, "%sD_%s" % (ROOT, b), "d" + d, "a%d" % a)
            s_, top, vac, dr, ck = read_dms(dirpath, b, d, a)
            S[(b, d, a)], TOPS[(b, d, a)], VACS[(b, d, a)], DRIFTS[(b, d, a)], CKS[(b, d, a)] = s_, top, vac, dr, ck
            n_ext = sum(1 for e in EXTINCT if e[:3] == (b, d, a))
            if len(s_) + n_ext != NBG: bad.append("파일 수 %s d%s a%d %d + 멸종 %d ≠ %d" % (b, d, a, len(s_), n_ext, NBG))
ck_n, ck_bad = 0, 0
for d in DENS:
    for a in AGES:
        seeds = set().union(*[set(CKS[(b, d, a)]) for b in NBASES])
        for s in seeds:
            vals = [CKS[(b, d, a)][s] for b in NBASES if s in CKS[(b, d, a)]]
            if len(vals) > 1: ck_n += 1; ck_bad += len(set(vals)) != 1
if ck_bad: bad.append("grow_checksum 이 밑 코드 사이에 다름 %d/%d" % (ck_bad, ck_n))
drift_all = {(d, a): [v for b in NBASES for v in DRIFTS[(b, d, a)].values()] for d in DENS for a in AGES}
drift_med = {k: (pct([abs(x) for x in v], .5) if v else float("nan")) for k, v in drift_all.items()}
drift_max = {k: (max(abs(x) for x in v) if v else float("nan")) for k, v in drift_all.items()}
g0_drift = all(v <= DRIFT_MAX for v in drift_med.values() if v == v)
print("== 해부 19 %s— 수명 사다리 × 밀도 × 짝 · 루트 %s · 수명 300 = %s · 밑 코드 %s · 배경 %d" % ("파일럿 " if PILOT else "판정 ", ROOT, REF, " ".join(NBASES), NBG))
print("  G0 설정·충실도·멸종 문제 %d · 밑 코드 사이 grow_checksum 같음 %d/%d · 정상 상태(racld ref 팔) |드리프트| 중앙값(최대): %s %s" % (len(bad), ck_n - ck_bad, ck_n, " · ".join("d%s a%d %.3f (%.3f)" % (d, a, drift_med[(d, a)], drift_max[(d, a)]) for d in DENS for a in AGES), "✅" if g0_drift else "❌"))
for x in bad[:12]: print("    ", x)
ext_bg = sorted(set((e[1], e[2], e[3]) for e in EXTINCT))
print("  멸종 배경(제외) %d: %s" % (len(ext_bg), " · ".join("d%s a%d s%d" % e for e in ext_bg) or "없음"))
OUT["extinct"] = [list(e) for e in EXTINCT]
OUT["G0"] = {"bad": bad, "ck_same": ck_n - ck_bad, "ck_n": ck_n, "drift_med": {"d%s_a%d" % k: v for k, v in drift_med.items()}, "drift_max": {"d%s_a%d" % k: v for k, v in drift_max.items()}, "drift_pass": g0_drift}
if bad and not PILOT: print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)

# ---- V5 서술 · P1 ----
print("\n== 칸마다 상위 코드 · 개체 수 · 빈 칸 몫 (칸 %d)" % CELLS)
vac_seed = {}
for d in DENS:
    for a in AGES:
        tops = {}; pops = []; vac = {}
        for b in NBASES:
            for s, (t, p) in TOPS[(b, d, a)].items(): tops[t] = tops.get(t, 0) + 1; pops.append(p)
            for s, v in VACS[(b, d, a)].items(): vac.setdefault(s, []).append(v)
        vac_seed[(d, a)] = {s: sum(v) / len(v) for s, v in vac.items()}
        if pops:
            vm = [x for x in vac_seed[(d, a)].values()]
            print("   밀도 %-2s 수명 %-3d: 상위 코드 %s · 개체 수 중앙 %d [%d, %d] · 빈 칸 몫 중앙 %.3f" % (d, a, dict(sorted(tops.items(), key=lambda kv: -kv[1])), pct(pops, .5), min(pops), max(pops), pct(vm, .5)))
            OUT["cell|d%s|a%d" % (d, a)] = {"tops": tops, "pop_median": pct(pops, .5), "pop_min": min(pops), "pop_max": max(pops), "vac_median": pct(vm, .5)}
def ratio_pair(d, a1, a2):
    v1, v2 = vac_seed.get((d, a1), {}), vac_seed.get((d, a2), {})
    ss = sorted(set(v1) & set(v2)); return boot_mean([v1[s] / v2[s] for s in ss if v2[s] > 0]), len(ss)
(r150, r150l, r150h), n150 = ratio_pair("32", 150, 300); (r600, r600l, r600h), n600 = ratio_pair("32", 600, 300)
p1a = r150l >= V150_MIN; p1b = r600h <= V600_MAX
print("   P1 밀도 32 빈 칸 비 v(150)/v(300) %s (하한 ≥ %.2f) %s · v(600)/v(300) %s (상한 ≤ %.2f) %s · 배경 %d/%d" % (fmt(r150, r150l, r150h), V150_MIN, "✅" if p1a else "❌", fmt(r600, r600l, r600h), V600_MAX, "✅" if p1b else "❌", n150, n600))
OUT["P1"] = {"r150": r150, "ci150": [r150l, r150h], "r600": r600, "ci600": [r600l, r600h], "pass_a": p1a, "pass_b": p1b, "pass": p1a and p1b}

# ---- 합친 |ΔΔ(0)| ----
def dd(b, d, a, x, y):
    S_ = S[(b, d, a)]; return {s: S_[s][y] - S_[s][x] for s in S_ if x in S_[s] and y in S_[s]}
A = {}   # (d, a) -> seed -> 합친 |ΔΔ(0)|
Ab = {}  # (b, d, a) -> seed -> 밑 코드 안 평균
sign = {}
for d in DENS:
    for a in AGES:
        per = {}
        for b in NBASES:
            codes = inserts(b, "n"); perb = {}
            for x, y in pairs_of(codes):
                v = dd(b, d, a, x, y)
                for s, val in v.items(): per.setdefault(s, []).append(abs(val)); perb.setdefault(s, []).append(abs(val))
                if d == "32" and a == 150 and v:
                    m, lo, hi = boot_mean([v[s] for s in sorted(v)]); sign[(b, x, y)] = (m, lo, hi, hi < 0)
            Ab[(b, d, a)] = {s: sum(x) / len(x) for s, x in perb.items()}
        A[(d, a)] = {s: sum(x) / len(x) for s, x in per.items()}
print("\n== 합친 |ΔΔ(0)| (시드마다 짝 평균 → 시드 부트스트랩)")
for d in DENS:
    row = []
    for a in AGES:
        m, lo, hi = boot_mean([A[(d, a)][s] for s in sorted(A[(d, a)])]); row.append("수명 %d %s (배경 %d)" % (a, fmt(m, lo, hi), len(A[(d, a)])))
        OUT["A|d%s|a%d" % (d, a)] = {"mean": m, "ci": [lo, hi], "n": len(A[(d, a)])}
    print("   밀도 %-2s: %s" % (d, " · ".join(row)))
def pdiff(d, a1, a2):
    x, y = A[(d, a1)], A[(d, a2)]; ss = sorted(set(x) & set(y)); return {s: x[s] - y[s] for s in ss}
def contrast(a1, a2):
    p32, p8 = pdiff("32", a1, a2), pdiff("8", a1, a2); ss = sorted(set(p32) & set(p8)); return boot_mean([p32[s] - p8[s] for s in ss]), len(ss)
d32_150 = pdiff("32", 150, 300); d32_600 = pdiff("32", 300, 600); d8_150 = pdiff("8", 150, 300); d8_600 = pdiff("8", 600, 300)
(v1, v1l, v1h) = boot_mean([d32_150[s] for s in sorted(d32_150)]); (v1b, v1bl, v1bh) = boot_mean([d32_600[s] for s in sorted(d32_600)])
(v2, v2l, v2h), n2 = contrast(150, 300); (v2b, v2bl, v2bh), n2b = contrast(300, 600)
(v3a, v3al, v3ah) = boot_mean([d8_150[s] for s in sorted(d8_150)]); (v3b, v3bl, v3bh) = boot_mean([d8_600[s] for s in sorted(d8_600)])
V1 = v1l > 0; V1b = v1bl > 0; V2 = v2l > 0; V2b = v2bl > 0
flat = abs(v3a) <= FLAT_TOL and abs(v3b) <= FLAT_TOL
print("\n== V1 — 밀도 32: 150 − 300 %s %s   ·   V1b 300 − 600 %s %s   (배경 %d · %d)" % (fmt(v1, v1l, v1h), "✅" if V1 else "❌", fmt(v1b, v1bl, v1bh), "✅" if V1b else "❌", len(d32_150), len(d32_600)))
print("== V2 — 대비 (150 − 300)₃₂ − (150 − 300)₈ %s %s   ·   V2b (300 − 600)₃₂ − (300 − 600)₈ %s %s   (배경 %d · %d)" % (fmt(v2, v2l, v2h), "✅" if V2 else "❌", fmt(v2b, v2bl, v2bh), "✅" if V2b else "❌", n2, n2b))
print("== V3 서술 — 밀도 8 줄: 150 − 300 %s · 600 − 300 %s → %s (|평균| ≤ %.2f)" % (fmt(v3a, v3al, v3ah), fmt(v3b, v3bl, v3bh), "평평" if flat else "평평하지 않음", FLAT_TOL))
# V4 회복 몫
a8_300, a32_300, a32_150 = A[("8", 300)], A[("32", 300)], A[("32", 150)]
ss = sorted(set(a8_300) & set(a32_300) & set(a32_150))
rho_v = [(a32_150[s] - a32_300[s]) / (a8_300[s] - a32_300[s]) for s in ss if (a8_300[s] - a32_300[s]) > 0]
(rho, rhol, rhoh) = boot_mean(rho_v)
gap = boot_mean([a8_300[s] - a32_300[s] for s in ss])
print("== V4 서술 — 회복 몫 ρ = (150 − 300)₃₂ / (300₈ − 300₃₂) %s (배경 %d · 틈 300₈ − 300₃₂ %s)" % (fmt(rho, rhol, rhoh), len(rho_v), fmt(*gap)))
# V5 밑 코드별 · 부호
nb_pass = 0; per_base = {}
for b in NBASES:
    x, y = Ab[(b, "32", 150)], Ab[(b, "32", 300)]; ss = sorted(set(x) & set(y))
    m, lo, hi = boot_mean([x[s] - y[s] for s in ss]); ok = lo > 0; nb_pass += ok
    x8, y8 = Ab[(b, "8", 150)], Ab[(b, "8", 300)]; s8 = sorted(set(x8) & set(y8)); m8, lo8, hi8 = boot_mean([x8[s] - y8[s] for s in s8])
    per_base[b] = {"d32_150_300": m, "ci": [lo, hi], "pass": ok, "d8_150_300": m8, "ci8": [lo8, hi8], "n": len(ss)}
    print("   V5 %-8s 밀도 32 150 − 300 %s %s · 밀도 8 150 − 300 %s · 배경 %d" % (b, fmt(m, lo, hi), "✅" if ok else "❌", fmt(m8, lo8, hi8), len(ss)))
nsign = sum(1 for v in sign.values() if v[3])
print("   V5 밑 코드 %d/%d 에서 밀도 32 150 − 300 하한 > 0 · 밀도 32 수명 150 짝 ΔΔ(0) 상한 < 0: %d/%d" % (nb_pass, len(NBASES), nsign, len(sign)))
for (b, x, y), (m, lo, hi, ok) in sign.items():
    if not ok: print("      부호 실패 %s %s→%s %s" % (b, x, y, fmt(m, lo, hi)))
d32_span = pdiff("32", 600, 150); d8_span = pdiff("8", 600, 150)
(v6a, v6al, v6ah) = boot_mean([d32_span[s] for s in sorted(d32_span)]); (v6b, v6bl, v6bh) = boot_mean([d8_span[s] for s in sorted(d8_span)])
print("== V6 서술(🟡 판정 아님) — (600 − 150): 밀도 32 %s · 밀도 8 %s" % (fmt(v6a, v6al, v6ah), fmt(v6b, v6bl, v6bh)))
OUT["V6"] = {"d32_600_150": v6a, "ci32": [v6al, v6ah], "d8_600_150": v6b, "ci8": [v6bl, v6bh]}
OUT.update({"V1": {"diff": v1, "ci": [v1l, v1h], "pass": V1}, "V1b": {"diff": v1b, "ci": [v1bl, v1bh], "pass": V1b},
            "V2": {"contrast": v2, "ci": [v2l, v2h], "pass": V2, "n": n2}, "V2b": {"contrast": v2b, "ci": [v2bl, v2bh], "pass": V2b, "n": n2b},
            "V3": {"d8_150_300": v3a, "ci_a": [v3al, v3ah], "d8_600_300": v3b, "ci_b": [v3bl, v3bh], "flat": flat},
            "V4": {"rho": rho, "ci": [rhol, rhoh], "n": len(rho_v), "gap": gap[0], "gap_ci": [gap[1], gap[2]]},
            "V5": {"per_base": per_base, "n_base_pass": nb_pass, "sign_pass": nsign, "sign_n": len(sign)}})
print("\n== 종합: G0 %s · P1 %s · V1 %s · V1b %s · V2 %s · V2b %s · (서술) V3 %s · ρ %.2f · V5 %d/%d" % ("✅" if (not bad and g0_drift) else "❌", "✅" if (p1a and p1b) else "❌", "✅" if V1 else "❌", "✅" if V1b else "❌", "✅" if V2 else "❌", "✅" if V2b else "❌", "평평" if flat else "안 평평", rho, nb_pass, len(NBASES)))
OUT["bad"] = bad
fn = "_RESULT_pilot19.json" if PILOT else "_RESULT_pairs19.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
