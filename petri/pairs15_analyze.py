#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 15 판정 — 짝 안의 μ 사다리(B) · 세 규칙 × 짝(C) · 새 배경 재현(E)
사전등록: 해부15-사전등록-2026-09-27.md §1 · §4. 판정선(§4b)은 파일럿 뒤 고정 — 고정 전에는 PETRI_PILOT=1 로만 돈다.

  ΔΔ = s(고리 L+1 코드) − s(고리 L 코드)  (해부 13 과 같은 부호 · 기준 팔과 무관 · 배경 짝 · 부트스트랩 2,000)
  자료 출처(μ · 규칙별):
     원래 규칙 μ 0~1%     inv15B_<밑>/muX (목록 모드 · --wt 밑 코드 · 5 μ 전부 새로 돌린다)
     회귀 검사만          첫 가족 μ 0 inv12nbr_* · 둘째 가족 μ 0 · 0.003 inv13nbr_* (이웃 모드 · --wt racld) — 목록 팔 = 이웃 팔(0~3열)
     재료 먼저 · 기억     inv15C_<밑>/<ff|rem>/muX
     새 배경              inv15E_*(Δ) · inv15Eib_*(invbud R · T)

  G0   목록 모드 팔 기록(0~3열) = 이웃 모드 팔 기록(μ 0 여섯 · 둘째 가족 μ 0.003) · C 의 `racld` ref 팔 = inv7ff · inv9mem 의 ref 팔
  B1   짝마다 ΔΔ(μ) 의 μ 기울기(배경별 최소제곱 · 5 μ) 하한 > 0 — 60 짝
  B2   서술: 교차 μ*(선형 보간)
  C1   규칙마다 **짝 60 · 배경 20 을 합친** (ΔΔ_규칙(0) − ΔΔ_원래(0)) 의 평균: 구간이 0 을 품고 |평균| ≤ C1_TOL  (배경별로 짝 평균 → 배경 재추출)
  C2   기억 규칙의 합친 Δμ = ΔΔ(0.003) − ΔΔ(0): 구간이 0 을 품고 |평균| ≤ C2_TOL · 원래 규칙의 합친 Δμ 하한 > 0 · 대비(원래 − 기억) 하한 > 0
       (짝별 값은 서술 — 파일럿 §5b: 짝 하나의 ΔΔ 는 배경당 로그비 하나라 잡음이 크다)
  C3   서술: 재료 먼저의 같은 양
  E1   새 배경 N1(상한 < 0) · N2(틱당 [−1.0, −0.2]) — 해부 13 판정선
  E2   서술: 새 배경 N3(0.15) 실패 쌍 · 해부 13 실패 쌍과 겹침
  E3   새 배경 A1(하한 > 0) · A2a(|ln R 비| ≤ 0.10) · A2b(|T/R 차| ≤ 0.10) · 자기 검산 — 해부 14 판정선
실행: python3 pairs15_analyze.py [petri 폴더]  (환경: PETRI_ROOT PETRI_NBASES PETRI_NBG PETRI_PILOT PETRI_REF12 PETRI_REF13 PETRI_REF7 PETRI_REF9)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv15")
NBASES = os.environ.get("PETRI_NBASES", "racld rascld rsacld racldx acld rascled").split()
NBG = int(os.environ.get("PETRI_NBG", "20"))
PILOT = os.environ.get("PETRI_PILOT") == "1"
REF12 = os.environ.get("PETRI_REF12", "inv12"); REF13 = os.environ.get("PETRI_REF13", "inv13")
REF7 = os.environ.get("PETRI_REF7", "inv7ff"); REF9 = os.environ.get("PETRI_REF9", "inv9mem")
FAM1 = {"racld", "rascld", "rsacld", "racldx"}
MUS = [0.0, 0.001, 0.003, 0.005, 0.01]
# ---- 판정선 (§4b 에서 고정 · 고정 전 None) ----
C1_TOL = 0.05         # §4b 09-27 파일럿 뒤 고정 — 합친 평균의 허용(배경 20 · 짝 60 이면 구간 반폭 ≈ 0.07 예상)
C2_TOL = 0.05
TOL_NEG, PT_LO, PT_HI = 0.15, -1.0, -0.2        # 해부 13 그대로
TAU, TR_PAIR, TR_LO, TR_HI = 0.10, 0.10, -0.5, 1.5   # 해부 14 그대로
N_BOOT = 2000
RNG = random.Random(20260929)
N3_FAIL13 = {("acnld", "ancld"), ("acnld", "nacld"), ("nrascled", "rasclend"), ("ranscled", "rasclend"), ("rascledn", "rasclend"),
             ("rasclend", "rasclned"), ("rasclend", "rnascled"), ("rascldn", "rnascld"), ("rasclnd", "rnascld")}
N3_FAIL13 |= {(b, a) for a, b in N3_FAIL13}
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT, "C1_TOL": C1_TOL, "C2_TOL": C2_TOL}


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def s_of(arm):
    r = arm["series"]; n0, m0 = r[0][1], r[0][2]; nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot_mean(v):
    m = sum(v) / len(v)
    bs = [sum(v[RNG.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fmt(m, lo, hi):
    return "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def pairs_of(codes):
    """(짧은 코드, 긴 코드) — 이웃 등급 · 틱 차 1 인 짝만(해부 13 과 같다)"""
    out = []
    Ls = sorted(set(codes.values()))
    for L1, L2 in zip(Ls, Ls[1:]):
        if L2 - L1 != 1: continue
        for a in codes:
            for b in codes:
                if codes[a] == L1 and codes[b] == L2: out.append((a, b))
    return out


def same_loop_pairs(codes):
    cl = sorted(codes); out = []
    for i, a in enumerate(cl):
        for b in cl[i + 1:]:
            if codes[a] == codes[b]: out.append((a, b))
    return out


# ---------------- 읽기: s 값 표 ----------------
# S[(b, rule, mu)] = {seed: {code: s}} · SER[(b, rule, mu)] = {seed: {code: (series, stopped)}}
S, SER, REFSER = {}, {}, {}


def read_nbr(root, b, mu):
    d, ser = {}, {}
    for f in sorted(glob.glob(os.path.join(BASE, "%snbr_%s" % (root, b), mtag(mu), "dms_mat_racld_s*_nbr_%s.json" % b))):
        z = json.load(open(f)); s = z["seed"]
        if abs(z["mu_assay"] - mu) > 1e-12 or z.get("nbr_of") != b or not z["checksum_match"] or not z["fidelity"]["same"]:
            bad.append("이웃 기록 설정 %s %s %d" % (root, b, s)); continue
        d[s] = {a["key"]: s_of(a) for a in z["arms"] if a["kind"] == "mut"}
        ser[s] = {a["key"]: (a["series"], a["stopped"]) for a in z["arms"] if a["kind"] == "mut"}
    return d, ser


def read_list(dirpath, b, mu, rule):
    d, ser, refser = {}, {}, {}
    codes = inserts(b, "n")
    fs = sorted(glob.glob(os.path.join(BASE, dirpath, "dms_mat_%s_s*_list.json" % b)))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        okflag = (bool(z.get("find_first")) == (rule == "ff")) and (bool(z.get("remember_die")) == (rule == "rem"))
        if abs(z["mu_assay"] - mu) > 1e-12 or z["wt"] != b or not okflag: bad.append("목록 설정 %s %d" % (dirpath, s)); continue
        if not z["checksum_match"]: bad.append("배경 체크섬 %s %d" % (dirpath, s))
        if not z["fidelity"]["same"]: bad.append("충실도 %s %d" % (dirpath, s))
        muts = [a for a in z["arms"] if a["kind"] == "mut"]
        if [a["key"] for a in muts] != list(codes) or len(z["arms"]) != 11 + len(codes): bad.append("팔 목록 %s %d" % (dirpath, s))
        d[s] = {a["key"]: s_of(a) for a in muts}
        ser[s] = {a["key"]: (a["series"], a["stopped"]) for a in muts}
        refser[s] = [(a["series"], a["stopped"]) for a in z["arms"] if a["kind"] == "ref"][0]
    return d, ser, refser, len(fs)


for b in NBASES:
    fam1 = b in FAM1
    # 원래 규칙 — 5 μ 전부 목록 모드(inv15B). 이웃 기록은 G0 대조에만 쓴다.
    for mu in MUS:
        d, ser, _, n = read_list("%sB_%s/%s" % (ROOT, b, mtag(mu)), b, mu, "orig")
        if n != NBG: bad.append("B 파일 수 %s %s %d/%d" % (b, mu, n, NBG))
        S[(b, "orig", mu)] = d
        refroot = (REF12 if fam1 else REF13) if (mu == 0.0 or (not fam1 and mu == 0.003)) else None
        if refroot:
            dn, sern = read_nbr(refroot, b, mu)
            chk = mm = 0
            for sd_ in ser:
                for c in ser[sd_]:
                    if sd_ in sern and c in sern[sd_]:
                        chk += 1
                        if [r[:4] for r in ser[sd_][c][0]] != [r[:4] for r in sern[sd_][c][0]] or ser[sd_][c][1] != sern[sd_][c][1]: mm += 1
            OUT["G0|list=nbr|%s|%s" % (b, mu)] = {"checked": chk, "mismatch": mm}
            if mm: bad.append("G0 목록≠이웃 %s %s (%d/%d)" % (b, mu, mm, chk))
            if not PILOT and chk == 0: notes.append("G0 기준 팔 없음 %s %s" % (b, mu))
    # 규칙 세계
    for rule in ("ff", "rem"):
        for mu in (0.0, 0.003):
            d, ser, refser, n = read_list("%sC_%s/%s/%s" % (ROOT, b, rule, mtag(mu)), b, mu, rule)
            if n != NBG: bad.append("C 파일 수 %s %s %s %d/%d" % (b, rule, mu, n, NBG))
            S[(b, rule, mu)] = d
            if b == "racld":   # G0: ref 팔(racld 끼움) = 해부 7 · 9 의 ref 팔
                root = REF7 if rule == "ff" else REF9
                chk = mm = 0
                for s, (series, stopped) in refser.items():
                    cands = glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s%05d_*.json" % s)) + \
                            glob.glob(os.path.join(BASE, root, "dms_mat_racld_s%05d_*.json" % s))
                    if not cands: continue
                    z = json.load(open(cands[0])); ra = [a for a in z["arms"] if a["kind"] == "ref"][0]
                    chk += 1
                    if [r[:3] for r in ra["series"]] != [r[:3] for r in series] or ra["stopped"] != stopped: mm += 1
                OUT["G0|ref|%s|%s" % (rule, mu)] = {"checked": chk, "mismatch": mm}
                if mm: bad.append("G0 규칙 ref 팔 불일치 %s %s (%d/%d)" % (rule, mu, mm, chk))
    # 새 배경 Δ
    d, ser, _, n = read_list("%sE_%s/mu0" % (ROOT, b), b, 0.0, "orig")
    if n != NBG: bad.append("E 파일 수 %s %d/%d" % (b, n, NBG))
    S[(b, "E", 0.0)] = d

# 새 배경 invbud
IB = {}
KEYS_SUM = ["copy_ok", "copy_stall", "copy_del", "copy_phase"]
for b in NBASES:
    codes = inserts(b, "n"); IB[b] = {c: {} for c in codes}
    fs = sorted(glob.glob(os.path.join(BASE, "%sEib_%s" % (ROOT, b), "mu0", "ib_inv_s*_mu0.json")))
    if len(fs) != NBG: bad.append("Eib 파일 수 %s %d/%d" % (b, len(fs), NBG))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if not z["checksum_match"] or list(z["codes"]) != list(codes): bad.append("Eib 설정 %s %d" % (b, s))
        for a in z["arms"]:
            if a["kind"] != "mut": continue
            d = {k: sum(bn[k][1] for bn in a["bins"]) for k in KEYS_SUM}
            if a["stopped"] >= 0 and a["series"][-1][2] == 0: continue
            if d["copy_ok"] == 0: continue
            IB[b][a["key"]][s] = {"R": (d["copy_ok"] + d["copy_stall"] + d["copy_del"]) / d["copy_ok"], "T": d["copy_phase"] / d["copy_ok"], **d}

print("== 해부 15 %s— 루트 %s · 밑 코드 %s · 배경 %d" % ("파일럿 " if PILOT else "판정 ", ROOT, " ".join(NBASES), NBG))
for k, v in OUT.items():
    if k.startswith("G0|"): print("  %s: 대조 %d · 불일치 %d" % (k, v["checked"], v["mismatch"]))
print("  점검 문제 %d" % len(bad))
for x in bad[:15]: print("    ", x)
if bad and not PILOT:
    print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)
if (C1_TOL is None or C2_TOL is None) and not PILOT:
    print("🚨 C1_TOL · C2_TOL 이 고정되지 않았다(§4b) — 판정하지 않는다"); sys.exit(1)


def dd_by_seed(b, rule, mu, a, c):
    """ΔΔ = s(긴 코드 c) − s(짧은 코드 a) · 배경마다"""
    d = S.get((b, rule, mu), {})
    return {s: d[s][c] - d[s][a] for s in d if a in d[s] and c in d[s]}


# ---------------- B ----------------
print("\n== B1 — 짝마다 ΔΔ(μ) 의 μ 기울기 하한 > 0 (배경별 최소제곱 · μ %s)" % MUS)
b1_all, nb1 = True, 0
for b in NBASES:
    codes = inserts(b, "n")
    print("  밑 코드 `%s`" % b)
    for a, c in pairs_of(codes):
        per = {mu: dd_by_seed(b, "orig", mu, a, c) for mu in MUS}
        seeds = sorted(set.intersection(*[set(per[mu]) for mu in MUS])) if all(per[mu] for mu in MUS) else []
        if not seeds: print("     %s→%s: 자료 없음" % (a, c)); continue
        mx = sum(MUS) / len(MUS)
        slopes = []
        for s in seeds:
            ys = [per[mu][s] for mu in MUS]; my = sum(ys) / len(ys)
            slopes.append(sum((mu - mx) * (y - my) for mu, y in zip(MUS, ys)) / sum((mu - mx) ** 2 for mu in MUS))
        m, lo, hi = boot_mean(slopes); ok = lo > 0; b1_all &= ok; nb1 += 1
        means = {mu: sum(per[mu][s] for s in seeds) / len(seeds) for mu in MUS}
        # B2 교차 μ*
        cross = None
        for m1, m2 in zip(MUS, MUS[1:]):
            y1, y2 = means[m1], means[m2]
            if y1 < 0 <= y2: cross = m1 + (0 - y1) * (m2 - m1) / (y2 - y1); break
        print("     %-8s→%-8s 기울기/μ %s %s · ΔΔ(μ) %s · μ* %s" % (a, c, fmt(m, lo, hi), "✅" if ok else "❌",
              " ".join("%+.2f" % means[mu] for mu in MUS), ("%.4f" % cross) if cross is not None else "구간 밖"))
        OUT["B1|%s|%s|%s" % (b, a, c)] = {"slope": m, "ci": [lo, hi], "pass": ok, "dd_by_mu": means, "cross": cross, "n_bg": len(seeds)}
print("  B1 종합: %s (%d 짝)" % ("✅" if b1_all else "❌", nb1)); OUT["B1"] = b1_all

# ---------------- C ----------------
print("\n== C — 규칙 × 짝 (합친 판정: C1 |평균(ΔΔ_규칙(0) − ΔΔ_원래(0))| ≤ %s · 구간이 0 을 품음 / C2 기억 |평균 Δμ| ≤ %s · 원래 Δμ 하한 > 0 · 대비 하한 > 0)" % (C1_TOL, C2_TOL))
c1_abs, c2_abs = [], []
POOL = {}   # (rule, what) -> {seed: [values over pairs]}
def push(rule, what, s_, v):
    POOL.setdefault((rule, what), {}).setdefault(s_, []).append(v)
for b in NBASES:
    codes = inserts(b, "n")
    print("  밑 코드 `%s` (짝별 서술)" % b)
    for a, c in pairs_of(codes):
        o0 = dd_by_seed(b, "orig", 0.0, a, c); o3 = dd_by_seed(b, "orig", 0.003, a, c)
        row = []
        for rule in ("ff", "rem"):
            r0 = dd_by_seed(b, rule, 0.0, a, c); r3 = dd_by_seed(b, rule, 0.003, a, c)
            sd = sorted(set(o0) & set(r0))
            if not sd: continue
            for s_ in sd: push(rule, "c1", s_, r0[s_] - o0[s_])
            m1, lo1, hi1 = boot_mean([r0[s_] - o0[s_] for s_ in sd]); c1_abs.append(abs(m1))
            sd2 = sorted(set(r0) & set(r3))
            for s_ in sd2: push(rule, "dmu", s_, r3[s_] - r0[s_])
            m2, lo2, hi2 = boot_mean([r3[s_] - r0[s_] for s_ in sd2]) if sd2 else (float("nan"),) * 3
            if rule == "rem": c2_abs.append(abs(m2))
            row.append("%s: C1 %s · Δμ %s" % (rule, fmt(m1, lo1, hi1), fmt(m2, lo2, hi2)))
            OUT["C|%s|%s|%s|%s" % (b, rule, a, c)] = {"c1": m1, "c1_ci": [lo1, hi1], "dmu": m2, "dmu_ci": [lo2, hi2], "n_bg": len(sd)}
        sdo = sorted(set(o0) & set(o3))
        if sdo:
            for s_ in sdo: push("orig", "dmu", s_, o3[s_] - o0[s_])
            mo, loo, hio = boot_mean([o3[s_] - o0[s_] for s_ in sdo])
            row.append("orig Δμ %s" % fmt(mo, loo, hio))
            OUT["C|%s|orig|%s|%s" % (b, a, c)] = {"dmu": mo, "dmu_ci": [loo, hio]}
        print("     %-8s→%-8s " % (a, c) + " | ".join(row))
if c1_abs: print("  짝별 C1 |차| 중앙값 %.3f · 최대 %.3f · 90분위 %.3f" % (pct(c1_abs, .5), max(c1_abs), pct(c1_abs, .9))); OUT["C1_abs"] = {"median": pct(c1_abs, .5), "max": max(c1_abs), "p90": pct(c1_abs, .9)}
if c2_abs: print("  짝별 C2(기억) |Δμ| 중앙값 %.3f · 최대 %.3f · 90분위 %.3f" % (pct(c2_abs, .5), max(c2_abs), pct(c2_abs, .9))); OUT["C2_abs"] = {"median": pct(c2_abs, .5), "max": max(c2_abs), "p90": pct(c2_abs, .9)}
def pooled(rule, what):
    d = POOL.get((rule, what), {})
    v = [sum(x) / len(x) for s_, x in sorted(d.items()) if x]
    return boot_mean(v) if v else (float("nan"),) * 3, len(v)
c1_all, c2_all, c2o_all, c2c_all = True, True, True, True
print("  합친 판정(배경별 짝 평균 → 배경 재추출):")
for rule in ("ff", "rem"):
    (m, lo, hi), n = pooled(rule, "c1")
    ok = (lo <= 0 <= hi) and abs(m) <= C1_TOL
    if rule == "rem": c1_all &= ok
    print("     C1 %s: 평균 %s · 배경 %d %s%s" % (rule, fmt(m, lo, hi), n, "✅" if ok else "❌", "" if rule == "rem" else " (ff 는 서술)"))
    OUT["C1|%s" % rule] = {"m": m, "ci": [lo, hi], "pass": ok, "n_bg": n}
(mr, lor, hir), nr = pooled("rem", "dmu"); okr = (lor <= 0 <= hir) and abs(mr) <= C2_TOL; c2_all &= okr
(mo, loo, hio), no = pooled("orig", "dmu"); oko = loo > 0; c2o_all &= oko
(mf, lof, hif), nf = pooled("ff", "dmu")
# 대비 원래 − 기억 (배경 짝)
dr = POOL.get(("rem", "dmu"), {}); do = POOL.get(("orig", "dmu"), {})
cs = sorted(set(dr) & set(do)); (mc, loc, hic) = boot_mean([sum(do[s_]) / len(do[s_]) - sum(dr[s_]) / len(dr[s_]) for s_ in cs]) if cs else (float("nan"),) * 3
okc = loc > 0; c2c_all &= okc
print("     C2 기억 Δμ: 평균 %s · 배경 %d %s" % (fmt(mr, lor, hir), nr, "✅" if okr else "❌"))
print("     C2 원래 Δμ: 평균 %s · 배경 %d %s" % (fmt(mo, loo, hio), no, "✅" if oko else "❌"))
print("     C2 대비(원래 − 기억): %s %s" % (fmt(mc, loc, hic), "✅" if okc else "❌"))
print("     C3 서술 재료 먼저 Δμ: 평균 %s · 배경 %d" % (fmt(mf, lof, hif), nf))
OUT["C2|rem"] = {"m": mr, "ci": [lor, hir], "pass": okr}; OUT["C2|orig"] = {"m": mo, "ci": [loo, hio], "pass": oko}; OUT["C2|contrast"] = {"m": mc, "ci": [loc, hic], "pass": okc}; OUT["C3|ff"] = {"m": mf, "ci": [lof, hif]}
print("  C1 종합: %s · C2 종합: 기억 %s · 원래 %s · 대비 %s" % ("✅" if c1_all else "❌", "✅" if c2_all else "❌", "✅" if c2o_all else "❌", "✅" if c2c_all else "❌"))
OUT["C1"] = c1_all; OUT["C2"] = c2_all and c2o_all and c2c_all; OUT["C2_orig"] = c2o_all

# ---------------- E ----------------
print("\n== E1 · E2 — 새 배경 N1 · N2 · N3(해부 13 판정선)")
e1_all = True; n3f = []
for b in NBASES:
    codes = inserts(b, "n")
    for a, c in pairs_of(codes):
        dd = dd_by_seed(b, "E", 0.0, a, c)
        if not dd: continue
        v = [dd[s] for s in sorted(dd)]; m, lo, hi = boot_mean(v)
        ok = hi < 0 and PT_LO <= m <= PT_HI; e1_all &= ok
        print("   %s %-8s→%-8s ΔΔ %s · 배경 %d %s" % (b, a, c, fmt(m, lo, hi), len(v), "✅" if ok else "❌"))
        OUT["E1|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "pass": ok, "n_bg": len(v)}
    for a, c in same_loop_pairs(codes):
        d = S.get((b, "E", 0.0), {}); v = [d[s][c] - d[s][a] for s in d if a in d[s] and c in d[s]]
        if not v: continue
        m, lo, hi = boot_mean(v); fail = abs(m) > TOL_NEG
        if fail: n3f.append((b, a, c, m))
        OUT["E2|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "fail": fail, "fail13": (a, c) in N3_FAIL13}
print("  E1 종합: %s" % ("✅" if e1_all else "❌")); OUT["E1"] = e1_all
print("  E2 서술: 새 배경 N3 실패 %d 쌍 — %s" % (len(n3f), " · ".join("%s %s/%s %+.3f%s" % (b, a, c, m, "★" if (a, c) in N3_FAIL13 else "") for b, a, c, m in n3f)))
OUT["E2_fail"] = n3f

print("\n== E3 — 새 배경 A1 · A2a · A2b · 자기 검산(해부 14 판정선)")
e3_a1, e3_a2a, e3_a2b, e3_tr = True, True, True, True
for b in NBASES:
    codes = inserts(b, "n"); Ls = sorted(set(codes.values()))
    pooled = {}
    for c in codes:
        ds = IB[b][c].values()
        ok_ = sum(d["copy_ok"] for d in ds); st = sum(d["copy_stall"] + d["copy_del"] for d in ds); ph = sum(d["copy_phase"] for d in ds)
        if ok_: pooled[c] = ((ok_ + st) / ok_, ph / ok_)
    for c, (R, T) in pooled.items():
        if not (TR_LO <= T / R - codes[c] <= TR_HI): e3_tr = False; print("   ❌ 자기 검산 %s T/R−L %+.2f" % (c, T / R - codes[c]))
    seeds = sorted(set.intersection(*[set(IB[b][c]) for c in codes if IB[b][c]])) if any(IB[b][c] for c in codes) else []
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        v = [sum(math.log(IB[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(IB[b][c][s]["R"]) for c in c2) / len(c2) for s in seeds]
        if not v: continue
        m, lo, hi = boot_mean(v); ok = lo > 0; e3_a1 &= ok
        print("   %s A1 ln R(%d) − ln R(%d) = %s %s" % (b, L1, L2, fmt(m, lo, hi), "✅" if ok else "❌"))
        OUT["E3A1|%s|%d|%d" % (b, L1, L2)] = {"d": m, "ci": [lo, hi], "pass": ok}
    for a, c in same_loop_pairs(codes):
        sd = sorted(set(IB[b][a]) & set(IB[b][c]))
        if not sd or a not in pooled or c not in pooled: continue
        m, lo, hi = boot_mean([math.log(IB[b][a][s]["R"]) - math.log(IB[b][c][s]["R"]) for s in sd])
        trd = abs(pooled[a][1] / pooled[a][0] - pooled[c][1] / pooled[c][0])
        oka, okb = abs(m) <= TAU, trd <= TR_PAIR; e3_a2a &= oka; e3_a2b &= okb
        OUT["E3A2|%s|%s|%s" % (b, a, c)] = {"d": m, "ci": [lo, hi], "tr_diff": trd, "pass_a": oka, "pass_b": okb}
        if not (oka and okb): print("   ❌ %s %s vs %s ln R 비 %s · T/R 차 %.3f" % (b, a, c, fmt(m, lo, hi), trd))
print("  E3 종합: A1 %s · A2a %s · A2b %s · 자기 검산 %s" % ("✅" if e3_a1 else "❌", "✅" if e3_a2a else "❌", "✅" if e3_a2b else "❌", "✅" if e3_tr else "❌"))
OUT.update({"E3_A1": e3_a1, "E3_A2a": e3_a2a, "E3_A2b": e3_a2b, "E3_tr": e3_tr})

print("\n== 종합")
if PILOT:
    print("  파일럿 — 판정 없음. C1 · C2 문턱 후보는 위 |차| 분포 줄.")
else:
    print("  B1 %s · C1 %s · C2 %s · E1 %s · E3 %s" % tuple("✅" if x else "❌" for x in (OUT["B1"], OUT["C1"], OUT["C2"], OUT["E1"], e3_a1 and e3_a2a and e3_a2b and e3_tr)))
OUT["bad"] = bad; OUT["notes"] = notes
fn = "_RESULT_pilot15.json" if PILOT else "_RESULT_pairs15.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
