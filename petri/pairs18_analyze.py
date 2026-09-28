#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 18 판정 — 재료가 넉넉하면(밀도 32) 고리 한 틱의 두 얼굴 중 무엇이 남나. 사전등록 해부18-사전등록-2026-09-28.md §4.
  🔀 §4b 설계 변경 뒤: 두 밀도 모두 순수 racld 공동체(--ancestor racld --grow-mu 0) · inv18ib_<밑>/d<D>/mu0/ib_inv_s*_mu0[_d32]_gmu0.json · inv18D_<밑>/d<D>/mu0|mu0p003/dms_mat_<밑>_s*_list[_d32]_gmu0.json
  밀도 8 대조 = 같은 실험의 순수 공동체 밀도 8 (해부 14 · 15 의 진화 공동체 값과 섞지 않는다)
  G0   같은 시드 · 밑 코드의 dms μ0 · dms μ0.3% · invbud 의 grow_checksum 이 같다 · dms μ0 mut 팔 [t, n, m] = invbud 팔 · 배경 멸종 없음
  D1   (§4b 고정) 밑 코드마다 이웃 등급 ln R 차가 밀도 32 에서 밀도 8 보다 **작다**: 시드 짝 차(8 − 32) 하한 > 0 — 여섯 전부 · 비(32/8)는 서술
  D1b  밀도 32 코드 40 의 R 중앙값 ≤ R_MAX(1.4)   (밀도 8 의 R 중앙값도 적는다)
  D2a  대비(밀도 8 합친 Δμ − 밀도 32 합친 Δμ) 하한 > 0   ·   D2b 밀도 32 합친 Δμ: 구간이 0 을 품고 |평균| ≤ D2_TOL
  D3   밀도 32 μ 0 ΔΔ 60 짝 상한 < 0 · 합친 |ΔΔ(0)|₃₂ − |ΔΔ(0)|₈ 하한 > 0
  D4   서술: 밀도 32 같은-고리 |ln R 비| 최대 · N3(0.15) 실패 수
  D5   서술: 밀도 32 배경 상위 코드 · 개체 수
실행: python3 pairs18_analyze.py [petri]  (환경 PETRI_ROOT PETRI_NBASES PETRI_NBG PETRI_PILOT PETRI_REF14 PETRI_REF15 PETRI_DENS)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv18")
NBASES = os.environ.get("PETRI_NBASES", "racld rascld rsacld racldx acld rascled").split()
NBG = int(os.environ.get("PETRI_NBG", "20"))
PILOT = os.environ.get("PETRI_PILOT") == "1"
DENS = os.environ.get("PETRI_DENS", "32"); DREF = os.environ.get("PETRI_DREF", "8")
GT = os.environ.get("PETRI_GTAG", "_gmu0")
def suffix(d): return ("" if d == "8" else "_d" + d) + GT
D1_TOL, R_MAX, D2_TOL, TOL_NEG = 0.05, 1.4, 0.05, 0.15   # §4b 09-28 01:55 고정: D1 은 시드 짝 대비(8 − 32) 하한 > 0 · R_MAX 1.4(파일럿 중앙 1.27 · 최대 1.39) · D2a 대비 · D2b 절대 0.05
N_BOOT = 2000
RNG = random.Random(20261002)
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT, "dens": DENS}


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


def same_loop_pairs(codes):
    cl = sorted(codes); return [(a, b) for i, a in enumerate(cl) for b in cl[i + 1:] if codes[a] == codes[b]]


def read_ib(dirpath, suffix):
    """b -> code -> seed -> {R, T} · 팔 기록 · grow_checksum"""
    R, SER, CK = {}, {}, {}
    for f in sorted(glob.glob(os.path.join(dirpath, "ib_inv_s*_mu0%s.json" % suffix))):
        z = json.load(open(f)); s = z["seed"]
        CK[s] = z["grow_checksum"]
        for a in z["arms"]:
            if a["kind"] != "mut": continue
            SER.setdefault(s, {})[a["key"]] = ([r[:3] for r in a["series"]], a["stopped"])
            if a["stopped"] >= 0 and a["series"][-1][2] == 0: continue
            ok = sum(bn["copy_ok"][1] for bn in a["bins"]); st = sum(bn["copy_stall"][1] + bn["copy_del"][1] for bn in a["bins"]); ph = sum(bn["copy_phase"][1] for bn in a["bins"])
            if ok: R.setdefault(a["key"], {})[s] = {"R": (ok + st) / ok, "T": ph / ok}
    return R, SER, CK


def read_dms(dirpath, b, suffix):
    S, SER, CK, TOP = {}, {}, {}, {}
    for f in sorted(glob.glob(os.path.join(dirpath, "dms_mat_%s_s*_list%s.json" % (b, suffix)))):
        z = json.load(open(f)); s = z["seed"]
        if z["wt"] != b or not z["fidelity"]["same"]: bad.append("dms 설정 %s %d" % (dirpath, s)); continue
        if z["grow_extinct"] >= 0: bad.append("배경 멸종 %s %d" % (dirpath, s)); continue
        CK[s] = z["grow_checksum"]; TOP[s] = (z["grow_top"][0][0], z["grow_pop"])
        S[s] = {a["key"]: (s_of(a), a["stopped"] >= 0 and a["series"][-1][2] == 0) for a in z["arms"] if a["kind"] == "mut"}
        SER[s] = {a["key"]: ([r[:3] for r in a["series"]], a["stopped"]) for a in z["arms"] if a["kind"] == "mut"}
    return S, SER, CK, TOP


R32, R8, S32, S8 = {}, {}, {}, {}
g0 = {"ck": 0, "ck_bad": [], "ser": 0, "ser_bad": 0}
TOPS = {}; TOPS8 = {}
for b in NBASES:
    codes = inserts(b, "n")
    r32, ser_ib, ck_ib = read_ib(os.path.join(BASE, "%sib_%s" % (ROOT, b), "d" + DENS, "mu0"), suffix(DENS))
    r8, _, _ = read_ib(os.path.join(BASE, "%sib_%s" % (ROOT, b), "d" + DREF, "mu0"), suffix(DREF))
    R32[b], R8[b] = r32, r8
    if len(ck_ib) != NBG: bad.append("IB32 파일 수 %s %d/%d" % (b, len(ck_ib), NBG))
    for mu, tag in ((0.0, "mu0"), (0.003, "mu0p003")):
        s32, ser_d, ck_d, top = read_dms(os.path.join(BASE, "%sD_%s" % (ROOT, b), "d" + DENS, tag), b, suffix(DENS))
        s8, _, _, top8 = read_dms(os.path.join(BASE, "%sD_%s" % (ROOT, b), "d" + DREF, tag), b, suffix(DREF))
        if mu == 0.0: TOPS8[b] = top8
        S32[(b, mu)], S8[(b, mu)] = s32, s8
        if len(ck_d) != NBG: bad.append("D32 파일 수 %s %s %d/%d" % (b, mu, len(ck_d), NBG))
        for s in ck_d:
            if s in ck_ib:
                g0["ck"] += 1
                if ck_d[s] != ck_ib[s]: g0["ck_bad"].append((b, mu, s))
            if mu == 0.0 and s in ser_ib:
                for c in ser_d[s]:
                    if c in ser_ib[s]:
                        g0["ser"] += 1
                        if ser_d[s][c] != ser_ib[s][c]: g0["ser_bad"] += 1
        if mu == 0.0: TOPS[b] = top
if g0["ck_bad"] or g0["ser_bad"]: bad.append("G0 %s %d" % (g0["ck_bad"][:3], g0["ser_bad"]))
print("== 해부 18 %s— 밀도 %s × 짝 · 루트 %s · 밑 코드 %s · 배경 %d" % ("파일럿 " if PILOT else "판정 ", DENS, ROOT, " ".join(NBASES), NBG))
print("  G0 grow_checksum 대조 %d · 불일치 %d · dms=invbud 팔 %d · 불일치 %d" % (g0["ck"], len(g0["ck_bad"]), g0["ser"], g0["ser_bad"]))
print("  점검 문제 %d" % len(bad))
for x in bad[:10]: print("    ", x)
if bad and not PILOT: print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)

# ---- D5 서술 ----
print("\n== D5 서술 — 밀도 %s 배경의 상위 코드 · 개체 수" % DENS)
for b in NBASES:
    tops = TOPS.get(b, {})
    cnt = {}
    for s, (t, p) in tops.items(): cnt[t] = cnt.get(t, 0) + 1
    pops = [p for _, p in tops.values()]
    tops8 = TOPS8.get(b, {}); cnt8 = {}
    for s, (t, p) in tops8.items(): cnt8[t] = cnt8.get(t, 0) + 1
    pops8 = [p for _, p in tops8.values()]
    if pops: print("   %-8s 밀도 %s: 상위 코드 %s · 개체 수 중앙 %d [%d, %d]   |   밀도 %s: %s · 개체 %s" % (b, DENS, dict(sorted(cnt.items(), key=lambda kv: -kv[1])), pct(pops, .5), min(pops), max(pops), DREF, dict(sorted(cnt8.items(), key=lambda kv: -kv[1])), ("%d" % pct(pops8, .5)) if pops8 else "—"))
    OUT["D5|%s" % b] = {"tops": cnt, "pop_median": pct(pops, .5) if pops else None, "tops8": cnt8, "pop8_median": pct(pops8, .5) if pops8 else None}

# ---- D1 · D1b · D4 ----
print("\n== D1 — 이웃 등급 ln R 차(짧은 − 긴)가 밀도 %s 에서 밀도 %s 보다 작은가(시드 짝 차 하한 > 0) · D1b 밀도 %s R 중앙값 ≤ %.1f" % (DENS, DREF, DENS, R_MAX))
d1_all = True; allR = []
for b in NBASES:
    codes = inserts(b, "n"); Ls = sorted(set(codes.values()))
    seeds = sorted(set.intersection(*[set(R32[b].get(c, {})) for c in codes])) if all(R32[b].get(c) for c in codes) else []
    for c in codes:
        if R32[b].get(c): allR.append(sum(v["R"] for v in R32[b][c].values()) / len(R32[b][c]))
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        v = [sum(math.log(R32[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(R32[b][c][s]["R"]) for c in c2) / len(c2) for s in seeds]
        if not v: print("   %s: 자료 없음" % b); continue
        m, lo, hi = boot_mean(v)
        s8 = [s for s in seeds if all(s in R8[b].get(c, {}) for c in c1 + c2)]
        r8 = {s: sum(math.log(R8[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(R8[b][c][s]["R"]) for c in c2) / len(c2) for s in s8}
        v32 = {s: sum(math.log(R32[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(R32[b][c][s]["R"]) for c in c2) / len(c2) for s in seeds}
        dv = [r8[s] - v32[s] for s in s8 if s in v32]
        md, lod, hid = boot_mean(dv) if dv else (float("nan"),) * 3
        m8 = sum(r8.values()) / len(r8) if r8 else float("nan")
        ok = lod > 0; d1_all &= ok
        print("   %-8s ln R(%d) − ln R(%d): 밀도 %s %s · 밀도 %s %+.3f · 차(8 − 32) %s · 비 %.2f · 배경 %d %s" % (b, L1, L2, DENS, fmt(m, lo, hi), DREF, m8, fmt(md, lod, hid), (m / m8) if m8 else float("nan"), len(dv), "✅" if ok else "❌"))
        OUT["D1|%s|%d|%d" % (b, L1, L2)] = {"d32": m, "ci32": [lo, hi], "d8": m8, "diff": md, "ci_diff": [lod, hid], "ratio": (m / m8) if m8 else None, "pass": ok}
    mx = 0.0
    for a, c in same_loop_pairs(codes):
        sd = sorted(set(R32[b].get(a, {})) & set(R32[b].get(c, {})))
        if sd: mx = max(mx, abs(sum(math.log(R32[b][a][s]["R"]) - math.log(R32[b][c][s]["R"]) for s in sd) / len(sd)))
    OUT["D4|%s" % b] = {"max_same_lnR": mx}
    print("      D4 서술: 같은-고리 최대 |ln R 비| %.3f" % mx)
allR8 = [sum(v["R"] for v in R8[b][c].values()) / len(R8[b][c]) for b in NBASES for c in inserts(b, "n") if R8[b].get(c)]
Rmed8 = pct(allR8, .5) if allR8 else float("nan")
Rmed = pct(allR, .5) if allR else float("nan"); d1b = Rmed <= R_MAX
print("  밀도 %s R 중앙값 %.3f (범위 %.3f~%.3f)" % (DREF, Rmed8, min(allR8) if allR8 else float("nan"), max(allR8) if allR8 else float("nan")))
print("  D1 종합: %s · D1b: 코드 %d 의 R 중앙값 %.3f (범위 %.3f~%.3f) %s" % ("✅" if d1_all else "❌", len(allR), Rmed, min(allR) if allR else float("nan"), max(allR) if allR else float("nan"), "✅" if d1b else "❌"))
OUT["D1"] = d1_all; OUT["D1b"] = {"median": Rmed, "pass": d1b, "min": min(allR) if allR else None, "max": max(allR) if allR else None, "median8": Rmed8}   # 09-28 02:00 판정 전 출력 추가(그림용 · 판정선 무관)

# ---- D2 · D3 ----
def dd(Sd, b, mu, a, c):
    d = Sd.get((b, mu), {}); return {s: d[s][c][0] - d[s][a][0] for s in d if a in d[s] and c in d[s]}
print("\n== D2 · D3 — 밀도 %s 짝의 μ 0 ΔΔ 와 μ 상승 (대조 = 같은 순수 공동체 밀도 %s)" % (DENS, DREF))
d3_all, n3f = True, 0
per_seed = {"dmu32": {}, "dmu8": {}, "abs32": {}, "abs8": {}}
for b in NBASES:
    codes = inserts(b, "n")
    for a, c in pairs_of(codes):
        d0 = dd(S32, b, 0.0, a, c); d3 = dd(S32, b, 0.003, a, c); e0 = dd(S8, b, 0.0, a, c); e3 = dd(S8, b, 0.003, a, c)
        if d0:
            m, lo, hi = boot_mean([d0[s] for s in sorted(d0)]); ok = hi < 0; d3_all &= ok
            OUT["D3|%s|%s|%s" % (b, a, c)] = {"dd0": m, "ci": [lo, hi], "pass": ok}
            if not ok: print("   ❌ D3 %s %s→%s ΔΔ(0) %s" % (b, a, c, fmt(m, lo, hi)))
        for s in set(d0) & set(d3): per_seed["dmu32"].setdefault(s, []).append(d3[s] - d0[s])
        for s in set(e0) & set(e3): per_seed["dmu8"].setdefault(s, []).append(e3[s] - e0[s])
        for s in d0: per_seed["abs32"].setdefault(s, []).append(abs(d0[s]))
        for s in e0: per_seed["abs8"].setdefault(s, []).append(abs(e0[s]))
    for a, c in same_loop_pairs(codes):
        v = dd(S32, b, 0.0, a, c)
        if v and abs(sum(v.values()) / len(v)) > TOL_NEG: n3f += 1
def pooled(key):
    d = per_seed[key]; v = [sum(x) / len(x) for s, x in sorted(d.items()) if x]
    return (boot_mean(v) if v else (float("nan"),) * 3), len(v)
(m32, l32, h32), n32 = pooled("dmu32"); (m8, l8, h8), n8 = pooled("dmu8")
cs = sorted(set(per_seed["dmu32"]) & set(per_seed["dmu8"]))
(mc, lc, hc) = boot_mean([sum(per_seed["dmu8"][s]) / len(per_seed["dmu8"][s]) - sum(per_seed["dmu32"][s]) / len(per_seed["dmu32"][s]) for s in cs]) if cs else (float("nan"),) * 3
d2a = lc > 0; d2b = (l32 <= 0 <= h32) and abs(m32) <= D2_TOL; d2 = d2a
print("   D2a 대비(밀도 %s Δμ − 밀도 %s Δμ) %s %s   ·   D2b 밀도 %s 합친 Δμ %s (배경 %d) %s   ·   밀도 %s 합친 Δμ %s" % (DREF, DENS, fmt(mc, lc, hc), "✅" if d2a else "❌", DENS, fmt(m32, l32, h32), n32, "✅" if d2b else "❌", DREF, fmt(m8, l8, h8)))
(a32, al32, ah32), _ = pooled("abs32"); (a8, al8, ah8), _ = pooled("abs8")
cs2 = sorted(set(per_seed["abs32"]) & set(per_seed["abs8"]))
(mk, lk, hk) = boot_mean([sum(per_seed["abs32"][s]) / len(per_seed["abs32"][s]) - sum(per_seed["abs8"][s]) / len(per_seed["abs8"][s]) for s in cs2]) if cs2 else (float("nan"),) * 3
npairs = sum(1 for k in OUT if k.startswith("D3|")); npass = sum(1 for k, v in OUT.items() if k.startswith("D3|") and v["pass"])
d3 = d3_all and lk > 0
print("   D3 μ 0 ΔΔ 상한 < 0: %d/%d 짝 %s · 합친 |ΔΔ(0)| 밀도 %s %s 대 밀도 8 %s · 차(%s − 8) %s %s" % (npass, npairs, "✅" if d3_all else "❌", DENS, fmt(a32, al32, ah32), fmt(a8, al8, ah8), DENS, fmt(mk, lk, hk), "✅" if lk > 0 else "❌"))
print("   D4 서술: 밀도 %s N3(0.15) 실패 %d 쌍" % (DENS, n3f))
OUT.update({"D2": {"m32": m32, "ci32": [l32, h32], "m8": m8, "ci8": [l8, h8], "contrast": mc, "ci_c": [lc, hc], "pass_a": d2a, "pass_b": d2b, "pass": d2},
            "D3": {"npass": npass, "npairs": npairs, "abs32": a32, "abs8": a8, "diff": mk, "ci_diff": [lk, hk], "pass": d3}, "D4_n3fail": n3f})
print("\n== 종합: G0 %s · D1 %s · D1b %s · D2a %s · D2b %s · D3a %s · D3b %s" % ("✅" if not (g0["ck_bad"] or g0["ser_bad"]) else "❌", "✅" if d1_all else "❌", "✅" if d1b else "❌", "✅" if d2a else "❌", "✅" if d2b else "❌", "✅" if d3_all else "❌", "✅" if lk > 0 else "❌"))
OUT["bad"] = bad
fn = "_RESULT_pilot18.json" if PILOT else "_RESULT_pairs18.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
