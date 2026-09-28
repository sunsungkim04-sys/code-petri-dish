#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 16 판정 — 긴 꼬리(고리가 코드의 작은 몫일 때). 사전등록 해부16-사전등록-2026-09-27.md §4 — 문턱은 해부 13 · 14 것 그대로.
  G0   dms 목록 모드 mut 팔 [t, n, m] = invbud mut 팔 [t, n, m] (같은 배경 · 코드 · 난수) · 불일치 0
  자기 검산  T/R − L ∈ [−0.5, 1.5]
  F1   밑 코드마다 ln R(4) − ln R(5) 하한 > 0                       (해부 14 A1)
  F2   같은-고리 쌍 |ln R 비| ≤ 0.10 · |T/R 차| ≤ 0.10 틱             (해부 14 A2a · A2b)
  F3   모든 4→5 짝 ΔΔ 상한 < 0                                       (해부 13 N1)
  F4   서술: 틱당 ΔΔ 중앙값 · (R̄ − 1) × ℓ · N3(0.15) 실패 쌍 · 사라진 팔
실행: python3 pairs16_analyze.py [petri 폴더]  (환경: PETRI_ROOT PETRI_NBASES PETRI_NBG PETRI_PILOT)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv16")
NBASES = os.environ.get("PETRI_NBASES", "racldxx racldxxxx").split()   # 꼬리 8 은 파일럿에서 계통 소멸로 제외(사전등록 §2)
NBG = int(os.environ.get("PETRI_NBG", "20"))
PILOT = os.environ.get("PETRI_PILOT") == "1"
TAU, TR_PAIR, TR_LO, TR_HI, TOL_NEG = 0.10, 0.10, -0.5, 1.5, 0.15
N_BOOT = 2000
RNG = random.Random(20260930)
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT}


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


KEYS = ["copy_ok", "copy_stall", "copy_del", "copy_phase"]
S, IB, DEAD, g0 = {}, {}, {}, {"checked": 0, "mismatch": 0}
for b in NBASES:
    codes = inserts(b, "n"); S[b] = {}; IB[b] = {c: {} for c in codes}; DEAD[b] = {c: [] for c in codes}
    ser_d = {}
    fs = sorted(glob.glob(os.path.join(BASE, "%sD_%s" % (ROOT, b), "mu0", "dms_mat_%s_s*_list.json" % b)))
    if len(fs) != NBG: bad.append("D 파일 수 %s %d/%d" % (b, len(fs), NBG))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if z["wt"] != b or abs(z["mu_assay"]) > 1e-12 or not z["checksum_match"] or not z["fidelity"]["same"]: bad.append("D 설정 %s %d" % (b, s)); continue
        muts = [a for a in z["arms"] if a["kind"] == "mut"]
        if [a["key"] for a in muts] != list(codes): bad.append("D 코드 목록 %s %d" % (b, s))
        S[b][s] = {}
        for a in muts:
            dead = a["stopped"] >= 0 and a["series"][-1][2] == 0
            if dead: DEAD[b][a["key"]].append(s)
            S[b][s][a["key"]] = (s_of(a), dead)
            ser_d[(s, a["key"])] = ([r[:3] for r in a["series"]], a["stopped"])
    fs = sorted(glob.glob(os.path.join(BASE, "%sib_%s" % (ROOT, b), "mu0", "ib_inv_s*_mu0.json")))
    if len(fs) != NBG: bad.append("IB 파일 수 %s %d/%d" % (b, len(fs), NBG))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if not z["checksum_match"] or list(z["codes"]) != list(codes): bad.append("IB 설정 %s %d" % (b, s))
        for a in z["arms"]:
            if a["kind"] != "mut": continue
            key = (s, a["key"])
            if key in ser_d:
                g0["checked"] += 1
                if ser_d[key] != ([r[:3] for r in a["series"]], a["stopped"]): g0["mismatch"] += 1
            if a["stopped"] >= 0 and a["series"][-1][2] == 0: continue
            d = {k: sum(bn[k][1] for bn in a["bins"]) for k in KEYS}
            if d["copy_del"] != 0: bad.append("IB μ0 빠뜨림 %s %d %s" % (b, s, a["key"]))
            if d["copy_ok"] == 0: continue
            IB[b][a["key"]][s] = {"R": (d["copy_ok"] + d["copy_stall"]) / d["copy_ok"], "T": d["copy_phase"] / d["copy_ok"]}
if g0["mismatch"]: bad.append("G0 dms≠invbud %d/%d" % (g0["mismatch"], g0["checked"]))
print("== 해부 16 %s— 루트 %s · 밑 코드 %s · 배경 %d" % ("파일럿 " if PILOT else "판정 ", ROOT, " ".join(NBASES), NBG))
print("  G0 dms 팔 = invbud 팔: 대조 %d · 불일치 %d" % (g0["checked"], g0["mismatch"]))
excluded = []
for b in DEAD:
    for c, v in DEAD[b].items():
        if v:
            print("  사라진 팔 %s %s: 배경 %s%s" % (b, c, v, " → 코드 제외" if len(v) * 2 > NBG else ""))
            if len(v) * 2 > NBG: excluded.append((b, c))
print("  점검 문제 %d" % len(bad))
for x in bad[:10]: print("    ", x)
if bad and not PILOT: print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)

f1_all, f2_all, f3_all, tr_ok = True, True, True, True
for b in NBASES:
    codes = {c: L for c, L in inserts(b, "n").items() if (b, c) not in excluded}
    print("\n-- 밑 코드 `%s` (길이 %d · 코드 %d)" % (b, len(b), len(codes)))
    pooled = {}
    for c in codes:
        ds = list(IB[b][c].values())
        if not ds: continue
        ok = sum(d["R"] for d in ds) / len(ds); pooled[c] = (ok, sum(d["T"] for d in ds) / len(ds))
        dev = pooled[c][1] / pooled[c][0] - codes[c]
        if not (TR_LO <= dev <= TR_HI): tr_ok = False; print("   ❌ 자기 검산 %s T/R−L %+.2f" % (c, dev))
    Ls = sorted(set(codes.values()))
    seeds = sorted(set.intersection(*[set(IB[b][c]) for c in codes if IB[b][c]])) if any(IB[b][c] for c in codes) else []
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        v = [sum(math.log(IB[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(IB[b][c][s]["R"]) for c in c2) / len(c2) for s in seeds]
        if not v: continue
        m, lo, hi = boot_mean(v); ok = lo > 0; f1_all &= ok
        print("   F1 ln R(%d) − ln R(%d) = %s · 배경 %d %s" % (L1, L2, fmt(m, lo, hi), len(v), "✅" if ok else "❌"))
        OUT["F1|%s|%d|%d" % (b, L1, L2)] = {"d": m, "ci": [lo, hi], "pass": ok}
    same = [(a, c) for i, a in enumerate(sorted(codes)) for c in sorted(codes)[i + 1:] if codes[a] == codes[c]]
    mx = 0.0; mxt = 0.0
    for a, c in same:
        sd = sorted(set(IB[b][a]) & set(IB[b][c]))
        if not sd or a not in pooled or c not in pooled: continue
        m = sum(math.log(IB[b][a][s]["R"]) - math.log(IB[b][c][s]["R"]) for s in sd) / len(sd)
        trd = abs(pooled[a][1] / pooled[a][0] - pooled[c][1] / pooled[c][0])
        mx = max(mx, abs(m)); mxt = max(mxt, trd)
        if abs(m) > TAU or trd > TR_PAIR: f2_all = False; print("   ❌ F2 %s vs %s ln R 비 %+.3f · T/R 차 %.3f" % (a, c, m, trd))
    print("   F2 같은-고리 %d 쌍 · 최대 |ln R 비| %.3f · 최대 T/R 차 %.3f %s" % (len(same), mx, mxt, "✅" if (mx <= TAU and mxt <= TR_PAIR) else "❌"))
    OUT["F2|%s" % b] = {"n": len(same), "max_lnR": mx, "max_tr": mxt}
    # F3 · F4
    pairs = [(a, c) for a in codes for c in codes if codes[a] == 4 and codes[c] == 5]
    dds, n3f = [], 0
    for a, c in pairs:
        v = [S[b][s][c][0] - S[b][s][a][0] for s in S[b] if a in S[b][s] and c in S[b][s]]
        if not v: continue
        m, lo, hi = boot_mean(v); ok = hi < 0; f3_all &= ok; dds.append(m)
        OUT["F3|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "pass": ok, "n_bg": len(v)}
        if not ok: print("   ❌ F3 %s→%s ΔΔ %s" % (a, c, fmt(m, lo, hi)))
    for a, c in same:
        v = [S[b][s][c][0] - S[b][s][a][0] for s in S[b] if a in S[b][s] and c in S[b][s]]
        if v and abs(sum(v) / len(v)) > TOL_NEG: n3f += 1
    Rbar = sum(p[0] for p in pooled.values()) / len(pooled) if pooled else float("nan")
    print("   F3 4→5 짝 %d · 상한<0 %d %s · F4 서술: 틱당 ΔΔ 중앙값 %s · 범위 %s ~ %s · (R̄−1)·ℓ %.1f · N3(0.15) 실패 %d/%d 쌍" % (
        len(pairs), sum(1 for k, v in OUT.items() if k.startswith("F3|%s|" % b) and v["pass"]), "✅" if all(v["pass"] for k, v in OUT.items() if k.startswith("F3|%s|" % b)) else "❌",
        ("%+.3f" % pct(dds, .5)) if dds else "—", ("%+.3f" % min(dds)) if dds else "—", ("%+.3f" % max(dds)) if dds else "—", (Rbar - 1) * (len(b) + 1), n3f, len(same)))
    OUT["F4|%s" % b] = {"dd_median": pct(dds, .5) if dds else None, "dd_min": min(dds) if dds else None, "dd_max": max(dds) if dds else None, "Rbar": Rbar, "n3_fail": n3f, "n_same": len(same)}
print("\n== 종합: G0 %s · 자기 검산 %s · F1 %s · F2 %s · F3 %s" % ("✅" if not g0["mismatch"] else "❌", "✅" if tr_ok else "❌", "✅" if f1_all else "❌", "✅" if f2_all else "❌", "✅" if f3_all else "❌"))
OUT.update({"G0": g0, "tr": tr_ok, "F1": f1_all, "F2": f2_all, "F3": f3_all, "excluded": excluded, "bad": bad})
fn = "_RESULT_pilot16.json" if PILOT else "_RESULT_pairs16.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
