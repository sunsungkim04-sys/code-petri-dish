#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 14 판정 — 짝 코드 전수(밑 코드 여섯의 n 삽입 40 코드)의 글자당 주사위 횟수를 직접 센다
사전등록: 해부14-사전등록-2026-09-27.md §1 · §4. 판정선(§4b)은 파일럿 뒤에 고정한다 — 고정 전에는 PETRI_PILOT=1 로만 돈다.

  IB (invbud.js · μ 0)  끼운 계통(표지 1)의 팔 3,000틱 합:
     R = (copy_ok + copy_stall + copy_del) / copy_ok      글자 하나에 주사위 몇 번  (μ 0 에서 copy_del = 0 — 점검)
     T = copy_phase / copy_ok                             글자 하나에 몇 틱
     S = (R − 1) × L                                      글자당 헛손질 틱   (L = loop_ticks — 해부 13 판정기와 같은 함수)
  Q  (dms.js --wt <밑 코드> --arms list --kids 1 · μ 0.3%)  끼운 코드 그대로인 부모의 출생 b_exact 와 자식 분포(ref 팔은 밑 코드 — 안 쓴다):
     U_copy = Σ kids / b_exact  ·  U_all = (Σ kids + Σ cosmic) / b_exact (해부 12 정의)

  G0a  IB 모든 팔의 [t, n, m] 기록이 같은 배경 · 코드의 μ 0 이웃 팔(inv12nbr_* · inv13nbr_*) 기록과 같다
  G0b  Q 팔 기록(0~3열)이 같은 배경 · 코드 · μ 0.3% 이웃 팔(inv13nbr_*)과 같다(둘째 가족만 — 첫 가족은 μ 0.3% 이웃 기록이 없다) · 4열은 --wt 에 따라 다른 장부라 뺀다(09-27 개정 1회)
  자기 검산  모든 코드 T/R − L ∈ [TR_LO, TR_HI]  — 벗어나면 멈춤(고리 배정 의심)
  A1   밑 코드마다 이웃 고리 등급 (L, L+1): 배경마다 등급 평균 ln R 의 차 ln R_L − ln R_{L+1} — 하한 > 0
  A1b  서술: 등급 평균 S 의 max/min
  A2a  같은 고리의 두 배치 |ln(R_i/R_j)| ≤ TAU(0.10) — 모든 같은-고리 쌍(N3 실패 아홉 쌍 표시)
  A2b  같은 고리의 두 배치 |T/R_i − T/R_j| ≤ TR_PAIR(0.10 틱) — 자리는 시도당 틱(= 고리 길이)을 못 바꾼다
  A2c  서술: 밑 코드마다 max 같은-고리 |ln R 비| / min 이웃 등급 효과
  A3   서술: N3 실패 아홉 쌍의 R · T · (ln T 비 − ln R 비) · 출생/노출 · 사망/노출(원인별) ln 비
  A4   밑 코드마다 이웃 등급의 U_copy: 배경마다 등급 평균 ln U 의 차 — 하한 > 0
  A4b  서술: 40 코드에서 ln U_copy 대 ln R 의 기울기 · Spearman
실행: python3 pairs14_analyze.py [petri 폴더]   (환경: PETRI_ROOT · PETRI_NBASES · PETRI_NBG · PETRI_REF1 · PETRI_REF2 · PETRI_PILOT)
  → _RESULT_pairs14.json (파일럿이면 _RESULT_pilot14.json)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv14")
NBASES = os.environ.get("PETRI_NBASES", "racld rascld rsacld racldx acld rascled").split()
NBG = int(os.environ.get("PETRI_NBG", "20"))
MUQ = float(os.environ.get("PETRI_MUQ", "0.003"))
REF1 = os.environ.get("PETRI_REF1", "inv12")     # 첫 가족 μ 0 이웃 기록의 루트(빈 문자열이면 기준 없음)
REF2 = os.environ.get("PETRI_REF2", "inv13")     # 둘째 가족 μ 0 · μ 0.3% 이웃 기록의 루트
PILOT = os.environ.get("PETRI_PILOT") == "1"
FAM1 = {"racld", "rascld", "rsacld", "racldx"}
# ---- 판정선 (§4b 에서 고정 · 고정 전 None) ----
TAU = 0.10            # A2a: |ln(R_i/R_j)| ≤ TAU — §4b 09-27 파일럿(시드 24 · 25) 뒤 고정: 파일럿 최소 고리 효과 0.203 의 절반. 파일럿 같은-고리 최대 0.052
TR_PAIR = 0.10        # A2b: |T/R_i − T/R_j| ≤ TR_PAIR 틱 — 파일럿 등급 안 퍼짐 ≤ 0.05
TR_LO, TR_HI = -0.5, 1.5
S_BAND = 1.25         # A1b 서술용
N_BOOT = 2000
RNG = random.Random(20260928)
N3_FAIL = {("acnld", "ancld"), ("acnld", "nacld"), ("nrascled", "rasclend"), ("ranscled", "rasclend"), ("rascledn", "rasclend"),
           ("rasclend", "rasclned"), ("rasclend", "rnascled"), ("rascldn", "rnascld"), ("rasclnd", "rnascld")}
N3_FAIL |= {(b, a) for a, b in N3_FAIL}
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT, "tau": TAU}


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l")
    si = code.rfind("s", 0, li)
    start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]
        L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


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


def ref_root(b):
    return REF1 if b in FAM1 else REF2


KEYS_SUM = ["copy_ok", "copy_stall", "copy_del", "copy_phase", "births", "deaths", "death_age", "death_prestep", "death_other",
            "stepped", "not_stepped", "unknown"]

# ---------------- IB 읽기 ----------------
IB = {}          # b -> code -> seed -> dict(sums + R, T)
DEAD = {}        # b -> code -> [seed]
g0a = {"checked": 0, "mismatch": [], "no_ref": 0}
for b in NBASES:
    codes = inserts(b, "n")
    IB[b] = {c: {} for c in codes}; DEAD[b] = {c: [] for c in codes}
    fs = sorted(glob.glob(os.path.join(BASE, "%sib_%s" % (ROOT, b), "mu0", "ib_inv_s*_mu0.json")))
    if len(fs) != NBG:
        bad.append("IB 파일 수 %s %d/%d" % (b, len(fs), NBG))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if z["mode"] != "inv" or abs(z["mu"]) > 1e-12: bad.append("IB 설정 %s %d" % (b, s))
        if not z["checksum_match"]: bad.append("IB 배경 체크섬 %s %d" % (b, s))
        if list(z["codes"]) != list(codes): bad.append("IB 코드 목록 %s %d" % (b, s))
        if len(z["arms"]) != 1 + len(codes): bad.append("IB 팔 수 %s %d (%d)" % (b, s, len(z["arms"])))
        rr = ref_root(b)
        ref = None
        if rr:
            rp = os.path.join(BASE, "%snbr_%s" % (rr, b), "mu0", "dms_mat_racld_s%05d_nbr_%s.json" % (s, b))
            if os.path.exists(rp):
                ref = json.load(open(rp))
        if ref is None:
            g0a["no_ref"] += 1
        for a in z["arms"]:
            if ref is not None:
                cand = [x for x in ref["arms"] if x["key"] == a["key"] and (x["kind"] == "mut") == (a["kind"] == "mut")]
                cand = [x for x in cand if x["kind"] != "neutral"]
                if not cand:
                    g0a["mismatch"].append("%s s%d %s: 기준 팔 없음" % (b, s, a["key"]))
                else:
                    x = cand[0]
                    g0a["checked"] += 1
                    if [r[:3] for r in x["series"]] != [r[:3] for r in a["series"]] or x["stopped"] != a["stopped"]:
                        g0a["mismatch"].append("%s s%d %s" % (b, s, a["key"]))
            if a["kind"] != "mut":
                continue
            c = a["key"]
            d = {k: sum(bn[k][1] for bn in a["bins"]) for k in KEYS_SUM}
            if a["stopped"] >= 0 and a["series"][-1][2] == 0:
                DEAD[b][c].append(s)
                continue
            if d["copy_del"] != 0: bad.append("IB μ0 인데 빠뜨림 %s %d %s" % (b, s, c))
            expo = d["stepped"] + d["not_stepped"]
            if expo and d["unknown"] / expo > 0.01: bad.append("IB 모름 몫 %s %d %s %.3f" % (b, s, c, d["unknown"] / expo))
            if d["copy_ok"] == 0: bad.append("IB copy_ok 0 %s %d %s" % (b, s, c)); continue
            d["R"] = (d["copy_ok"] + d["copy_stall"] + d["copy_del"]) / d["copy_ok"]
            d["T"] = d["copy_phase"] / d["copy_ok"]
            d["expo"] = expo
            IB[b][c][s] = d

# ---------------- Q 읽기 ----------------
Q = {}           # b -> code -> seed -> dict(b_exact, kids, cosmic, U_copy, U_all)
g0b = {"checked": 0, "mismatch": [], "no_ref": 0}
for b in NBASES:
    codes = inserts(b, "n")
    Q[b] = {c: {} for c in codes}
    fs = sorted(glob.glob(os.path.join(BASE, "%skids_%s" % (ROOT, b), mtag(MUQ), "dms_mat_%s_s*_list.json" % b)))
    if len(fs) != NBG:
        bad.append("Q 파일 수 %s %d/%d" % (b, len(fs), NBG))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if abs(z["mu_assay"] - MUQ) > 1e-12 or not z.get("kids_on") or z["wt"] != b: bad.append("Q 설정 %s %d" % (b, s))
        if not z["checksum_match"]: bad.append("Q 배경 체크섬 %s %d" % (b, s))
        if not z["fidelity"]["same"]: bad.append("Q 충실도 %s %d" % (b, s))
        muts = [a for a in z["arms"] if a["kind"] == "mut"]
        if [a["key"] for a in muts] != list(codes): bad.append("Q 코드 목록 %s %d" % (b, s))
        if len(z["arms"]) != 11 + len(codes): bad.append("Q 팔 수 %s %d (%d)" % (b, s, len(z["arms"])))
        ref = None
        if b not in FAM1 and REF2:
            rp = os.path.join(BASE, "%snbr_%s" % (REF2, b), mtag(MUQ), "dms_mat_racld_s%05d_nbr_%s.json" % (s, b))
            if os.path.exists(rp):
                ref = json.load(open(rp))
        if ref is None:
            g0b["no_ref"] += 1
        for a in muts:
            if ref is not None:
                cand = [x for x in ref["arms"] if x["kind"] == "mut" and x["key"] == a["key"]]
                if not cand:
                    g0b["mismatch"].append("%s s%d %s: 기준 팔 없음" % (b, s, a["key"]))
                else:
                    g0b["checked"] += 1
                    # 09-27 판정 전 점검 개정 1회: --exact 1 의 다섯째 칸(나머지 중 WT 그대로)은 --wt 가 다르면 당연히 다르다(실측: 0~3열 260/260 같음 · 4열만 다름) → 0~3열만 대조
                    if [r[:4] for r in cand[0]["series"]] != [r[:4] for r in a["series"]] or cand[0]["stopped"] != a["stopped"]:
                        g0b["mismatch"].append("%s s%d %s" % (b, s, a["key"]))
            kd = a["kids"]
            tot = kd["b_exact"] + kd["b_other"] + kd["unknown"]
            if tot and kd["unknown"] / tot > 0.01: bad.append("Q 모름 몫 %s %d %s" % (b, s, a["key"]))
            if kd["b_exact"] == 0:
                continue
            nk = sum(n for _, n in kd["kids"]); nc = sum(n for _, n in kd["cosmic"])
            Q[b][a["key"]][s] = {"b_exact": kd["b_exact"], "kids": nk, "cosmic": nc,
                                 "U_copy": nk / kd["b_exact"], "U_all": (nk + nc) / kd["b_exact"]}

# ---------------- 점검 ----------------
print("== 해부 14 %s— 짝 코드 전수의 글자당 주사위 횟수 · 루트 %s · 밑 코드 %s · 배경 %d" % ("파일럿 " if PILOT else "판정 ", ROOT, " ".join(NBASES), NBG))
print("  G0a IB 회귀: 대조 %d 팔 · 불일치 %d · 기준 없는 파일 %d" % (g0a["checked"], len(g0a["mismatch"]), g0a["no_ref"]))
for x in g0a["mismatch"][:5]: print("     ✗", x)
print("  G0b Q 회귀: 대조 %d 팔 · 불일치 %d · 기준 없는 파일 %d" % (g0b["checked"], len(g0b["mismatch"]), g0b["no_ref"]))
for x in g0b["mismatch"][:5]: print("     ✗", x)
if g0a["mismatch"] or g0b["mismatch"]:
    bad.append("G0 회귀 불일치 %d + %d" % (len(g0a["mismatch"]), len(g0b["mismatch"])))
if not PILOT and (g0a["no_ref"] or (g0b["no_ref"] and any(b not in FAM1 for b in NBASES))):
    notes.append("기준 기록 없는 파일 G0a %d · G0b %d" % (g0a["no_ref"], g0b["no_ref"]))
ndead = sum(len(v) for b in DEAD for v in DEAD[b].values())
print("  끼운 계통이 사라진 팔 %d" % ndead)
for b in DEAD:
    for c, v in DEAD[b].items():
        if v: print("     %s %s: 배경 %s" % (b, c, v))
print("  점검 문제 %d" % len(bad))
for x in bad[:12]: print("    ", x)
if bad and not PILOT:
    print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)
if TAU is None and not PILOT:
    print("🚨 TAU 가 고정되지 않았다(§4b) — 판정하지 않는다"); sys.exit(1)


def pooled(b, c):
    ds = IB[b][c].values()
    ok = sum(d["copy_ok"] for d in ds); st = sum(d["copy_stall"] for d in ds); de = sum(d["copy_del"] for d in ds); ph = sum(d["copy_phase"] for d in ds)
    return ((ok + st + de) / ok if ok else float("nan"), ph / ok if ok else float("nan"))


# ---------------- 자기 검산 T/R − L ----------------
print("\n== 계측기 자기 검산 — T/R − L ∈ [%.1f, %.1f]" % (TR_LO, TR_HI))
tr_fail = []
for b in NBASES:
    codes = inserts(b, "n")
    for c, L in codes.items():
        if not IB[b][c]: continue
        R, T = pooled(b, c)
        dev = T / R - L
        flag = "✅" if TR_LO <= dev <= TR_HI else "❌"
        if flag == "❌": tr_fail.append((b, c, dev))
        print("   %-9s L %d · R %.3f · T %.2f · T/R %.2f · 차 %+.2f %s" % (c, L, R, T, T / R, dev, flag))
        OUT["code|%s|%s" % (b, c)] = {"L": L, "R": R, "T": T, "S": (R - 1) * L, "n_bg": len(IB[b][c])}
if tr_fail:
    print("🚨 자기 검산 실패 %d 코드 — 고리 배정 의심 · 멈춤" % len(tr_fail))
    if not PILOT: sys.exit(1)

# ---------------- A1 · A1b ----------------
print("\n== A1 — 이웃 고리 등급의 배경별 등급 평균 ln R 차 (L 짧은 쪽 − 긴 쪽) · 하한 > 0")
a1_all = True
for b in NBASES:
    codes = inserts(b, "n")
    Ls = sorted(set(codes.values()))
    seeds = sorted(set.intersection(*[set(IB[b][c]) for c in codes if IB[b][c]])) if any(IB[b][c] for c in codes) else []
    print("  밑 코드 `%s` · 등급 %s · 공통 배경 %d" % (b, Ls, len(seeds)))
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        v = []
        for s in seeds:
            m1 = sum(math.log(IB[b][c][s]["R"]) for c in c1) / len(c1)
            m2 = sum(math.log(IB[b][c][s]["R"]) for c in c2) / len(c2)
            v.append(m1 - m2)
        if not v: continue
        m, lo, hi = boot_mean(v)
        ok = lo > 0; a1_all &= ok
        print("     ln R(%d) − ln R(%d) = %s %s" % (L1, L2, fmt(m, lo, hi), "✅" if ok else "❌"))
        OUT["A1|%s|%d|%d" % (b, L1, L2)] = {"d": m, "ci": [lo, hi], "pass": ok, "n_bg": len(v)}
    # A1b
    Sc = {}
    for L in Ls:
        vals = [OUT["code|%s|%s" % (b, c)]["S"] for c in codes if codes[c] == L and ("code|%s|%s" % (b, c)) in OUT]
        if vals: Sc[L] = sum(vals) / len(vals)
    if Sc:
        ratio = max(Sc.values()) / min(Sc.values())
        print("     A1b 서술: 등급 평균 S = %s · max/min %.3f %s" % (" · ".join("L%d %.2f" % (L, v) for L, v in sorted(Sc.items())), ratio, "(≤ %.2f)" % S_BAND if ratio <= S_BAND else "(> %.2f)" % S_BAND))
        OUT["A1b|%s" % b] = {"S": Sc, "ratio": ratio}
print("  A1 종합: %s" % ("✅ 모든 밑 코드 · 모든 이웃 등급" if a1_all else "❌"))
OUT["A1"] = a1_all

# ---------------- A2 · A3 ----------------
print("\n== A2 — 같은 고리의 두 배치: A2a |ln(R_i/R_j)| ≤ %s · A2b |T/R_i − T/R_j| ≤ %.2f 틱 · (★ = 해부 13 N3 실패 쌍)" % (TAU if TAU is not None else "TAU 미고정", TR_PAIR))
a2_all, a2_n, a2_fail, a2_abs = True, 0, 0, []
a2b_all, a2b_fail, a2b_max = True, 0, 0.0
A2C = {}
for b in NBASES:
    codes = inserts(b, "n")
    cl = sorted(codes)
    for i, ci in enumerate(cl):
        for cj in cl[i + 1:]:
            if codes[ci] != codes[cj]: continue
            seeds = sorted(set(IB[b][ci]) & set(IB[b][cj]))
            if not seeds: continue
            v = [math.log(IB[b][ci][s]["R"]) - math.log(IB[b][cj][s]["R"]) for s in seeds]
            m, lo, hi = boot_mean(v)
            star = "★" if (ci, cj) in N3_FAIL else " "
            a2_n += 1; a2_abs.append(abs(m))
            A2C.setdefault(b, []).append(abs(m))
            ki, kj = OUT["code|%s|%s" % (b, ci)], OUT["code|%s|%s" % (b, cj)]
            trd = abs(ki["T"] / ki["R"] - kj["T"] / kj["R"])
            okb = trd <= TR_PAIR; a2b_all &= okb; a2b_max = max(a2b_max, trd)
            if not okb: a2b_fail += 1
            if TAU is not None:
                ok = abs(m) <= TAU; a2_all &= ok
                if not ok: a2_fail += 1
                fl = "✅" if ok else "❌"
            else:
                fl = ""
            print("   %s %-9s vs %-9s (L %d) ln R 비 = %s · 배경 %d %s · T/R 차 %.3f %s" % (star, ci, cj, codes[ci], fmt(m, lo, hi), len(v), fl, trd, "✅" if okb else "❌"))
            OUT["A2|%s|%s|%s" % (b, ci, cj)] = {"d": m, "ci": [lo, hi], "n3fail": (ci, cj) in N3_FAIL, "n_bg": len(v), "pass": (abs(m) <= TAU) if TAU is not None else None, "tr_diff": trd, "pass_b": okb}
if a2_abs:
    print("  같은-고리 쌍 %d · |ln 비| 중앙값 %.4f · 최대 %.4f · 90분위 %.4f" % (a2_n, pct(a2_abs, 0.5), max(a2_abs), pct(a2_abs, 0.9)))
    OUT["A2_abs"] = {"n": a2_n, "median": pct(a2_abs, 0.5), "max": max(a2_abs), "p90": pct(a2_abs, 0.9)}
if TAU is not None:
    print("  A2a 종합: %s (%d/%d 쌍)" % ("✅" if a2_all else "❌", a2_n - a2_fail, a2_n))
    OUT["A2"] = a2_all
print("  A2b 종합: %s (%d/%d 쌍 · 최대 T/R 차 %.3f 틱)" % ("✅" if a2b_all else "❌", a2_n - a2b_fail, a2_n, a2b_max))
OUT["A2b"] = a2b_all; OUT["A2b_max"] = a2b_max
for b in NBASES:
    eff = [v["d"] for k, v in OUT.items() if k.startswith("A1|%s|" % b)]
    if eff and A2C.get(b):
        print("  A2c 서술 `%s`: max 같은-고리 |ln R 비| %.3f / min 이웃 등급 효과 %.3f = %.2f" % (b, max(A2C[b]), min(eff), max(A2C[b]) / min(eff) if min(eff) else float("nan")))
        OUT["A2c|%s" % b] = {"max_same": max(A2C[b]), "min_eff": min(eff)}

print("\n== A3 서술 — N3 실패 쌍의 ln 비 (앞 코드 / 뒤 코드): R · T · (T−R) · 출생/노출 · 사망/노출 · 나이 · 차례 전 · 그 밖")
for b in NBASES:
    codes = inserts(b, "n")
    for (ci, cj) in sorted(p for p in N3_FAIL if p[0] < p[1] and p[0] in codes and p[1] in codes):
        seeds = sorted(set(IB[b][ci]) & set(IB[b][cj]))
        if not seeds: continue
        row = {}
        for name, f in [("R", lambda d: d["R"]), ("T", lambda d: d["T"]), ("T-R", lambda d: d["T"] / d["R"]), ("births", lambda d: d["births"] / d["expo"]), ("deaths", lambda d: d["deaths"] / d["expo"]),
                        ("age", lambda d: d["death_age"] / d["expo"]), ("prestep", lambda d: d["death_prestep"] / d["expo"]), ("other", lambda d: d["death_other"] / d["expo"])]:
            v = []
            for s in seeds:
                x, y = f(IB[b][ci][s]), f(IB[b][cj][s])
                if x > 0 and y > 0: v.append(math.log(x) - math.log(y))
            row[name] = boot_mean(v) if v else (float("nan"),) * 3
        print("   %-9s / %-9s: " % (ci, cj) + " · ".join("%s %s" % (k, fmt(*row[k])) for k in row))
        OUT["A3|%s|%s|%s" % (b, ci, cj)] = {k: {"d": row[k][0], "ci": [row[k][1], row[k][2]]} for k in row}

# ---------------- A4 · A4b ----------------
print("\n== A4 — 이웃 고리 등급의 배경별 등급 평균 ln U_copy 차 (짧은 쪽 − 긴 쪽) · 하한 > 0   (μ %s)" % MUQ)
a4_all = True
for b in NBASES:
    codes = inserts(b, "n")
    Ls = sorted(set(codes.values()))
    for c in codes:
        if Q[b][c]:
            be = sum(d["b_exact"] for d in Q[b][c].values()); nk = sum(d["kids"] for d in Q[b][c].values()); nc = sum(d["cosmic"] for d in Q[b][c].values())
            OUT["code|%s|%s" % (b, c)].update({"U_copy": nk / be, "U_all": (nk + nc) / be, "b_exact": be}) if ("code|%s|%s" % (b, c)) in OUT else None
    seeds = sorted(set.intersection(*[set(Q[b][c]) for c in codes if Q[b][c]])) if any(Q[b][c] for c in codes) else []
    print("  밑 코드 `%s` · 공통 배경 %d · U_copy(합): %s" % (b, len(seeds), " · ".join("%s %.3f" % (c, OUT["code|%s|%s" % (b, c)]["U_copy"]) for c in codes if ("code|%s|%s" % (b, c)) in OUT and "U_copy" in OUT["code|%s|%s" % (b, c)])))
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        v = []
        for s in seeds:
            try:
                m1 = sum(math.log(Q[b][c][s]["U_copy"]) for c in c1) / len(c1)
                m2 = sum(math.log(Q[b][c][s]["U_copy"]) for c in c2) / len(c2)
            except ValueError:
                continue
            v.append(m1 - m2)
        if not v: continue
        m, lo, hi = boot_mean(v)
        ok = lo > 0; a4_all &= ok
        print("     ln U(%d) − ln U(%d) = %s %s" % (L1, L2, fmt(m, lo, hi), "✅" if ok else "❌"))
        OUT["A4|%s|%d|%d" % (b, L1, L2)] = {"d": m, "ci": [lo, hi], "pass": ok, "n_bg": len(v)}
print("  A4 종합: %s" % ("✅" if a4_all else "❌"))
OUT["A4"] = a4_all

pts = [(math.log(v["R"]), math.log(v["U_copy"])) for k, v in OUT.items() if k.startswith("code|") and "U_copy" in v and v["U_copy"] > 0]
if len(pts) >= 3:
    n = len(pts); mx = sum(p[0] for p in pts) / n; my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts); sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    slope = sxy / sxx if sxx else float("nan")
    def rank(a):
        o = sorted(range(len(a)), key=lambda i: a[i]); r = [0] * len(a)
        for i, j in enumerate(o): r[j] = i
        return r
    rx, ry = rank([p[0] for p in pts]), rank([p[1] for p in pts])
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry)); rho = 1 - 6 * d2 / (n * (n * n - 1))
    print("\n== A4b 서술 — %d 코드에서 ln U_copy 대 ln R: 기울기 %.3f · Spearman %.3f" % (n, slope, rho))
    OUT["A4b"] = {"n": n, "slope": slope, "spearman": rho}

# ---------------- 종합 ----------------
print("\n== 종합")
if PILOT:
    print("  파일럿 — 판정 없음. TAU 후보를 위한 같은-고리 쌍 |ln R 비| 분포는 위 A2 요약 줄.")
else:
    print("  A1 %s · A2a %s · A2b %s · A4 %s · 자기 검산 %s · G0 %s" % ("✅" if OUT["A1"] else "❌", "✅" if OUT.get("A2") else "❌", "✅" if OUT.get("A2b") else "❌", "✅" if OUT["A4"] else "❌",
                                                          "✅" if not tr_fail else "❌", "✅" if not (g0a["mismatch"] or g0b["mismatch"]) else "❌"))
for x in notes: print("  🟡", x)
OUT["notes"] = notes; OUT["bad"] = bad; OUT["tr_fail"] = tr_fail
fn = "_RESULT_pilot14.json" if PILOT else "_RESULT_pairs14.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False)
print("기록 →", fn)
