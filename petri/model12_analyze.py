#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 12 판정 — 자유 매개변수 없는 해석 모형이 침입 곡선을 재현하는가
사전등록: 해부12-사전등록-2026-09-20.md §2~§4. 판정선 · 모형 · G 는 결과 전에 고정했다. 배경 단위 짝 부트스트랩 2,000회(P · Q 가 같은 배경).

입력
  P  이웃 DMS(μ 0): 코드 C ∈ {racld · rascld · rsacld · racldx} 의 한 글자 이웃 m 마다 Δ_m(배경별 · 같은 작업의 기준 팔 평균과의 차)
       ω_m = 0                      (그 팔의 표지 계통이 배경 절반 넘게에서 0 이 됨 = 불임)
           = exp(중앙값 Δ_m / G)    (그 밖 · 세대당 상대 성장)
  Q  침입 팔(μ 5칸 · --kids 1): 끼운 코드 그대로인 부모의 출생 수 B_c 와 자식 코드 분포 · 복사 없는 변화(우주선)
       u_c[m] = (자식 m + 우주선 m) / B_c  · U_c = Σ u_c[m]   (WT 는 기준 팔 11개를 합친다)
       실측 Δ_c(μ) = 팔의 s − 기준 팔 평균 s   (해부 6C 와 같은 정의 · 팔 기록이 inv6 과 같은지 대조)
모형 M1 (세대 G = 3000 / 420 = 7.14 · 변종은 제 코드를 그대로 남긴다고 본다)
  W_c = exp(Δ_c(0) / G) · a = W_c (1 − U_c)
  N_c = a^G + Σ_m u_c[m] · W_c · (a^G − ω_m^G) / (a − ω_m)        (표에 없는 자식 = 불임 ω 0)
  예측 Δ_c(μ) = ln N_c(μ) − ln N_WT(μ)                            (μ 0 에선 실측과 같아진다 — 뜻대로)
  M2(서술): ω_m 에 (1 − U_c) 를 곱한다(변종도 부모만큼 오류를 낸다)
판정  허용 tol = max(0.10, 0.25 × |실측 Δ_c(μ) − 실측 Δ_c(0)|)
  P1 씨앗 · P2 rsacld · P3 racldx: μ 0.1% · 0.3% 에서 |예측 − 실측| ≤ tol (둘 다)
  P4 비 [Δ_rsacld(0.3%) − Δ_rsacld(0)] / [Δ_rascld(0.3%) − Δ_rascld(0)]: |예측 − 실측| ≤ 0.10
  서술: 0.5% · 1% (M1 · M2) · U_c · 표에 없는 자식 몫 · 우주선 몫 · 씨앗의 변이 흐름 구성
실행: python3 model12_analyze.py [base]
"""
import glob
import json
import math
import os
import random
import statistics
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv12")
INV6 = os.environ.get("PETRI_ROOT_ORIG", "inv6")
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CHECK_INV6 = os.environ.get("PETRI_CHECK_INV6", "1") == "1"
NCODES = ["racld", "rascld", "rsacld", "racldx"]
CODES = ["rascld", "rsacld", "racldx"]
G = 3000 / 420
JUDGE_MUS = [0.001, 0.003]
N_BOOT = 2000
RNG = random.Random(20260921)
bad = []


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def ci(v):
    return [pct(v, .025), pct(v, .975)]


def fmt(m, c, d=3):
    return ("%+." + str(d) + "f [%+." + str(d) + "f, %+." + str(d) + "f]") % (m, c[0], c[1])


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def common(z, what):
    if not z["fidelity"]["same"]: bad.append("충실도 " + what)
    if not z["checksum_match"]: bad.append("배경 체크섬 " + what)
    if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 " + what)
    if z["grow_extinct"] >= 0 or z["grow_top"][0][0] != "racld": bad.append("배경 " + what)


# ---- P: 이웃 DMS
DM = {}     # (code) -> {m: {bg: (delta, stopped)}}
for code in NCODES:
    DM[code] = {}
    fs = sorted(glob.glob(os.path.join(BASE, "%snbr_%s" % (ROOT, code), "mu0", "dms_mat_racld_s*_nbr_%s.json" % code)))
    keysets = []
    for f in fs:
        z = json.load(open(f)); s = z["seed"]; common(z, "P %s %d" % (code, s))
        if z["mu_assay"] != 0 or z.get("nbr_of") != code or z["ticks"] != 3000: bad.append("P 설정 %s %d" % (code, s))
        wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt) / len(wt)
        ks = set()
        for a in z["arms"]:
            if a["kind"] == "mut":
                DM[code].setdefault(a["key"], {})[s] = (s_of(a) - base, a["stopped"] >= 0 and a["series"][-1][2] == 0)
                ks.add(a["key"])
        keysets.append(ks)
    if len(fs) != NBG or any(k != keysets[0] for k in keysets): bad.append("P 파일 수 · 이웃 목록 %s (%d)" % (code, len(fs)))

# ---- Q: 자식 세기
Q = {}
for mu in MUS:
    for f in sorted(glob.glob(os.path.join(BASE, ROOT + "kids", mtag(mu), "dms_mat_racld_s*_list.json"))):
        z = json.load(open(f)); s = z["seed"]; common(z, "Q %s %d" % (mu, s))
        if abs(z["mu_assay"] - mu) > 1e-12 or not z.get("kids_on") or len(z["arms"]) != 14: bad.append("Q 설정 %s %d" % (mu, s))
        if CHECK_INV6:
            o = json.load(open(os.path.join(BASE, INV6, mtag(mu), os.path.basename(f))))
            if len(o["arms"]) != len(z["arms"]) or any(x["series"] != y["series"] or x["key"] != y["key"] for x, y in zip(o["arms"], z["arms"])):
                bad.append("Q 팔 기록이 inv6 과 다르다 %s %d" % (mu, s))
        Q[(mu, s)] = z
BGS = sorted({k[1] for k in Q})
BGS = [s for s in BGS if all((mu, s) in Q for mu in MUS) and all(s in next(iter(DM[c].values())) for c in NCODES if DM[c])]
print("== 해부 12 판정 — 해석 모형 M1 (G = %.2f 세대 · 자유 매개변수 없음)" % G)
print("  P 이웃 수: " + " · ".join("%s %d" % (c, len(DM[c])) for c in NCODES) + " · Q 파일 %d/%d · 배경 %d · 점검 문제 %d%s" % (len(Q), len(MUS) * NBG, len(BGS), len(bad), "" if CHECK_INV6 else " · (inv6 대조 끔)"))
tb = sum(a["kids"]["b_exact"] + a["kids"]["b_other"] + a["kids"]["unknown"] for z in Q.values() for a in z["arms"])
unk = sum(a["kids"]["unknown"] for z in Q.values() for a in z["arms"])
print("  표지 계통 출생 %d · 부모를 못 찾은 몫 %.4f%% (선 1%%)" % (tb, 100 * unk / max(1, tb)))
if bad or len(BGS) != NBG or unk / max(1, tb) >= 0.01:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:12]:
        print("    ", b)
    sys.exit(1)


def omega_table(ix):
    """부트스트랩 표본(배경 색인)에서 자식 코드 → ω. 부모 코드의 이웃 표를 먼저, 없으면 다른 표를 쓴다."""
    T = {}
    for code in NCODES:
        t = {}
        for m, per in DM[code].items():
            vals = [per[BGS[i]] for i in ix]
            if sum(1 for v in vals if v[1]) * 2 > len(vals): t[m] = 0.0
            else: t[m] = math.exp(statistics.median(v[0] for v in vals) / G)
        T[code] = t
    return T


def flows(code, mu, ix):
    """끼운 코드 그대로인 부모의 출생당 자식 코드 흐름 u[m]. WT(racld)는 기준 팔 11개를 합친다."""
    B = 0; K = {}; cos = 0
    for i in ix:
        z = Q[(mu, BGS[i])]
        arms = [a for a in z["arms"] if (a["kind"] in ("ref", "neutral") if code == "racld" else (a["kind"] == "mut" and a["key"] == code))]
        for a in arms:
            B += a["kids"]["b_exact"]
            for k, n in a["kids"]["kids"]: K[k] = K.get(k, 0) + n
            for k, n in a["kids"]["cosmic"]: K[k] = K.get(k, 0) + n; cos += n
    return B, K, cos


def meas(code, mu, ix):
    v = []
    for i in ix:
        z = Q[(mu, BGS[i])]
        wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt) / len(wt)
        v.append(next(s_of(a) for a in z["arms"] if a["kind"] == "mut" and a["key"] == code) - base)
    return sum(v) / len(v)


def lnN(code, mu, ix, T, W, m2=False, detail=None):
    B, K, cos = flows(code, mu, ix)
    if B == 0:
        return float("nan")
    U = sum(K.values()) / B
    a = W * (1 - U); tot = a ** G
    for m, n in K.items():
        if m == "racld": w = 1.0
        elif m in T[code]: w = T[code][m]
        else:
            w = next((T[c][m] for c in NCODES if m in T[c]), None)
            if w is None:
                w = 0.0
                if detail is not None: detail["untabled"] = detail.get("untabled", 0) + n
        if m2: w *= (1 - U)
        u = n / B
        tot += u * W * (G * a ** (G - 1) if abs(a - w) < 1e-9 else (a ** G - w ** G) / (a - w))
        if detail is not None:
            cls = "불임" if w == 0 else ("WT 이상" if w >= 1 else ("부모 이상" if w >= W else "느림"))
            detail.setdefault("cls", {}); detail["cls"][cls] = detail["cls"].get(cls, 0) + n
    if detail is not None: detail.update({"U": U, "B": B, "cosmic": cos, "mut": sum(K.values())})
    return math.log(tot)


def predict(code, mu, ix, m2=False):
    T = omega_table(ix)
    W = math.exp(meas(code, 0.0, ix) / G)
    return lnN(code, mu, ix, T, W, m2) - lnN("racld", mu, ix, T, 1.0, m2)


ALL = list(range(len(BGS)))
IX = [[RNG.randrange(len(BGS)) for _ in BGS] for _ in range(N_BOOT)]
OUT = {"G": G, "bg": BGS}
print("\n-- 변이 흐름 (μ 0.3% · 끼운 코드 그대로인 부모의 출생당)")
T0 = omega_table(ALL)
for code in ["racld"] + CODES:
    d = {}
    lnN(code, 0.003, ALL, T0, 1.0 if code == "racld" else math.exp(meas(code, 0.0, ALL) / G), detail=d)
    cl = d.get("cls", {}); mt = max(1, d["mut"])
    print("  %-6s U %.4f (출생 %d) · 우주선 몫 %.3f · 표에 없는 자식 몫 %.3f · 구성: %s" % (code, d["U"], d["B"], d["cosmic"] / mt, d.get("untabled", 0) / mt, " · ".join("%s %.2f" % (k, v / mt) for k, v in sorted(cl.items(), key=lambda kv: -kv[1]))))
    OUT["flow|" + code] = {"U": d["U"], "cosmic_share": d["cosmic"] / mt, "untabled_share": d.get("untabled", 0) / mt, "cls": {k: v / mt for k, v in cl.items()}}
print("  이웃 표의 불임 몫: " + " · ".join("%s %d/%d" % (c, sum(1 for v in T0[c].values() if v == 0), len(T0[c])) for c in NCODES))

print("\n-- 예측 대 실측 Δ (M1 · 괄호는 M2)")
print("      코드 ·      μ |   실측   |   예측 M1  (M2)   | 예측 − 실측 [구간] | 허용 | ")
verd = {}
for code in CODES:
    m0 = meas(code, 0.0, ALL)
    ok = True
    for mu in MUS[1:]:
        me = meas(code, mu, ALL); p1 = predict(code, mu, ALL); p2 = predict(code, mu, ALL, True)
        bs = [predict(code, mu, ix) - meas(code, mu, ix) for ix in IX]
        tol = max(0.10, 0.25 * abs(me - m0)); judged = mu in JUDGE_MUS
        hit = abs(p1 - me) <= tol
        if judged: ok &= hit
        print("  %8s · %6g | %+.3f | %+.3f (%+.3f) | %s | %.3f %s" % (code, mu, me, p1, p2, fmt(p1 - me, ci(bs)), tol, ("✅" if hit else "❌") if judged else "(서술)"))
        OUT["cell|%s|%g" % (code, mu)] = {"meas": me, "M1": p1, "M2": p2, "diff": [p1 - me, ci(bs)], "tol": tol, "judged": judged, "hit": hit}
    verd[code] = ok
    print("  %8s · μ 0 실측 Δ %+.3f (모형은 여기서 실측과 같다 — 뜻대로)" % (code, m0))


def ratio(fn, ix):
    return (fn("rsacld", 0.003, ix) - meas("rsacld", 0.0, ix)) / (fn("rascld", 0.003, ix) - meas("rascld", 0.0, ix))


rp, rm = ratio(predict, ALL), ratio(meas, ALL)
rb = [ratio(predict, ix) - ratio(meas, ix) for ix in IX]
p4 = abs(rp - rm) <= 0.10
print("\n  P1 씨앗 %s · P2 rsacld %s · P3 racldx %s" % tuple("✅" if verd[c] else "❌" for c in CODES))
print("  P4 rsacld/씨앗 하락 비(μ 0.3%%): 예측 %.3f · 실측 %.3f · 차 %s → %s" % (rp, rm, fmt(rp - rm, ci(rb)), "✅" if p4 else "❌"))
allok = all(verd.values()) and p4
print("  → " + ("자유 매개변수 없는 모형이 낮은 μ 의 침입 곡선 셋과 코드 사이 비를 재현한다" if allok else "모형 M1 은 사전등록한 허용 안에서 전부를 재현하지는 못한다 — 어긋난 칸을 그대로 적는다"))
OUT["verdict"] = {"P1": verd["rascld"], "P2": verd["rsacld"], "P3": verd["racldx"], "P4": p4, "ratio": {"pred": rp, "meas": rm, "diff_ci": ci(rb)}}
tag = os.environ.get("PETRI_OUT", "_RESULT_model12")
json.dump(OUT, open(os.path.join(BASE, tag + ".json"), "w"), ensure_ascii=False, indent=1)
print("\n기록 → %s.json" % tag)
