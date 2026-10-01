#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 22 판정기 — 축소 모형 D 씨앗 20. 사전등록 해부22-사전등록-2026-10-01.md §4.
집계 정의는 원래 판 요약기(reduced_d_analyze.py)와 같다: 코드마다 시드별 R 의 평균 → S = (R̄ − 1)·L.
불확실성: 시드(배경) 부트스트랩 2,000회(난수 22) · 백분위 95% 구간. 접시 값(_RESULT_pairs14.json)은 상수로 둔다.

  G0  구현 sha256 = d9473aa0…(원래 판) · 시드 목록 · 파일마다 설정 = 원래 판 설정 · arm 104 · 코드 13 · 끼우기 m0 > 0 · 배경 개체 > 0
  G1  (--orig 를 주면) 공통 시드의 arm 이 원래 판 _RESULT_reduced_d.json 과 같은가 — 정수 칸 일치 · 실수 칸 상대 1e-12 · 정규화 sha256
  P1  D2 중앙 비 M = median_c(S_c / S접시_c) · 구간 ⊂ (0.50, 1.00)                  [주 판정]
  P2  가족 대비 F = 평균 비(racld 가족 6) − 평균 비(rascld 가족 7) · 하한 > 0         [주 판정]
  P3  A1 꼴 ln R̄(4) − ln R̄(5) (racld) · ln R̄(2) − ln R̄(3) (rascld) · 둘 다 하한 > 0  [부 판정]
  P4  D4 기울기 ln U(0.3%) 대 ln R(0) · 구간 ⊂ [0.484, 0.884] (접시 0.684 ± 0.20)      [부 판정]
  서술: 코드별 비와 구간 · 범위 · D1 max/min · 원래 판 2 시드 값과의 차 · 자유 글자 몫 · D3(ΔΔ · Δμ · 기억 규칙 R) — 판정 아님
실행: python3 red22_analyze.py --root red22 --seeds 2201-2220 --dish _RESULT_pairs14.json [--orig _RESULT_reduced_d.json] [--out _RESULT_reduced22.json]
"""
import argparse
import glob
import hashlib
import json
import math
import os
import random
import sys

IMPL_SHA = "d9473aa0c1faaecaeb36df11d6cb3e6a62194f86e65c34ab0ec5a258da70ea43"
EXPECT = {"cells": 1000, "density": 8, "grow": 6000, "ticks": 3000, "seeds": 1, "codes_of": "racld,rascld",
          "mus": "0,0.003", "rules": "orig,rem", "resident": "racld", "frac": 0.1}
FAM = {"racld": ["nracld", "rnacld", "rancld", "racnld", "raclnd", "racldn"],
       "rascld": ["nrascld", "rnascld", "ranscld", "rasncld", "rascnld", "rasclnd", "rascldn"]}
CLASS = {"racld": (4, 5), "rascld": (2, 3)}
ORIG2 = {"M": 0.63, "range": (0.55, 0.91), "A1|racld": 0.261, "A1|rascld": 0.232, "D4": 0.639}  # 원래 판 _RESULT_reduced_d.txt
P1_BAND = (0.50, 1.00); P4_BAND = (0.684 - 0.20, 0.684 + 0.20)
NBOOT = 2000; BOOT_SEED = 22

ap = argparse.ArgumentParser()
ap.add_argument("--root", default="red22"); ap.add_argument("--seeds", default="2201-2220")
ap.add_argument("--dish", default="_RESULT_pairs14.json"); ap.add_argument("--orig", default=None)
ap.add_argument("--out", default=None)
a = ap.parse_args()
if "-" in a.seeds:
    lo, hi = map(int, a.seeds.split("-")); SEEDS = list(range(lo, hi + 1))
else:
    SEEDS = [int(x) for x in a.seeds.split(",")]
dish = json.load(open(a.dish))
problems = []
OUT = {"root": a.root, "seeds": SEEDS}

# ---------------- G0 ----------------
imp = open(os.path.join(a.root, "impl.sha256")).read().split()[0] if os.path.exists(os.path.join(a.root, "impl.sha256")) else None
if imp != IMPL_SHA: problems.append("impl.sha256 %s ≠ 원래 판 %s" % (imp, IMPL_SHA))
here = hashlib.sha256(open("reduced_d.py", "rb").read()).hexdigest() if os.path.exists("reduced_d.py") else None
if here != IMPL_SHA: problems.append("현재 reduced_d.py sha256 %s ≠ 원래 판" % here)
files = sorted(glob.glob(os.path.join(a.root, "s*.json")))
runs = {}
for f in files:
    z = json.load(open(f)); runs[z["args"]["seed0"]] = z
if sorted(runs) != sorted(SEEDS): problems.append("시드 목록 불일치: 파일 %s · 기대 %s" % (sorted(runs), SEEDS))
codes = None
for s, z in sorted(runs.items()):
    for k, v in EXPECT.items():
        if z["args"].get(k) != v: problems.append("시드 %d 설정 %s = %r ≠ %r" % (s, k, z["args"].get(k), v))
    if codes is None: codes = z["codes"]
    if z["codes"] != codes: problems.append("시드 %d 코드 목록이 다르다" % s)
    if len(z["arm"]) != 2 * 2 * 13: problems.append("시드 %d arm %d ≠ 52" % (s, len(z["arm"])))
    if any(r["m0"] <= 0 for r in z["arm"]): problems.append("시드 %d 끼우기 0 인 arm 있음" % s)
    if z["cells"][0]["pop"] <= 0: problems.append("시드 %d 배경 멸종" % s)
    # 검토 개정 10-01: 판정에 쓰는 R 이 NaN(복제 0) 이거나 ≤ 0 이면 mean() 이 그 시드를 조용히 뺀다 → G0 에서 멈춘다
    nbadR = sum(1 for r in z["arm"] if not (isinstance(r["R"], (int, float)) and r["R"] == r["R"] and r["R"] > 0))
    if nbadR: problems.append("시드 %d R 이 NaN/≤0 인 arm %d 줄" % (s, nbadR))
sig = {}
for s, z in runs.items():
    sig.setdefault(json.dumps([(r["R"], r["mT"]) for r in z["arm"]]), []).append(s)
for v in sig.values():
    if len(v) > 1: problems.append("시드 %s 의 arm 이 서로 똑같다 — 시드가 안 먹었나" % v)
if codes is not None and sorted(codes) != sorted(FAM["racld"] + FAM["rascld"]): problems.append("코드 13 목록 불일치")
print("== 해부 22 — 축소 모형 D · 루트 %s · 시드 %d 개(%d~%d) · 구현 %s…" % (a.root, len(SEEDS), min(SEEDS), max(SEEDS), (imp or "없음")[:8]))
print("G0 점검: %s" % ("문제 0" if not problems else "문제 %d" % len(problems)))
for p in problems: print("   ✗", p)
if problems: sys.exit(2)
for s in SEEDS:
    c = runs[s]["cells"][0]
    print("   시드 %d 정상 개체 %d · 소비 글자 중 자유 %.1f%%" % (s, c["pop"], 100 * c["free_used_frac"]))
OUT["cells"] = [runs[s]["cells"][0] for s in SEEDS]


# ---------------- G1 (원래 판과의 동일성) ----------------
def canon(row):
    return json.dumps({k: (("%.12g" % v) if isinstance(v, float) else v) for k, v in sorted(row.items())}, sort_keys=True, ensure_ascii=False)


if a.orig and not (set(r["seed"] for r in json.load(open(a.orig))["arm"]) & set(SEEDS)):
    print("G1 원래 판 대조: 공통 시드 없음 — 동일성은 마른 실행(시드 101 102)에서 확인한다")
elif a.orig:
    zo = json.load(open(a.orig)); common = sorted(set(r["seed"] for r in zo["arm"]) & set(SEEDS))
    ro = [r for r in zo["arm"] if r["seed"] in common]; rn = [r for s in common for r in runs[s]["arm"]]
    bad = 0
    for x, y in zip(ro, rn):
        for k in x:
            if isinstance(x[k], float) and not (math.isnan(x[k]) and math.isnan(y[k])):
                if abs(x[k] - y[k]) > 1e-12 * max(1.0, abs(x[k])): bad += 1
            elif not isinstance(x[k], float) and x[k] != y[k]: bad += 1
    ho = hashlib.sha256("\n".join(canon(r) for r in ro).encode()).hexdigest()
    hn = hashlib.sha256("\n".join(canon(r) for r in rn).encode()).hexdigest()
    cells_same = all(zo["cells"][i] == runs[s]["cells"][0] for i, s in enumerate(common)) if common == sorted(r["seed"] for r in zo["cells"]) else None
    print("G1 원래 판 대조: 공통 시드 %s · arm %d 대 %d · 다른 칸 %d · 배경 요약 같음 %s · 정규화 sha256 원래 %s… · 이번 %s… → %s" % (
        common, len(ro), len(rn), bad, cells_same, ho[:12], hn[:12], "같음 ✅" if (ho == hn and bad == 0 and len(ro) == len(rn) and len(ro) > 0) else "다름 ✗"))
    OUT["G1"] = {"common": common, "n_orig": len(ro), "n_new": len(rn), "bad": bad, "sha_orig": ho, "sha_new": hn, "cells_same": cells_same}

# ---------------- 통계 ----------------
ARM = {}
for s in SEEDS:
    for r in runs[s]["arm"]:
        ARM[(s, r["rule"], round(r["mu"], 6), r["code"])] = r
L = {c: codes[c] for c in codes}
Sdish = {c: dish["code|%s|%s" % (b, c)]["S"] for b in FAM for c in FAM[b]}


def mean(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return sum(v) / len(v) if v else float("nan")


def median(v):
    v = sorted(v); n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def Rbar(ss, c, rule="orig", mu=0.0): return mean([ARM[(s, rule, mu, c)]["R"] for s in ss])
def Ubar(ss, c): return mean([ARM[(s, "orig", 0.003, c)]["U"] for s in ss])


def stats(ss):
    o = {}
    rat = {c: (Rbar(ss, c) - 1) * L[c] / Sdish[c] for c in L}
    o["ratio"] = rat; o["M"] = median(list(rat.values()))
    o["F"] = mean([rat[c] for c in FAM["racld"]]) - mean([rat[c] for c in FAM["rascld"]])
    for b, (l1, l2) in CLASS.items():
        r1 = mean([Rbar(ss, c) for c in FAM[b] if L[c] == l1]); r2 = mean([Rbar(ss, c) for c in FAM[b] if L[c] == l2])
        o["A1|" + b] = math.log(r1) - math.log(r2)
        Sc = {l: mean([(Rbar(ss, c) - 1) * L[c] for c in FAM[b] if L[c] == l]) for l in (l1, l2)}
        o["D1|" + b] = max(Sc.values()) / min(Sc.values())
    pts = [(math.log(Rbar(ss, c)), math.log(Ubar(ss, c))) for c in L if Ubar(ss, c) > 0]
    mx = mean([p[0] for p in pts]); my = mean([p[1] for p in pts])
    o["D4"] = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts); o["D4_n"] = len(pts)
    for b, (l1, l2) in CLASS.items():
        for rule in ("orig", "rem"):
            dd = {}
            for mu in (0.0, 0.003):
                v = []
                for s in ss:
                    s1 = mean([ARM[(s, rule, mu, c)]["s"] for c in FAM[b] if L[c] == l1])
                    s2 = mean([ARM[(s, rule, mu, c)]["s"] for c in FAM[b] if L[c] == l2])
                    v.append(s2 - s1)
                dd[mu] = mean(v)
            o["DD0|%s|%s" % (b, rule)] = dd[0.0]; o["DD3|%s|%s" % (b, rule)] = dd[0.003]; o["Dmu|%s|%s" % (b, rule)] = dd[0.003] - dd[0.0]
        o["Rrem|" + b] = mean([Rbar(ss, c, "rem", 0.0) for c in FAM[b]])
    return o


point = stats(SEEDS)
rng = random.Random(BOOT_SEED)
boots = [stats([SEEDS[int(rng.random() * len(SEEDS))] for _ in SEEDS]) for _ in range(NBOOT)]


NAN_DROP = {}


def ci(key, sub=None):
    v = sorted((b[key][sub] if sub else b[key]) for b in boots)
    v = [x for x in v if not math.isnan(x)]
    if len(v) < len(boots): NAN_DROP[key + ("|" + sub if sub else "")] = len(boots) - len(v)  # 검토 개정: 조용히 버린 표본을 센다
    return (v[int(0.025 * len(v))], v[min(len(v) - 1, int(0.975 * len(v)))])


def f3(x): return "%+.3f" % x


res = {}
M = point["M"]; cM = ci("M")
res["P1"] = {"M": M, "ci": cM, "pass": cM[0] > P1_BAND[0] and cM[1] < P1_BAND[1]}
F = point["F"]; cF = ci("F")
res["P2"] = {"F": F, "ci": cF, "pass": cF[0] > 0}
a1 = {b: (point["A1|" + b], ci("A1|" + b)) for b in CLASS}
res["P3"] = {b: {"d": a1[b][0], "ci": a1[b][1]} for b in CLASS}; res["P3"]["pass"] = all(a1[b][1][0] > 0 for b in CLASS)
cD4 = ci("D4")
res["P4"] = {"slope": point["D4"], "ci": cD4, "n": point["D4_n"], "pass": cD4[0] >= P4_BAND[0] and cD4[1] <= P4_BAND[1]}

print("\n== 주 판정")
print("P1 D2 중앙 비 M = %.3f [%.3f, %.3f] · 판정선 구간 ⊂ (%.2f, %.2f) → %s   (원래 판 2 시드 %.2f · 차 %s)" % (
    M, cM[0], cM[1], P1_BAND[0], P1_BAND[1], "✅" if res["P1"]["pass"] else "❌", ORIG2["M"], f3(M - ORIG2["M"])))
print("P2 가족 대비 F = 평균 비 racld − rascld = %s [%s, %s] · 하한 > 0 → %s" % (f3(F), f3(cF[0]), f3(cF[1]), "✅" if res["P2"]["pass"] else "❌"))
print("\n== 부 판정")
for b, (l1, l2) in CLASS.items():
    print("P3 %s: ln R̄(%d) − ln R̄(%d) = %s [%s, %s] (원래 판 %s · 접시 %s)" % (b, l1, l2, f3(a1[b][0]), f3(a1[b][1][0]), f3(a1[b][1][1]), f3(ORIG2["A1|" + b]),
          f3(dish["A1|%s|%d|%d" % (b, l1, l2)]["d"]) if ("A1|%s|%d|%d" % (b, l1, l2)) in dish else "—"))
print("P3 둘 다 하한 > 0 → %s" % ("✅" if res["P3"]["pass"] else "❌"))
print("P4 D4 기울기 %.3f [%.3f, %.3f] (코드 %d) · 판정선 구간 ⊂ [%.3f, %.3f] → %s   (원래 판 %.3f · 접시 %.3f)" % (
    point["D4"], cD4[0], cD4[1], point["D4_n"], P4_BAND[0], P4_BAND[1], "✅" if res["P4"]["pass"] else "❌", ORIG2["D4"], dish["A4b"]["slope"]))

d4n = [b["D4_n"] for b in boots]
nanU = sum(1 for s in SEEDS for r in runs[s]["arm"] if r["rule"] == "orig" and abs(r["mu"] - 0.003) < 1e-9 and r["U"] != r["U"])
print("   D4 부트스트랩 코드 수 %d ~ %d (13 이 아니면 U = 0 인 코드가 빠진 표본) · μ 0.3%% 원래 규칙 arm 중 U NaN %d 줄" % (min(d4n), max(d4n), nanU))
res["P4"]["boot_n"] = [min(d4n), max(d4n)]; res["P4"]["nanU"] = nanU
print("\n== 서술(판정 아님)")
print("코드별 S_잘섞임 / S_접시:")
for b in FAM:
    for c in FAM[b]:
        cc = ci("ratio", c)
        print("   %-8s L %d · 비 %.3f [%.3f, %.3f] · S %.2f · 접시 S %.2f" % (c, L[c], point["ratio"][c], cc[0], cc[1], (Rbar(SEEDS, c) - 1) * L[c], Sdish[c]))
rv = list(point["ratio"].values())
print("   범위 %.2f ~ %.2f (원래 판 %.2f ~ %.2f) · 가족 평균 racld %.3f · rascld %.3f" % (min(rv), max(rv), ORIG2["range"][0], ORIG2["range"][1],
      mean([point["ratio"][c] for c in FAM["racld"]]), mean([point["ratio"][c] for c in FAM["rascld"]])))
for b in CLASS:
    cc = ci("D1|" + b)
    print("D1 %s 등급 평균 S max/min %.3f [%.3f, %.3f] (≤ 1.15 이면 항등식 확인 · 서술)" % (b, point["D1|" + b], cc[0], cc[1]))
print("D3 (🚫 고리 길이가 Δ 를 예측한다는 읽기 아님 — 등급 평균의 로그비 차 서술):")
for b, (l1, l2) in CLASS.items():
    cr = ci("Rrem|" + b)
    print("   %s 기억 규칙 R(μ 0) %.3f [%.3f, %.3f]" % (b, point["Rrem|" + b], cr[0], cr[1]))
    for rule in ("orig", "rem"):
        k0, k3, km = "DD0|%s|%s" % (b, rule), "DD3|%s|%s" % (b, rule), "Dmu|%s|%s" % (b, rule)
        c0, c3, cm = ci(k0), ci(k3), ci(km)
        print("   %s L %d→%d %s: ΔΔ(0) %s [%s, %s] · ΔΔ(0.3%%) %s [%s, %s] · Δμ %s [%s, %s]" % (b, l1, l2, rule, f3(point[k0]), f3(c0[0]), f3(c0[1]),
              f3(point[k3]), f3(c3[0]), f3(c3[1]), f3(point[km]), f3(cm[0]), f3(cm[1])))
fu = [runs[s]["cells"][0]["free_used_frac"] for s in SEEDS]; pp = [runs[s]["cells"][0]["pop"] for s in SEEDS]
print("배경: 개체 %d ~ %d (평균 %.0f) · 소비 글자 중 자유 %.1f ~ %.1f%% (평균 %.1f%%)" % (min(pp), max(pp), mean(pp), 100 * min(fu), 100 * max(fu), 100 * mean(fu)))

res["describe"] = {"ratio": {c: {"r": point["ratio"][c], "ci": ci("ratio", c)} for c in L}, "range": [min(rv), max(rv)],
                   "D1": {b: {"v": point["D1|" + b], "ci": ci("D1|" + b)} for b in CLASS},
                   "D3": {k: {"v": point[k], "ci": ci(k)} for k in point if k.startswith(("DD", "Dmu", "Rrem"))}}
print("부트스트랩 NaN 으로 버린 표본: %s" % (NAN_DROP if NAN_DROP else "0"))
OUT["result"] = res; OUT["nboot"] = NBOOT; OUT["boot_seed"] = BOOT_SEED; OUT["nan_drop"] = NAN_DROP
if a.out:
    json.dump(OUT, open(a.out, "w"), indent=1, ensure_ascii=False)
    print("기록 →", a.out)
