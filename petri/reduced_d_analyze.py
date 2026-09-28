#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""축소 모형 D 요약 — 설계 축소모형-D-설계-2026-09-27.md §3. 판정 아님(D2 는 비를 적는다).
  D1  등급 사이 S = (R−1)·L 의 max/min (≤ 1.15 면 항등식 확인)
  D2  S_잘섞임 / S_접시 — 코드마다(접시 값은 _RESULT_pairs14.json 의 code|*)
  D3  기억 규칙: R → 1 · μ 0 ΔΔ 그대로(원래 대 기억 차) · Δμ = ΔΔ(0.003) − ΔΔ(0) 원래 대 기억
  D4  ln U 대 ln R 기울기
실행: python3 reduced_d_analyze.py <reduced_d_run.json> [_RESULT_pairs14.json]
"""
import json
import math
import sys

z = json.load(open(sys.argv[1]))
dish = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else {}
codes = z["codes"]; arms = z["arm"]
seeds = sorted(set(a["seed"] for a in arms))
print("== 축소 모형 D — 칸 %d · 성장 %d · 시드 %s · 코드 %d · 규칙 %s · μ %s" % (z["args"]["cells"], z["args"]["grow"], seeds, len(codes), z["args"]["rules"], z["args"]["mus"]))
for c in z["cells"]: print("   시드 %d 정상 개체 %d · 자유 글자 %d · 소비 글자 중 자유 %.1f%%" % (c["seed"], c["pop"], c["free"], 100 * c["free_used_frac"]))


def get(rule, mu, code):
    return [a for a in arms if a["rule"] == rule and abs(a["mu"] - mu) < 1e-12 and a["code"] == code]


def mean(v): return sum(v) / len(v) if v else float("nan")


# ---- D1 · D2 ----
print("\n== D1 · D2 — 원래 규칙 μ 0: R · T · S 와 접시 대비")
fam = {}
for code, L in codes.items():
    base = code.replace("n", "", 1)
    fam.setdefault(base, {}).setdefault(L, []).append(code)
ratios = []
for base, cls in fam.items():
    print("  밑 코드 `%s`" % base)
    Sc = {}
    for L in sorted(cls):
        for code in cls[L]:
            g = get("orig", 0.0, code)
            R = mean([a["R"] for a in g]); T = mean([a["T"] for a in g]); S = (R - 1) * L
            Sc.setdefault(L, []).append(S)
            dk = "code|%s|%s" % (base, code)
            dS = dish.get(dk, {}).get("S"); dR = dish.get(dk, {}).get("R")
            rat = S / dS if dS else float("nan")
            if dS: ratios.append(rat)
            print("     %-8s L %d · R %.2f · T %.1f · T/R−L %+.2f · S %.1f · 접시 S %s · 비 %s" % (code, L, R, T, T / R - L, S, ("%.1f" % dS) if dS else "—", ("%.2f" % rat) if dS else "—"))
    m = {L: mean(v) for L, v in Sc.items()}
    print("     D1: 등급 평균 S %s · max/min %.3f" % (" · ".join("L%d %.1f" % (L, v) for L, v in sorted(m.items())), max(m.values()) / min(m.values())))
    Ls = sorted(m)
    for L1, L2 in zip(Ls, Ls[1:]):
        r1 = mean([mean([a["R"] for a in get("orig", 0.0, c)]) for c in cls[L1]]); r2 = mean([mean([a["R"] for a in get("orig", 0.0, c)]) for c in cls[L2]])
        print("     A1 꼴: ln R(%d) − ln R(%d) = %+.3f (접시 %s)" % (L1, L2, math.log(r1) - math.log(r2), ("%+.3f" % dish["A1|%s|%d|%d" % (base, L1, L2)]["d"]) if ("A1|%s|%d|%d" % (base, L1, L2)) in dish else "—"))
if ratios: print("  D2: S_잘섞임 / S_접시 — 중앙값 %.2f · 범위 %.2f ~ %.2f (코드 %d)" % (sorted(ratios)[len(ratios) // 2], min(ratios), max(ratios), len(ratios)))

# ---- D3 ----
print("\n== D3 — 기억 규칙: R · μ 0 ΔΔ · Δμ")
for base, cls in fam.items():
    Ls = sorted(cls)
    if len(Ls) < 2: continue
    Rrem = mean([mean([a["R"] for a in get("rem", 0.0, c)]) for L in Ls for c in cls[L]])
    Rrem3 = mean([mean([a["R"] for a in get("rem", 0.003, c)]) for L in Ls for c in cls[L]])
    print("  밑 코드 `%s`: 기억 규칙 R 평균 μ 0 %.3f · μ 0.3%% %.3f" % (base, Rrem, Rrem3))
    for L1, L2 in zip(Ls, Ls[1:]):
        rows = []
        for rule in ("orig", "rem"):
            dd = {}
            for mu in (0.0, 0.003):
                v = []
                for s_ in seeds:
                    s1 = mean([a["s"] for c in cls[L1] for a in get(rule, mu, c) if a["seed"] == s_])
                    s2 = mean([a["s"] for c in cls[L2] for a in get(rule, mu, c) if a["seed"] == s_])
                    v.append(s2 - s1)
                dd[mu] = mean(v)
            rows.append("%s ΔΔ(0) %+.3f · ΔΔ(0.3%%) %+.3f · Δμ %+.3f" % (rule, dd[0.0], dd[0.003], dd[0.003] - dd[0.0]))
        print("     L %d→%d: " % (L1, L2) + " | ".join(rows))

# ---- D4 ----
pts = []
for code in codes:
    g0 = get("orig", 0.0, code); g3 = get("orig", 0.003, code)
    R = mean([a["R"] for a in g0]); U = mean([a["U"] for a in g3 if not math.isnan(a["U"])])
    if R > 0 and U > 0: pts.append((math.log(R), math.log(U)))
if len(pts) >= 3:
    n = len(pts); mx = mean([p[0] for p in pts]); my = mean([p[1] for p in pts])
    sxx = sum((p[0] - mx) ** 2 for p in pts); sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    print("\n== D4 — %d 코드 ln U(0.3%%) 대 ln R: 기울기 %.3f (접시 %s)" % (n, sxy / sxx, ("%.3f" % dish["A4b"]["slope"]) if "A4b" in dish else "—"))
