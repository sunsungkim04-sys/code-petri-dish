#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""🟡 탐색 2 (판정 아님) — 왜 사다리(−0.51/틱)와 전수 회귀(−0.10 ± 0.11)가 다른가.

가설: 사다리의 네 코드는 **더한 글자를 거의 고정하고 자리만 바꾼** 것이다(`s` 를 두 자리에, `n` 을 두 자리에).
전수 회귀는 글자 종류 · 자리 · 기능 파괴가 전부 섞인다. 그래서 같은 자료를 **짝지어** 다시 본다:

  같은 밑 코드에 **같은 글자**를 서로 다른 자리에 넣어 고리 길이만 다르게 한 쌍
    → ΔΔ = Δ(긴 고리) − Δ(짧은 고리) 를 배경별로 짝지어 잰다.
  고리 길이가 **같은** 쌍(자리만 다름)은 **음성 대조**다 — ΔΔ ≈ 0 이어야 한다.

새 실행 없음 · 판정 아님. 해부 13 의 설계와 예측을 여기서 뽑는다.
실행: python3 loop13_pairs.py [base]
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
NBASES = os.environ.get("PETRI_NBASES", "racld,rascld,rsacld,racldx").split(",")
N_BOOT = 2000
RNG = random.Random(20260927)


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
    """base 에 글자 ch 를 한 자리씩 넣은 코드들 → {코드: 고리 틱}"""
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]
        L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


# 배경별 Δ 를 코드마다 모은다 (배경 = 시드)
PER = {}
for nb in NBASES:
    for f in sorted(glob.glob(os.path.join(BASE, "%snbr_%s" % (ROOT, nb), "mu0", "dms_mat_racld_s*_nbr_%s.json" % nb))):
        z = json.load(open(f)); s = z["seed"]
        wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        b = sum(s_of(a) for a in wt) / len(wt)
        for a in z["arms"]:
            if a["kind"] == "mut":
                PER.setdefault(a["key"], {})[s] = (s_of(a) - b, a["stopped"] >= 0 and a["series"][-1][2] == 0)

BGS = sorted(set.intersection(*[set(v) for v in PER.values()])) if PER else []
print("== 🟡 탐색 2 — 같은 글자를 다른 자리에 넣은 **짝**으로 본 고리 길이 (판정 아님)")
print("  코드 %d · 공통 배경 %d" % (len(PER), len(BGS)))
IX = [[RNG.randrange(len(BGS)) for _ in BGS] for _ in range(N_BOOT)]


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def paired(a, b):
    """Δ(a) − Δ(b) 의 배경 짝 평균과 95% 구간. 둘 다 살아 있는 배경만."""
    bg = [s for s in BGS if s in PER.get(a, {}) and s in PER.get(b, {})]
    if len(bg) < 10:
        return None
    d = [PER[a][s][0] - PER[b][s][0] for s in bg]
    m = sum(d) / len(d)
    bs = [sum(d[i] for i in ix) / len(ix) for ix in IX]
    dead_a = sum(1 for s in bg if PER[a][s][1]); dead_b = sum(1 for s in bg if PER[b][s][1])
    return m, pct(bs, .025), pct(bs, .975), len(bg), dead_a, dead_b


for base in NBASES:
    print("\n  밑 코드 `%s` (고리 %d) — 같은 글자 · 다른 자리" % (base, loop_ticks(base)))
    for ch in ("n", "x", "s"):
        cand = inserts(base, ch)
        have = {k: L for k, L in cand.items() if k in PER}
        byL = {}
        for k, L in have.items():
            byL.setdefault(L, []).append(k)
        if len(byL) < 2:
            print("    +%s : 고리 길이가 한 종류뿐(%s) — 짝 못 만듦" % (ch, sorted(byL)))
            continue
        Ls = sorted(byL)
        print("    +%s : 고리 %s" % (ch, {L: byL[L] for L in Ls}))
        for i in range(len(Ls) - 1):
            lo, hi = Ls[i], Ls[i + 1]
            for ka in byL[hi]:
                for kb in byL[lo]:
                    r = paired(ka, kb)
                    if r:
                        m, l, h, n, da, db = r
                        print("      ΔΔ %-7s(%d) − %-7s(%d) = %+.3f [%+.3f, %+.3f] · 배경 %d · 불임배경 %d/%d · 틱당 %+.3f"
                              % (ka, hi, kb, lo, m, l, h, n, da, db, m / (hi - lo)))
        # 음성 대조: 같은 고리 길이인데 자리만 다른 쌍
        for L in Ls:
            g = byL[L]
            if len(g) >= 2:
                r = paired(g[0], g[1])
                if r:
                    m, l, h, n, da, db = r
                    print("      [음성] 고리 같음 %s vs %s (%d틱) = %+.3f [%+.3f, %+.3f] · 불임배경 %d/%d"
                          % (g[0], g[1], L, m, l, h, da, db))

print("\n  🚩 이 값들이 해부 11 C 의 −0.510/틱 과 얼마나 맞는지가 해부 13 예측의 출처다.")
print("     음성 대조(고리 같고 자리만 다름)가 0 에서 멀면, '고리 길이' 가 아니라 '자리' 가 재고 있는 것이다.")
