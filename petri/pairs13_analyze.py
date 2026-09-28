#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 13 판정 — 둘째 코드 가족에서 짝지은 고리 사다리가 서는가
사전등록: 해부13-사전등록-2026-09-27.md §3~§4. 판정선은 결과 전에 고정했다(예측 구간의 출처는 §2 사후 재분석).

  밑 코드 B 에 **같은 글자를 다른 자리에** 넣은 쌍만 본다. 더한 글자는 `n`(쉬기 · 기능 없음)으로 고정.
    ΔΔ = Δ(고리 L+1 인 코드) − Δ(고리 L 인 코드)  · 배경 짝 · 부트스트랩 2,000
  N1 모든 `n` 짝의 ΔΔ 상한 < 0
  N2 틱당 점추정이 [−1.0, −0.2] 안 — 모든 짝
  N3 음성 대조(고리 같고 자리만 다른 `n` 쌍) |ΔΔ| ≤ 0.15 — 모든 쌍
  N4 서술 — `x` 삽입(전부 같은 고리)의 모든 쌍
  N5 서술 — μ 0.3% 의 같은 짝들
실행: python3 pairs13_analyze.py [base]
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv13")
NBASES = os.environ.get("PETRI_NBASES", "acld,rascled").split(",")
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.003]
NBG = int(os.environ.get("PETRI_NBG", "20"))
JUDGE_MU = 0.0
TOL_NEG = 0.15
PER_TICK_LO, PER_TICK_HI = -1.0, -0.2
N_BOOT = 2000
RNG = random.Random(20260927)
bad = []
OUT = {}


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


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


PER = {}      # (base, mu) -> {code: {seed: (delta, dead)}}
SELF = {}     # (base, mu) -> {seed: (delta, dead)}  밑 코드 자신
for b in NBASES:
    for mu in MUS:
        d, sf = {}, {}
        fs = sorted(glob.glob(os.path.join(BASE, "%snbr_%s" % (ROOT, b), mtag(mu), "dms_mat_racld_s*_nbr_%s.json" % b)))
        keysets = []
        for f in fs:
            z = json.load(open(f)); s = z["seed"]
            if abs(z["mu_assay"] - mu) > 1e-12 or z.get("nbr_of") != b: bad.append("설정 %s %s %d" % (b, mu, s))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (b, mu, s))
            if not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (b, mu, s))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (b, mu, s))
            if z["grow_extinct"] >= 0 or z["grow_top"][0][0] != "racld": bad.append("배경 %s %s %d" % (b, mu, s))
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base_s = sum(s_of(a) for a in wt) / len(wt)
            ks = set()
            for a in z["arms"]:
                if a["kind"] != "mut":
                    continue
                v = (s_of(a) - base_s, a["stopped"] >= 0 and a["series"][-1][2] == 0)
                ks.add(a["key"])
                if a["key"] == b: sf[s] = v
                else: d.setdefault(a["key"], {})[s] = v
            keysets.append(ks)
        if len(fs) != NBG or (keysets and any(k != keysets[0] for k in keysets)):
            bad.append("파일 수 · 이웃 목록 %s %s (%d)" % (b, mu, len(fs)))
        PER[(b, mu)] = d; SELF[(b, mu)] = sf

print("== 해부 13 판정 — 둘째 코드 가족에서 짝지은 고리 사다리")
for b in NBASES:
    print("  밑 코드 `%s`(고리 %s) · 이웃 %d · μ %s" % (b, loop_ticks(b), len(PER[(b, MUS[0])]), MUS))
print("  점검 문제 %d" % len(bad))
if bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for x in bad[:10]: print("    ", x)
    sys.exit(1)

verdict = {}
for b in NBASES:
    sf = SELF[(b, JUDGE_MU)]
    dead_self = sum(1 for v in sf.values() if v[1])
    bgs = sorted(sf)
    print("\n-- 밑 코드 `%s` (고리 %d) · 배경 %d · 밑 코드 자신이 불임인 배경 %d" % (b, loop_ticks(b), len(bgs), dead_self))
    if len(bgs) < max(1, NBG // 2) or dead_self * 2 > len(bgs):
        print("   ⏭️ N0 — 밑 코드가 배경 절반 넘게에서 불임이다. **이 밑 코드는 판정하지 않는다**(사전등록 §4)")
        verdict[b] = None
        OUT["skip|" + b] = {"bg": len(bgs), "dead_self": dead_self}
        continue
    print("   밑 코드 자신 Δ 중앙값 %+.3f" % pct([v[0] for v in sf.values()], .5))
    IX = [[RNG.randrange(len(bgs)) for _ in bgs] for _ in range(N_BOOT)]

    def paired(ka, kb, mu):
        P = PER[(b, mu)]
        g = [s for s in bgs if s in P.get(ka, {}) and s in P.get(kb, {})]
        if len(g) < max(1, NBG // 2):
            return None
        d = [P[ka][s][0] - P[kb][s][0] for s in g]
        m = sum(d) / len(d)
        bs = [sum(d[i % len(d)] for i in ix) / len(ix) for ix in IX]
        return m, pct(bs, .025), pct(bs, .975), len(g), sum(1 for s in g if P[ka][s][1]), sum(1 for s in g if P[kb][s][1])

    res = {"n1": [], "n2": [], "n3": []}
    for ch, tag in (("n", "판정"), ("x", "서술")):
        cand = {k: L for k, L in inserts(b, ch).items() if k in PER[(b, JUDGE_MU)]}
        byL = {}
        for k, L in cand.items(): byL.setdefault(L, []).append(k)
        Ls = sorted(byL)
        print("   +%s (%s) · 고리 %s" % (ch, tag, {L: sorted(byL[L]) for L in Ls}))
        for i in range(len(Ls) - 1):
            lo, hi = Ls[i], Ls[i + 1]
            for ka in sorted(byL[hi]):
                for kb in sorted(byL[lo]):
                    r = paired(ka, kb, JUDGE_MU)
                    if not r: continue
                    m, l, h, n, da, db = r
                    per = m / (hi - lo)
                    ok1 = h < 0; ok2 = PER_TICK_LO <= per <= PER_TICK_HI
                    if ch == "n": res["n1"].append(ok1); res["n2"].append(ok2)
                    print("     ΔΔ %-9s(%d) − %-9s(%d) = %+.3f [%+.3f, %+.3f] · 배경 %d · 틱당 %+.3f %s"
                          % (ka, hi, kb, lo, m, l, h, n, per, ("✅" if (ok1 and ok2) else ("상한≥0" if not ok1 else "틱당 밖")) if ch == "n" else ""))
                    OUT["pair|%s|%s|%s|%s" % (b, ch, ka, kb)] = {"dd": m, "ci": [l, h], "per_tick": per, "bg": n}
        for L in Ls:
            g = sorted(byL[L])
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    r = paired(g[i], g[j], JUDGE_MU)
                    if not r: continue
                    m, l, h, n, da, db = r
                    ok3 = abs(m) <= TOL_NEG
                    if ch == "n": res["n3"].append(ok3)
                    print("     [음성] %-9s vs %-9s (%d틱) = %+.3f [%+.3f, %+.3f] %s"
                          % (g[i], g[j], L, m, l, h, ("✅" if ok3 else "❌ 0.15 밖") if ch == "n" else ""))
                    OUT["neg|%s|%s|%s|%s" % (b, ch, g[i], g[j])] = {"dd": m, "ci": [l, h]}
    n1 = bool(res["n1"]) and all(res["n1"]); n2 = bool(res["n2"]) and all(res["n2"]); n3 = bool(res["n3"]) and all(res["n3"])
    print("   N1 %s(%d짝) · N2 %s · N3 %s(%d쌍)" % ("✅" if n1 else "❌", len(res["n1"]), "✅" if n2 else "❌", "✅" if n3 else "❌", len(res["n3"])))
    verdict[b] = {"N1": n1, "N2": n2, "N3": n3}
    OUT["verdict|" + b] = verdict[b]

    for mu in MUS[1:]:
        print("   N5 서술 — μ %g" % mu)
        cand = {k: L for k, L in inserts(b, "n").items() if k in PER[(b, mu)]}
        byL = {}
        for k, L in cand.items(): byL.setdefault(L, []).append(k)
        Ls = sorted(byL)
        for i in range(len(Ls) - 1):
            lo, hi = Ls[i], Ls[i + 1]
            for ka in sorted(byL[hi])[:2]:
                for kb in sorted(byL[lo])[:2]:
                    r = paired(ka, kb, mu)
                    if r:
                        m, l, h, n, da, db = r
                        print("     ΔΔ %-9s(%d) − %-9s(%d) = %+.3f [%+.3f, %+.3f] · 틱당 %+.3f"
                              % (ka, hi, kb, lo, m, l, h, m / (hi - lo)))
                        OUT["pair_mu|%s|%g|%s|%s" % (b, mu, ka, kb)] = {"dd": m, "ci": [l, h]}

judged = {k: v for k, v in verdict.items() if v}
passed = [k for k, v in judged.items() if v["N1"] and v["N3"]]
print("\n== 종합")
if not judged:
    print("  두 밑 코드 모두 N0 에서 빠졌다 — 판정 불가. 보고서 K2 는 `racld` 가족으로 한정한다")
elif passed:
    print("  ✅ 짝지은 고리 효과가 `racld` 가족 밖에서도 재현된다 — 통과한 밑 코드: %s" % ", ".join(passed))
else:
    print("  ❌ 판정한 밑 코드 어디서도 N1 · N3 이 서지 않는다 — 보고서 K2 를 `racld` 가족으로 한정한다")
OUT["summary"] = {"judged": list(judged), "passed": passed}
tag = os.environ.get("PETRI_OUT", "_RESULT_pairs13")
json.dump(OUT, open(os.path.join(BASE, tag + ".json"), "w"), ensure_ascii=False, indent=1)
print("기록 → %s.json" % tag)
