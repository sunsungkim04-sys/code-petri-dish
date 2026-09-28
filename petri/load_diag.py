"""해부 3B 사후 진단 2 — 판정 아님. '표적 크기(길이)' 후보를 이미 있는 자료로 검사한다.

  복제 오류가 있으면 자식을 만들 때 베끼는 글자가 많을수록 돌연변이를 맞을 기회가 많다(부하 ∝ 길이).
  길이만 한 글자 늘리고 기능은 그대로인 이웃 = 무늬 `x` 끼어듦. 이 이웃의 Δ 가 μ 0 → μ 1% 에서 얼마나 나빠지나.
  대조로 `n`(쉬기 · 틱은 쓰지만 기능 없음) 끼어듦, 그리고 길이가 줄어드는 빠뜨림 이웃도 같이 본다.
  자료: 해부 1(dms_main · μ 0)과 해부 2A(dms_mu1 · μ 1%) — 같은 배경 · 같은 팔. 새로 돌리는 것 없음.
실행: python3 load_diag.py [petri 폴더]
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
WTS = [("space", "reasccld"), ("mat", "racld"), ("energy", "reascld")]
N_BOOT = 2000
RNG = random.Random(20260923)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(a):
    r = a["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def per_bg(folder, cond, wt):
    out, labs = {}, {}
    for f in sorted(glob.glob(os.path.join(BASE, folder, f"dms_{cond}_{wt}_s*_all.json"))):
        z = json.load(open(f))
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt_arms) / len(wt_arms)
        out[z["seed"]] = {a["key"]: s_of(a) - base for a in z["arms"] if a["kind"] == "mut"}
        labs = {a["key"]: a["labels"] for a in z["arms"] if a["kind"] == "mut"}
    return out, labs


for cond, wt in WTS:
    d0, labs = per_bg("dms_main", cond, wt)
    d1, _ = per_bg("dms_mu1", cond, wt)
    seeds = sorted(set(d0) & set(d1))
    idxs = [[RNG.randrange(len(seeds)) for _ in range(len(seeds))] for _ in range(N_BOOT)]

    def ci(vals):
        n = len(vals)
        bm = [sum(vals[i] for i in ix) / n for ix in idxs]
        return sum(vals) / n, pct(bm, .025), pct(bm, .975)

    def group(pred):
        ks = [k for k, L in labs.items() if any(pred(l) for l in L)]
        if not ks:
            return None
        v0 = [sum(d0[s][k] for k in ks) / len(ks) for s in seeds]
        v1 = [sum(d1[s][k] for k in ks) / len(ks) for s in seeds]
        dd = [b - a for a, b in zip(v0, v1)]
        return len(ks), ci(v0), ci(v1), ci(dd)

    print(f"\n==================== {cond} · WT {wt}(길이 {len(wt)}) · 배경 {len(seeds)}")
    for name, pred in (("무늬 x 끼어듦(길이 +1 · 기능 그대로)", lambda l: l.startswith("ins:") and l.endswith(":x")),
                       ("쉬기 n 끼어듦(길이 +1 · 틱 씀)", lambda l: l.startswith("ins:") and l.endswith(":n")),
                       ("끼어듦 전부(길이 +1)", lambda l: l.startswith("ins:")),
                       ("바뀜 전부(길이 그대로)", lambda l: l.startswith("sub:")),
                       ("빠뜨림 전부(길이 −1)", lambda l: l.startswith("del:"))):
        g = group(pred)
        if not g:
            continue
        n, a, b, dd = g
        print(f"  {name:<28} 이웃 {n:>3} · Δ μ0 {a[0]:+.2f} [{a[1]:+.2f},{a[2]:+.2f}] · Δ μ1 {b[0]:+.2f} [{b[1]:+.2f},{b[2]:+.2f}] · "
              f"μ1−μ0 {dd[0]:+.2f} [{dd[1]:+.2f},{dd[2]:+.2f}]")
