"""해부 3B 사후 진단 — 판정 아님. 평평함 지표의 약점 하나를 보정해 본다.

  판정(§3)의 평평함은 **서로 다른 이웃 코드에 같은 무게**를 준다. 그런데 sim.js 의 복사 오류는 종류마다 빈도가 다르다:
  베끼는 글자마다 바뀜 μ(글자는 11종 균일) · 끼어듦 μ/3(글자 11종 균일) · 빠뜨림 μ/3.
  → 라벨 하나의 상대 무게: sub(p,글자) 1/11 · del(p) 1/3 · ins(p,글자) 1/33. 코드의 무게 = 그 코드를 내는 라벨 무게의 합.
  (1) 돌연변이 흐름으로 가중한 평평함 · (2) 종류별(바뀜 · 빠뜨림 · 끼어듦) 평균 Δ · (3) 빠뜨림 이웃 하나씩
실행: python3 flat_diag.py [petri 폴더]
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
SETS = {"racld": "dms_main", "rascld": "dms_flat"}
N_BOOT = 2000
RNG = random.Random(20260922)
W = {"sub": 1 / 11, "del": 1 / 3, "ins": 1 / 33}


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


D, LAB = {}, {}
for wt, folder in SETS.items():
    D[wt] = {}
    for f in sorted(glob.glob(os.path.join(BASE, folder, f"dms_mat_{wt}_s*_all.json"))):
        z = json.load(open(f))
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt_arms) / len(wt_arms)
        D[wt][z["seed"]] = {a["key"]: s_of(a) - base for a in z["arms"] if a["kind"] == "mut"}
        LAB[wt] = {a["key"]: a["labels"] for a in z["arms"] if a["kind"] == "mut"}
common = sorted(set(D["racld"]) & set(D["rascld"]))
idxs = [[RNG.randrange(len(common)) for _ in range(len(common))] for _ in range(N_BOOT)]


def ci(vals):
    n = len(vals)
    bm = [sum(vals[i] for i in ix) / n for ix in idxs]
    return sum(vals) / n, pct(bm, .025), pct(bm, .975)


print(f"배경 {len(common)} · 무게 sub {W['sub']:.3f} · del {W['del']:.3f} · ins {W['ins']:.3f}")
weighted = {}
for wt in SETS:
    keys = list(LAB[wt])
    w = {k: sum(W[l.split(':')[0]] for l in LAB[wt][k]) for k in keys}
    tot = sum(w.values())
    per_bg = [sum(w[k] * D[wt][s][k] for k in keys) / tot for s in common]
    weighted[wt] = per_bg
    m, lo, hi = ci(per_bg)
    print(f"(1) {wt:>7}: 돌연변이 흐름 가중 평평함 {m:+.3f} [{lo:+.3f}, {hi:+.3f}] · 이웃 {len(keys)}")
    for kind in ("sub", "del", "ins"):
        ks = [k for k in keys if any(l.startswith(kind) for l in LAB[wt][k])]
        vals = [sum(D[wt][s][k] for s in common) / len(common) for k in ks]
        print(f"(2)    {kind}: 이웃 {len(ks):>3} · 평균 Δ {sum(vals)/len(vals):+.2f} · 중앙 {pct(vals, .5):+.2f}")
d = ci([a - b for a, b in zip(weighted["rascld"], weighted["racld"])])
print(f"(1) 씨앗 − racld (가중): {d[0]:+.3f} [{d[1]:+.3f}, {d[2]:+.3f}]")
for wt in SETS:
    rows = []
    for k, labs in LAB[wt].items():
        for l in labs:
            if l.startswith("del"):
                p = int(l.split(":")[1])
                rows.append((p, k, sum(D[wt][s][k] for s in common) / len(common)))
    rows.sort()
    print(f"(3) {wt:>7} 빠뜨림 이웃: " + " · ".join(f"{p}→{k} {v:+.2f}" for p, k, v in rows))
