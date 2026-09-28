"""해부 3A 사후 진단 — 판정 아님. '더 나은 이웃이 왜 안 퍼지나' 를 기록으로 되짚는다.

  끝 우세 코드는 *가장 흔한 코드 하나* 다. 그 글자를 가진 코드가 여러 변형으로 흩어져 있으면 합은 클 수 있다.
  replay.js 의 samples 에는 그 글자를 **가진 개체 수**(has_h · has_e · has_x · no_c)가 그대로 있다 → 합쳐서 본다.
  본실험(시드 1~100) 생존 접시 · 조건마다 중앙값.
실행: python3 freq_diag.py [petri 폴더]
"""
import glob
import json
import math
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
CONDS = ["space", "energy", "mat"]
TICKS = [500, 1000, 2000, 3000, 5000, 10000, 25000, 50000]


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


data = {c: [] for c in CONDS}
for f in sorted(glob.glob(os.path.join(BASE, "main", "*.json"))):
    z = json.load(open(f))
    if z["extinct_at"] < 0:
        data[z["cond"]].append(z)

for c in CONDS:
    zs = data[c]
    H = zs[0]["header"]
    ih, ie, ix = H.index("has_h"), H.index("has_e"), H.index("has_x")
    print(f"\n== {c} · 생존 {len(zs)} — 그 글자를 가진 개체의 비율(중앙값)")
    for t in TICKS:
        h, e, x = [], [], []
        for z in zs:
            r = next((r for r in z["samples"] if r[0] == t), None)
            if not r or not r[1]:
                continue
            h.append(r[ih] / r[1]); e.append(r[ie] / r[1]); x.append(r[ix] / r[1])
        print(f"   틱 {t:>6,}: h {pct(h, .5):5.1%} · e {pct(e, .5):5.1%} · x(무늬 · 중립 대조) {pct(x, .5):5.1%}")
    fin = []
    for z in zs:
        tot = sum(cnt for _, cnt in z["final_counts"])
        fin.append(sum(cnt for k, cnt in z["final_counts"] if "h" in k) / tot)
    print(f"   끝 코드 전수에서 h 든 코드 비율 중앙 {pct(fin, .5):.1%} · 5% {pct(fin, .05):.1%} · 95% {pct(fin, .95):.1%}")
    kinds = []
    for z in zs:
        kinds.append(sum(1 for k, _ in z["final_counts"] if "h" in k))
    print(f"   끝에 h 든 서로 다른 코드 수 중앙 {pct(kinds, .5):,.0f}")
