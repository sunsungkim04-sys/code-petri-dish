"""되감기 2 설계용 탐색 — 본실험(시드 1~100) 기록으로 **측정이 되는지만** 본다. 판정 아님 · 순서의 답(τ)은 계산하지 않는다.

본실험 §8 의 설계 결함 둘을 고칠 후보를 잰다.
  Q2 기준자가 너무 촘촘했다 → 후보 기준자들이 포화되는지
    (a) 글자 구성 벡터(유전체당 글자 11종 평균 개수, 마지막 1,000틱 평균)의 유클리드 거리
        접시 안 = 같은 접시의 [49k,50k] 대 [39k,40k] · [24k,25k] / 접시 사이 = 끝끼리
    (b) 상위 20 코드 조성(개체 비율)의 Bray–Curtis — 접시 안 = 끝 대 40,000 · 25,000 / 사이 = 끝끼리
    (c) 우세 코드 정체의 접시 안 안정성 — 마지막 10,000틱 상위 기록 21개 중 끝 우세 코드와 같은 비율
  Q3 사건이 접시당 1~2개였다 → 사건 후보들이 접시당 몇 개 나오는지 · 접시 쌍이 공통 사건을 몇 개 갖는지
    (L) 글자 구성 사건: 글자마다 유전체당 평균 개수가 씨앗 값 ± 0.5 를 넘어 500틱(10 샘플) 유지
        개체 ≥ 0 · ≥ 1,000 두 판 · 무늬 x 의 사건률은 중립 대조
    (S) 우세 교대 사건: 씨앗이 아닌 코드가 상위 1위를 연속 2 기록(1,000틱) 차지한 첫 시각
실행: python3 rewind2_explore.py [petri 폴더, 기본 ~/petri]
"""
import glob
import json
import math
import os
import sys
from collections import Counter

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
LET = "nsracldjehx"
SEEDV = {ch: "rascld".count(ch) for ch in LET}
CONDS = ["space", "energy", "mat"]
HOLD = 10


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def q(xs):
    return f"중앙 {pct(xs, .5):.3f} [5% {pct(xs, .05):.3f} · 95% {pct(xs, .95):.3f}]"


data = {c: [] for c in CONDS}
for f in sorted(glob.glob(os.path.join(BASE, "main", "*.json"))):
    z = json.load(open(f))
    if z["extinct_at"] < 0 and z["opts"]["seed"] <= 100:
        data[z["cond"]].append(z)


def letter_vec(z, t_end):
    idx = [z["header"].index("n_" + ch) for ch in LET]
    rows = [r for r in z["samples"] if t_end - 1000 < r[0] <= t_end and r[1] > 0]
    return [sum(r[i] / r[1] for r in rows) / len(rows) for i in idx]


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def top_at(z, t):
    for tt, top in z["tops"]:
        if tt == t:
            return top
    return None


def pop_at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r[1]
    return None


def comp(z, t):
    top, pop = top_at(z, t), pop_at(z, t)
    return {k: c / pop for k, c in top}


def bray(a, b):
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.get(k, 0) + b.get(k, 0) for k in keys)
    return num / den if den else float("nan")


def letter_events(z, min_pop):
    H = z["header"]
    ev = {}
    for ch in LET:
        i = H.index("n_" + ch)
        for sign in ("+", "-"):
            if sign == "-" and SEEDV[ch] == 0:
                continue
            run, start = 0, None
            for r in z["samples"]:
                if r[1] <= 0 or r[1] < min_pop:
                    run = 0
                    continue
                v = r[i] / r[1]
                hit = v >= SEEDV[ch] + 0.5 if sign == "+" else v <= SEEDV[ch] - 0.5
                if hit:
                    if run == 0:
                        start = r[0]
                    run += 1
                    if run == HOLD:
                        ev[ch + sign] = start
                        break
                else:
                    run = 0
    return ev


def succession(z):
    ev, prev, run = {}, None, 0
    for t, top in z["tops"]:
        if not top:
            prev, run = None, 0
            continue
        k = top[0][0]
        run = run + 1 if k == prev else 1
        prev = k
        if run == 2 and k != "rascld" and k not in ev:
            ev[k] = t - 500
    return ev


def pair_common(evs, need=3):
    n = len(evs)
    pairs = hit = 0
    commons = []
    for i in range(n):
        for j in range(i + 1, n):
            c = len(set(evs[i]) & set(evs[j]))
            commons.append(c)
            pairs += 1
            hit += c >= need
    return pairs, hit, commons


for c in CONDS:
    zs = data[c]
    n = len(zs)
    print(f"\n==================== {c} · 생존 {n}")

    cov = [sum(cnt for _, cnt in top_at(z, 50000)) / pop_at(z, 50000) for z in zs]
    print(f"상위 20 코드가 끝 개체에서 차지하는 비율: {q(cov)}")

    ends = [letter_vec(z, 50000) for z in zs]
    w10 = [dist(letter_vec(z, 50000), letter_vec(z, 40000)) for z in zs]
    w25 = [dist(letter_vec(z, 50000), letter_vec(z, 25000)) for z in zs]
    btw = [dist(ends[i], ends[j]) for i in range(n) for j in range(i + 1, n)]
    print("(a) 글자 구성 거리")
    print(f"    접시 안 1만 틱: {q(w10)}")
    print(f"    접시 안 2.5만 틱: {q(w25)}")
    print(f"    접시 사이(끝끼리): {q(btw)}")

    cend = [comp(z, 50000) for z in zs]
    b10 = [bray(comp(z, 50000), comp(z, 40000)) for z in zs]
    b25 = [bray(comp(z, 50000), comp(z, 25000)) for z in zs]
    bb = [bray(cend[i], cend[j]) for i in range(n) for j in range(i + 1, n)]
    print("(b) 상위 20 조성 Bray–Curtis")
    print(f"    접시 안 1만 틱: {q(b10)}")
    print(f"    접시 안 2.5만 틱: {q(b25)}")
    print(f"    접시 사이(끝끼리): {q(bb)}")

    stab = []
    flips = []
    for z in zs:
        fin = z["final_counts"][0][0]
        recs = [top[0][0] for t, top in z["tops"] if t >= 40000 and top]
        stab.append(sum(k == fin for k in recs) / len(recs))
        flips.append(len(set(recs)))
    print(f"(c) 마지막 1만 틱 상위 기록 중 끝 우세 코드와 같은 비율: {q(stab)} · 서로 다른 1위 코드 수 분포 {dict(sorted(Counter(flips).items()))}")

    for mp in (0, 1000):
        evs = [letter_events(z, mp) for z in zs]
        cnt = Counter(len(e) for e in evs)
        rate = Counter(k for e in evs for k in e)
        pairs, hit, commons = pair_common(evs)
        print(f"(L) 글자 구성 사건 · 개체 ≥ {mp:,}: 접시당 개수 {dict(sorted(cnt.items()))}")
        print("    사건별 발생 접시: " + " · ".join(f"{k} {v}" for k, v in sorted(rate.items(), key=lambda kv: -kv[1])))
        print(f"    공통 사건 ≥ 3 인 접시 쌍 {hit}/{pairs} · 공통 사건 수 분포 {dict(sorted(Counter(commons).items()))}")
    evs = [succession(z) for z in zs]
    cnt = Counter(len(e) for e in evs)
    rate = Counter(k for e in evs for k in e)
    pairs, hit, commons = pair_common(evs)
    print(f"(S) 우세 교대 사건: 접시당 개수 {dict(sorted(cnt.items()))}")
    print("    흔한 교대 코드: " + " · ".join(f"{k} {v}" for k, v in rate.most_common(8)))
    print(f"    공통 사건 ≥ 3 인 접시 쌍 {hit}/{pairs} · 공통 사건 수 분포 {dict(sorted(Counter(commons).items()))}")
