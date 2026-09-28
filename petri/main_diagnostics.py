"""본실험 사후 진단 — 사전등록 밖, 판정에 쓰지 않는다. main_analyze.py 의 판정을 해석하려고 본 것.

  (1) 조건별 끝 우세 코드 상위 5 · 우세 코드의 끝 개체 비율
  (2) 접시당 일어난 사건 개수 분포 — Q3 에서 공통 사건 ≥ 3 인 쌍이 왜 적었나
  (3) 중립 문턱 θ 의 출처 — L2_x 최대 창이 몇 틱에 있었나 · 그때 개체 수
      + 🟡 사후 민감도: 개체 ≥ 1,000 인 기록만으로 L2 를 다시 재면 θ 와 사건 발생률이 어떻게 되나 (판정 아님)
  (4) 가장 큰 e 코드 · h 코드 비율의 시간 경과(생존 접시 중앙) · 특정 틱에 가장 흔한 '가장 큰 e 코드'
실행: python3 main_diagnostics.py [petri 폴더, 기본 ~/petri]  (먼저 main_analyze.py 를 돌려 _RESULT_main.json 이 있어야 한다)
"""
import glob
import json
import math
import os
import sys
from collections import Counter

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
RES = json.load(open(os.path.join(BASE, "_RESULT_main.json")))
CONDS = ["space", "energy", "mat"]
MIN_POP = 1000
TICKS_SHOW = [500, 1000, 1500, 2000, 2500, 3000, 5000, 10000, 25000, 50000]


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


data = {c: [] for c in CONDS}
for f in glob.glob(os.path.join(BASE, "main", "*.json")):
    z = json.load(open(f))
    if z["extinct_at"] < 0:
        data[z["cond"]].append(z)


def series(z, tag, min_pop=0):
    pop_at = {r[0]: r[1] for r in z["samples"]}
    col = z["special_header"].index(tag) + 1
    out = []
    for row in z["specials"]:
        pop = pop_at.get(row[0], 0)
        share = (row[col][1] / pop) if (row[col] and pop) else 0.0
        out.append((row[0], share if pop >= min_pop else 0.0, pop, row[col][0] if row[col] else None))
    return out


def l2_window(s):
    best, at = -1.0, None
    for i in range(len(s) - 1):
        v = min(s[i][1], s[i + 1][1])
        if v > best:
            best, at = v, i
    return best, (s[at][0] if at is not None else None), (s[at][2] if at is not None else None)


for c in CONDS:
    zs = data[c]
    n = len(zs)
    theta = RES["conds"][c]["theta"]
    print(f"\n==================== {c} · 생존 {n} · θ(사전등록) {theta:.2%}")

    doms = Counter(z["final_counts"][0][0] for z in zs)
    dom_share = [z["final_counts"][0][1] / sum(cnt for _, cnt in z["final_counts"]) for z in zs]
    print("(1) 끝 우세 코드 상위 5:", " · ".join(f"{k} {v}" for k, v in doms.most_common(5)))
    print(f"    우세 코드의 끝 개체 비율 중앙 {pct(dom_share, .5):.1%} [{pct(dom_share, .25):.1%}, {pct(dom_share, .75):.1%}] · 끝 개체 수 중앙 {pct([sum(cnt for _, cnt in z['final_counts']) for z in zs], .5):,.0f}")

    ev_n = Counter(RES["q3"][c]["with_seed"]["n_pairs_tau"] for _ in [0])  # 자리 표시 (아래에서 직접 센다)
    events_per_dish = []
    for z in zs:
        k = 0
        for tag in ("j", "noc", "e", "h"):
            if tag == "h" and c in ("space", "mat"):
                continue
            s = series(z, tag)
            if any(min(s[i][1], s[i + 1][1]) > theta for i in range(len(s) - 1)):
                k += 1
        k += 1 if any(True for _ in [0]) and any("rascld" not in {kk for kk, _ in top} for _, top in z["tops"]) else 0
        events_per_dish.append(k)
    print("(2) 접시당 사건 수(E_short 제외, E_seed 포함):", dict(sorted(Counter(events_per_dish).items())))

    wins = [l2_window(series(z, "x")) for z in zs]
    ticks = [t for _, t, _ in wins if t is not None]
    pops = [p for _, _, p in wins if p is not None]
    early = sum(t < 1000 for t in ticks)
    print(f"(3) L2_x 최대 창의 틱 중앙 {pct(ticks, .5):,.0f} · 1,000틱 이전 {early}/{len(ticks)} · 그때 개체 수 중앙 {pct(pops, .5):,.0f}")
    l2x_big = [l2_window(series(z, "x", MIN_POP))[0] for z in zs]
    theta2 = pct(l2x_big, .95)
    line = []
    for tag in ("j", "noc", "e", "h"):
        hits = sum(l2_window(series(z, tag, MIN_POP))[0] > theta2 for z in zs)
        line.append(f"{tag} {hits}/{n}")
    print(f"    🟡 사후 민감도(개체 ≥ {MIN_POP:,} 인 기록만): θ' {theta2:.2%} · 문턱 넘은 접시 " + " · ".join(line))

    for tag in ("e", "h"):
        byt = {t: [] for t in TICKS_SHOW}
        keys_at = {t: Counter() for t in (2500, 50000)}
        for z in zs:
            for tick, share, pop, key in series(z, tag):
                if tick in byt:
                    byt[tick].append(share)
                if tick in keys_at and key:
                    keys_at[tick][key] += 1
        print(f"(4) 가장 큰 {tag} 코드 비율 중앙: " + " · ".join(f"{t:,} {pct(v, .5):.1%}" for t, v in byt.items() if v))
        for t, cnt in keys_at.items():
            print(f"    {t:,}틱에 가장 흔한 '가장 큰 {tag} 코드': " + " · ".join(f"{k} {v}" for k, v in cnt.most_common(3)))
