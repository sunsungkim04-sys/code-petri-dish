"""되감기 2 사후 진단 — 사전등록 밖, 판정에 쓰지 않는다. 판정 결과를 읽으려고 본 것.

  (1) Q1 재현이 안 된 사건의 문턱 민감도 — 본실험(시드 1~100)과 되감기 2(101~200)의 조건별 θ 를 서로 바꿔 끼우면
      사건 발생 수가 어떻게 되나. L2 계산은 main_diagnostics.py 와 같다(가장 큰 단일 코드 비율의 연속 2 기록 최솟값의 최댓값).
  (2) 순서 층을 받친 코드 쌍 — 파생 쌍(한 글자 편집)이 아닌 교대 코드 쌍마다: 둘 다 거친 접시 수 · 그중 u 가 먼저인 접시 수
  (3) 길 층 — 조건별 흔한 길 상위 5
실행: python3 rewind2_diagnostics.py [petri 폴더, 기본 ~/petri]  (먼저 main_analyze.py 두 번과 rewind2_analyze.py 를 돌려 _RESULT_*.json 이 있어야 한다)
"""
import glob
import json
import os
import sys
from collections import Counter
from itertools import combinations

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
CONDS = ["space", "energy", "mat"]
SETS = {"본실험 1~100": ("main", "_RESULT_main.json"), "되감기2 101~200": ("rewind2", "_RESULT_rewind2_q1.json")}
CHECK = [("space", "j"), ("energy", "e"), ("energy", "h")]
SEED_CODE = "rascld"


def one_edit(a, b):
    if a == b:
        return True
    if abs(len(a) - len(b)) > 1:
        return False
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) == 1
    short, long_ = (a, b) if len(a) < len(b) else (b, a)
    return any(short == long_[:i] + long_[i + 1:] for i in range(len(long_)))


def l2(z, tag):
    pop_at = {r[0]: r[1] for r in z["samples"]}
    col = z["special_header"].index(tag) + 1
    s = []
    for row in z["specials"]:
        pop = pop_at.get(row[0], 0)
        s.append((row[col][1] / pop) if (row[col] and pop) else 0.0)
    return max((min(s[i], s[i + 1]) for i in range(len(s) - 1)), default=0.0)


def path(z):
    ev, order, prev, prev_t, run = {}, [], None, None, 0
    for t, top in z["tops"]:
        k = top[0][0] if top else None
        if k is None:
            prev, run = None, 0
            continue
        run = run + 1 if k == prev else 1
        if run == 2 and k != SEED_CODE and k not in ev:
            ev[k] = prev_t
            order.append(k)
        prev, prev_t = k, t
    return tuple(order), ev


surv, theta = {}, {}
for name, (folder, res) in SETS.items():
    R = json.load(open(os.path.join(BASE, res)))
    theta[name] = {c: R["conds"][c]["theta"] for c in CONDS}
    surv[name] = {c: [] for c in CONDS}
    for f in sorted(glob.glob(os.path.join(BASE, folder, "*.json"))):
        z = json.load(open(f))
        if z["extinct_at"] < 0:
            surv[name][z["cond"]].append(z)

print("(1) θ 바꿔 끼우기 — 사건 발생 접시 수")
for c, tag in CHECK:
    for name in SETS:
        vals = [l2(z, tag) for z in surv[name][c]]
        line = " · ".join(f"θ({other.split()[0]}) {theta[other][c]:.2%} → {sum(v > theta[other][c] for v in vals)}/{len(vals)}" for other in SETS)
        print(f"  {c:>6} E_{tag} · 자료 {name}: {line}")

for name in SETS:
    print(f"\n==================== 자료 {name}")
    for c in CONDS:
        paths = [path(z) for z in surv[name][c]]
        pc = Counter(p for p, _ in paths)
        print(f"(3) {c} 흔한 길: " + " · ".join(f"[{' → '.join(p) if p else '교대 없음'}] {k}" for p, k in pc.most_common(5)))
        pair_n, pair_first = Counter(), Counter()
        for _, ev in paths:
            for u, v in combinations(sorted(ev), 2):
                if one_edit(u, v) or ev[u] == ev[v]:
                    continue
                pair_n[(u, v)] += 1
                pair_first[(u, v)] += ev[u] < ev[v]
        top = pair_n.most_common(6)
        print(f"(2) {c} 파생 아닌 교대 쌍(둘 다 거친 접시 · 앞선 쪽): " + (" · ".join(
            f"{(u if pair_first[(u, v)] * 2 >= n else v)} 먼저 {max(pair_first[(u, v)], n - pair_first[(u, v)])}/{n} ({u}·{v})" for (u, v), n in top) if top else "없음"))
