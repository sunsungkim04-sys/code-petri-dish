"""파일럿 2 확인 — sim v0.3.0. 판정에 쓰지 않는다. 기능 있는 글자(j · e · c)의 퍼짐은 읽지 않는다.

  (1) 결정론: smoke2/a · smoke2/b 기록이 소요 시간 말고 같은가
  (2) 완료 · 재료 보존 — 계획한 조합마다 시드 20개가 다 모였나(분모는 계획 수)
  (3) 조합별 생존율 · Wilson 95% 구간
  (4) 무늬 x 가 중립 대조로 쓸 만한가 — 단일 코드 퍼짐 L2(가장 큰 단일 코드의 개체 비율을 연속 2 기록 = 500틱 이상 유지한 최댓값)
      자리만 · 재료 조건에서는 h · n 도 '쉬기' 와 같아 x 와 나란히 본다(x 는 틱을 안 써서 h · n 보다 가볍다)
  (5) 50,000틱 소요 시간
실행: python3 pilot2_check.py [petri 폴더, 기본 ~/petri]
"""
import glob
import json
import math
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
PLAN = [("space", 8), ("energy", 8), ("mat", 6), ("mat", 8), ("mat", 10)]
SEEDS = 20
MIN_TICKS = 5000


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def sustained_max(vals, k):
    if len(vals) < k:
        return 0.0
    return max(min(vals[i:i + k]) for i in range(len(vals) - k + 1))


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))] if xs else float("nan")


print("== (1) 결정론")
sa = glob.glob(os.path.join(BASE, "smoke2/a/*.json"))
sb = glob.glob(os.path.join(BASE, "smoke2/b/*.json"))
if len(sa) == 1 and len(sb) == 1:
    a, b = json.load(open(sa[0])), json.load(open(sb[0]))
    a.pop("elapsed_ms"); b.pop("elapsed_ms")
    print(f"  같은 시드 두 번 기록 동일: {a == b} · 체크섬 {a['checksum']} / {b['checksum']} · sim {a['sim_version']}")
else:
    print(f"  smoke2 파일 수 이상: a {len(sa)} · b {len(sb)}")

groups = {}
for f in glob.glob(os.path.join(BASE, "pilot2/*.json")):
    z = json.load(open(f))
    groups.setdefault((z["cond"], int(z["opts"]["density"])), []).append(z)

print("\n== (2)(3) 완료 · 재료 보존 · 생존율")
for key in PLAN:
    zs = groups.get(key, [])
    alive = sum(z["extinct_at"] < 0 for z in zs)
    lo, hi = wilson(alive, len(zs))
    versions = sorted({z["sim_version"] for z in zs})
    cons = sum(z["conservation_violations"] for z in zs)
    checks = sum(z["conservation_checks"] for z in zs)
    print(f"  {key[0]:>6} d{key[1]:<2} 완료 {len(zs)}/{SEEDS}{'  ← 미완' if len(zs) != SEEDS else ''} · 생존 {alive}/{len(zs)} "
          f"[{lo:.0%}, {hi:.0%}] · 보존 위반 {cons}/{checks} · sim {versions} · 최대 {max((z['elapsed_ms'] for z in zs), default=0) / 1000:.0f}초")

print("\n== (4) 중립 대조 — 단일 코드 퍼짐 L2 (5,000틱 전 멸종 접시 제외)")
print(f"  {'조건':>6} {'밀도':>4} {'접시':>4} | {'글자':>2} {'중앙':>6} {'95%':>6} {'최대':>6}")
for key in PLAN:
    zs = [z for z in groups.get(key, []) if not (0 <= z["extinct_at"] < MIN_TICKS)]
    letters = ("x",) if key[0] == "energy" else ("x", "h", "n")
    for letter in letters:
        vals = []
        for z in zs:
            pop_at = {r[0]: r[1] for r in z["samples"]}
            col = z["special_header"].index(letter) + 1
            shares = [(row[col][1] / pop_at[row[0]]) if (row[col] and pop_at.get(row[0])) else 0.0 for row in z["specials"]]
            vals.append(sustained_max(shares, 2))
        if vals:
            print(f"  {key[0]:>6} {key[1]:>4} {len(zs):>4} | {letter:>2} {q(vals, .5):>6.1%} {q(vals, .95):>6.1%} {max(vals):>6.1%}")
