"""파일럿 — 기능 없는 글자가 돌연변이 부하와 표류만으로 어디까지 퍼지나 (사건 문턱 보정용).

사전등록 초안의 내부 대조가 걸렸다: 자리만 조건에서 기능 없는 h 가 20접시 중 19개에서 200틱 만에 1% 를 넘었다.
그래서 사건 문턱을 **중립 글자로만** 다시 정한다. 기능 있는 글자(j · e · c)는 이 스크립트가 읽지 않는다.

중립 글자 = h 와 n. 자리만(space) · 재료(mat) 조건에서는 둘 다 '쉬기' 와 같다(에너지 조건에서는 h 가 기능을 갖고 n 은 비용이 싸서 제외).
재는 것 (접시마다):
  L1_h — h 가 든 코드 전체의 개체 비율, 연속 10 샘플(500틱) 유지된 최댓값
  L2_h · L2_n — h(n) 가 든 **단일 코드** 하나의 개체 비율, 상위 20 기록 연속 2회(500틱) 유지된 최댓값
             (상위 20 에 없으면 0 으로 센다 — 최댓값을 과소가 아니라 '20위 아래' 로만 모르는 것)
문턱 규칙(사용자 승인 전 제안): 후보 {1, 2, 5, 10, 20, 30, 50}% 중 중립 대조가 도달한 접시가 5% 이하인 가장 작은 값.
실행: python3 pilot_neutral_calibration.py [결과 폴더, 기본 ~/petri/pilot]
"""
import glob
import json
import os
import sys

DIR = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri/pilot")
MIN_TICKS = 5000          # 이보다 일찍 멸종한 접시는 표류를 볼 시간이 없어 뺀다(개수 인쇄)
CANDIDATES = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]


def sustained_max(vals, k):
    if len(vals) < k:
        return 0.0
    return max(min(vals[i:i + k]) for i in range(len(vals) - k + 1))


rows = []
skipped = 0
for f in sorted(glob.glob(os.path.join(DIR, "*.json"))):
    z = json.load(open(f))
    if z["cond"] not in ("space", "mat"):
        continue
    if 0 <= z["extinct_at"] < MIN_TICKS:
        skipped += 1
        continue
    ix = {h: i for i, h in enumerate(z["header"])}
    pop_at = {r[0]: r[ix["pop"]] for r in z["samples"]}
    shares_h = [r[ix["has_h"]] / r[ix["pop"]] if r[ix["pop"]] else 0.0 for r in z["samples"]]

    def single(letter):
        out = []
        for tick, top in z["tops"]:
            pop = pop_at.get(tick, 0)
            best = max((cnt for key, cnt in top if letter in key), default=0)
            out.append(best / pop if pop else 0.0)
        return out

    rows.append(dict(
        cond=z["cond"], density=int(z["opts"]["density"]), seed=z["opts"]["seed"],
        L1_h=sustained_max(shares_h, 10), L2_h=sustained_max(single("h"), 2), L2_n=sustained_max(single("n"), 2)))


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))]


print(f"폴더 {DIR} · 대상 접시 {len(rows)} (자리만·재료, {MIN_TICKS}틱 전 멸종 {skipped}개 제외)")
groups = [("space", None)] + [("mat", d) for d in (4, 6, 8, 10, 14)] + [("전체", "pool")]
print(f"{'조건':>6} {'밀도':>4} {'접시':>4} | {'L1_h 중앙':>9} {'90%':>6} {'최대':>6} | {'L2_h 중앙':>9} {'90%':>6} {'최대':>6} | {'L2_n 중앙':>9} {'90%':>6} {'최대':>6}")
for cond, d in groups:
    rs = rows if d == "pool" else [r for r in rows if r["cond"] == cond and (d is None or r["density"] == d)]
    if not rs:
        print(f"{cond:>6} {str(d):>4} {0:>4} | (없음)")
        continue
    cells = []
    for key in ("L1_h", "L2_h", "L2_n"):
        xs = [r[key] for r in rs]
        cells.append(f"{q(xs, .5):>9.1%} {q(xs, .9):>6.1%} {max(xs):>6.1%}")
    print(f"{cond:>6} {str(d if d is not None else '-'):>4} {len(rs):>4} | " + " | ".join(cells))

print("\n문턱 후보별 — 중립 대조가 그 비율에 도달한 접시 수 (전체 풀)")
for key in ("L1_h", "L2_h", "L2_n", "L2_hn"):
    xs = [max(r["L2_h"], r["L2_n"]) for r in rows] if key == "L2_hn" else [r[key] for r in rows]
    line = " · ".join(f"{c:.0%} {sum(x >= c for x in xs)}/{len(xs)}" for c in CANDIDATES)
    pick = next((c for c in CANDIDATES if sum(x >= c for x in xs) <= 0.05 * len(xs)), None)
    print(f"  {key:>6}: {line}  →  제안 문턱 {('%.0f%%' % (pick * 100)) if pick is not None else '후보 안에 없음'}")

json.dump(rows, open(os.path.join(DIR, "..", "_RESULT_pilot_neutral.json"), "w"), indent=1)
