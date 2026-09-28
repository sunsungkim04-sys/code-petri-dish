"""파일럿 요약 — 조합(조건 × 밀도)마다 멸종률 · 정착 시점 · 끝 상태 · 특수 글자 계통이 처음 1% 를 넘는 시각.

계획한 조합마다 시드 20개가 다 모였는지 먼저 확인하고, 모자라면 그 조합은 '미완' 으로 표시한다(분모는 계획 수).
실행: python3 pilot_summary.py [결과 폴더, 기본 ~/petri/pilot]
"""
import glob
import json
import os
import statistics as st
import sys
from collections import defaultdict

DIR = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri/pilot")
PLANNED_SEEDS = 20
PLAN = [("space", 4), ("energy", 4)] + [(c, d) for d in (4, 6, 8, 10, 14) for c in ("mat", "both")]


def med(xs):
    xs = [x for x in xs if x is not None]
    return f"{st.median(xs):,.0f}" if xs else "—"


rows = defaultdict(list)
files = glob.glob(os.path.join(DIR, "*.json"))
for f in files:
    z = json.load(open(f))
    dens = int(z["opts"]["density"]) if z["opts"]["mat"] else 4
    ix = {h: i for i, h in enumerate(z["header"])}
    S = z["samples"]

    def first_share(col, thr=0.01):
        for r in S:
            if r[ix["pop"]] > 0 and r[ix[col]] / r[ix["pop"]] >= thr:
                return r[0]
        return None

    mx = max(r[ix["pop"]] for r in S)
    plateau = next((r[0] for r in S if r[ix["pop"]] >= 0.9 * mx), None) if mx else None
    last = S[-1]
    rows[(z["cond"], dens)].append(dict(
        ext=z["extinct_at"], pop=last[ix["pop"]], kinds=last[ix["kinds"]], length=last[ix["mean_len"]],
        maxpop=mx, plateau=plateau, j=first_share("has_j"), e=first_share("has_e"), h=first_share("has_h"),
        cons=z["conservation_violations"], checks=z["conservation_checks"], sec=z["elapsed_ms"] / 1000))

print(f"파일 {len(files)}개 · 계획 {len(PLAN) * PLANNED_SEEDS}개 · 폴더 {DIR}")
print(f"{'조건':>6} {'밀도':>4} {'완료':>5} {'멸종':>6} {'멸종틱 중앙':>10} {'생존 끝개체':>10} {'코드종류':>8} {'평균길이':>8} "
      f"{'정착틱':>7} {'j≥1%':>12} {'e≥1%':>12} {'h≥1%':>12} {'보존위반':>8} {'최대초':>6}")
for key in PLAN:
    rs = rows.get(key, [])
    done = len(rs)
    tag = "" if done == PLANNED_SEEDS else "  ← 미완"
    ext = [r for r in rs if r["ext"] >= 0]
    alive = [r for r in rs if r["ext"] < 0]

    def reach(col):
        hit = [r[col] for r in rs if r[col] is not None]
        return f"{len(hit)}/{done} {med(hit)}"

    mlen = "{:.2f}".format(st.median([r["length"] for r in alive])) if alive else "—"
    print(f"{key[0]:>6} {key[1]:>4} {done:>2}/{PLANNED_SEEDS} {len(ext):>3}/{done:<2} {med([r['ext'] for r in ext]):>10} "
          f"{med([r['pop'] for r in alive]):>10} {med([r['kinds'] for r in alive]):>8} {mlen:>8} "
          f"{med([r['plateau'] for r in alive]):>7} {reach('j'):>12} {reach('e'):>12} {reach('h'):>12} "
          f"{sum(r['cons'] for r in rs):>4}/{sum(r['checks'] for r in rs):<5} {max((r['sec'] for r in rs), default=0):>6.0f}{tag}")
