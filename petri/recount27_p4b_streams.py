#!/usr/bin/env python3
"""해부27 감사 — P4b 상한의 재추출 흐름 민감도(판정 아님 · 판정선 불변).
mini27_out 원자료에서 loc · roll · ρ8 · μ0 · 끼운 팔의 ln[S(2)/S(4)] 를 배경마다 다시 세고
(정의 = prereg27_frozen.md §4 P4b · P4 자격), 서로 다른 sha256 흐름 N 개로 2,000회 재추출해
95% 상한이 ln 0.8 아래(❌)로 가는 비율을 센다. 쓰기: _RECOUNT27_p4b_streams.txt 만.
usage: /usr/bin/python3 recount27_p4b_streams.py  (petri/ 에서)"""
import json, os, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "mini27_out")
NB, NS = 2000, 500
cands = sorted(int(f.split("_s")[1].split(".")[0]) for f in os.listdir(OUT) if f.startswith("mini27_r8_s") and f.endswith(".json"))
bg = []
for sd in cands:
    d = json.load(open(os.path.join(OUT, "mini27_r8_s%d.json" % sd)))
    if all(d["by_geom"][g]["grow"]["qualified"] for g in ("glob", "loc", "place")):
        bg.append((sd, d))
    if len(bg) == 20: break
lS = []
for sd, d in bg:
    a = [x for x in d["by_geom"]["loc"]["arms"] if x["rule"] == "roll" and x["mu"] == 0.0 and x["kind"] == "mut"][0]
    c = a["counters"]
    p2, p4 = c["1_2"]["positions"], c["0_4"]["positions"]
    R2, R4 = c["1_2"]["draws"] / p2, c["0_4"]["draws"] / p4
    if p2 >= 100 and p4 >= 100 and R2 > 1 and R4 > 1:
        lS.append(math.log((R2 - 1) * 2) - math.log((R4 - 1) * 4))
lS = np.array(lS); a08 = math.log(0.8)
his = []
for k in range(NS):
    h = hashlib.sha256(("audit27|P4b|%d" % k).encode()).digest()
    r = np.random.default_rng(int.from_bytes(h[:8], "little"))
    bm = lS[r.integers(0, len(lS), size=(NB, len(lS)))].mean(1)
    his.append(np.percentile(bm, 97.5))
his = np.array(his)
L = ["배경 %s (자격 %d)" % (" ".join(str(s) for s, _ in bg), len(lS)),
     "loc·roll·ρ8·μ0 ln S비 점추정 %+.5f · 비 %.4f" % (lS.mean(), math.exp(lS.mean())),
     "흐름 %d 개 × 재추출 %d: 95%% 상한 평균 %+.5f · SD %.5f · 최소 %+.5f · 최대 %+.5f" % (NS, NB, his.mean(), his.std(ddof=1), his.min(), his.max()),
     "ln 0.8 = %+.5f · 상한 < ln 0.8(❌ 쪽) 비율 %.3f (%d/%d)" % (a08, (his < a08).mean(), int((his < a08).sum()), NS),
     "판정 아님 — 판정은 판정기 흐름(_RESULT_mini27.json P4.lnS_loc)에 잠긴 🟡"]
open(os.path.join(HERE, "_RECOUNT27_p4b_streams.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
