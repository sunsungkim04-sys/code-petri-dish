#!/usr/bin/env python3
"""p4win_posthoc_weights.py — P4 같은 창 재계수의 무게 감도 (사후 · 1차 결과 _RESULT_p4win.txt 를 본 뒤 썼다 · 판정 아님)

왜: 1차 정의(창마다 같은 무게)에서 Table S11 P4b 가 🟡(상한 − ln 0.8 = +0.0011) → ❌(−0.0173)로 바뀌었다.
그 행은 띠 경계에 있어 '같은 시각' 을 어떤 공통 무게로 합치느냐에 따라 달라질 수 있다. 공통 무게 둘을 더 본다.
  같은 창     창마다 같은 무게(1차 · p4win_analyze.py 와 같은 값이 나와야 한다 — 대조)
  L2 시간     두 계통의 창 R 을 L2 positions 몫으로 합친다 — L2 의 R 은 L2 합계 그대로, L4 의 R 을 L2 시각으로 옮긴다
  L4 시간     L4 positions 몫으로 합친다 — L4 의 R 은 L4 합계 그대로, L2 의 R 을 L4 시각으로 옮긴다
  P4b 는 무게 준 R 로 (R2 − 1)·2 / (R4 − 1)·4 를 만든다. 재추출 · 자격 배경 · 띠는 p4win_analyze.py 와 같다.
쓰는 법: python3 p4win_posthoc_weights.py p4w/h26 p4w/h27 > _RESULT_p4win_posthoc.txt
"""
import math
import os
import sys

import numpy as np

import p4win_analyze as A


def weighted(arm, W, mode):
    d2, p2 = W["1_2"][:, A.I_DRAW].astype(float), W["1_2"][:, A.I_POS].astype(float)
    d4, p4 = W["0_4"][:, A.I_DRAW].astype(float), W["0_4"][:, A.I_POS].astype(float)
    ok = (p2 >= A.MIN_POS_WIN) & (p4 >= A.MIN_POS_WIN)
    r2, r4 = d2[ok] / p2[ok], d4[ok] / p4[ok]
    if mode == "eq":
        lnd = float(np.mean(np.log(r2) - np.log(r4)))
        lns = float(np.mean(np.log((r2 - 1) * 2) - np.log((r4 - 1) * 4)))
        return lnd, lns
    v = p2[ok] if mode == "L2" else p4[ok]
    v = v / v.sum()
    R2, R4 = float((v * r2).sum()), float((v * r4).sum())
    return math.log(R2) - math.log(R4), math.log((R2 - 1) * 2) - math.log((R4 - 1) * 4)


def main():
    run26, run27 = sys.argv[1], sys.argv[2]
    print("== p4win_posthoc_weights — 무게 감도 (사후 · 1차 결과를 본 뒤 · 판정 아님)")
    print("   이 파일 sha256 %s · p4win_analyze %s" % (A.sha_file(os.path.abspath(__file__)), A.VERSION))
    bad = []
    data = {}
    for study, run in (("26", run26), ("27", run27)):
        A.check_manifest(run, bad)
        data[study], _ = A.check_run(run, os.path.join(A.HERE, A.STUDY[study]["ref"]), study,
                                     A.STUDY[study]["seeds"], bad, strict=True)
    if bad:
        print("** 가드 실패 %d — 내지 않는다" % len(bad))
        for b in bad[:20]:
            print("   ", b)
        return 2
    cells = (("26", "glob", "haebu26", "P4"), ("27", "loc", "haebu27", "P4|loc"), ("27", "glob", "haebu27", "P4|glob"))
    keep = {}
    for study, geom, pre, cell in cells:
        rows = []
        for s in A.STUDY[study]["seeds"]:
            arm, W = data[study][s][geom]
            if A.bg_stats(arm, W)["qual"]:
                rows.append((arm, W))
        idx = A.boot_idx(pre, cell, len(rows))
        print("\n-- h%s %s · 자격 배경 %d" % (study, geom, len(rows)))
        for mode, lab in (("eq", "같은 창(1차)"), ("L2", "L2 시간"), ("L4", "L4 시간")):
            v = np.array([weighted(a, W, mode) for a, W in rows])
            lnd, lns = v[:, 0], v[:, 1]
            cd, cs = A.ci(lnd[idx].mean(axis=1)), A.ci(lns[idx].mean(axis=1))
            keep[(study, geom, mode)] = (lnd, idx)
            print("   %-10s ln R2 − ln R4 %+.4f [%+.4f, %+.4f] · (R−1)L 비 %.3f [%.3f, %.3f] · 꼴 %s · 상한 − ln 0.8 %+.4f"
                  % (lab, lnd.mean(), cd[0], cd[1], math.exp(lns.mean()), math.exp(cs[0]), math.exp(cs[1]),
                     A.band(*cs), cs[1] - A.S_BAND[0]))
    print("\n-- h27 P4a [loc] − [glob] (독립 재추출 · 동결 판정기와 같은 행)")
    for mode, lab in (("eq", "같은 창(1차)"), ("L2", "L2 시간"), ("L4", "L4 시간")):
        (vl, ia), (vg, ib) = keep[("27", "loc", mode)], keep[("27", "glob", mode)]
        lo, hi = A.ci(vl[ia].mean(axis=1) - vg[ib].mean(axis=1))
        print("   %-10s %+.4f [%+.4f, %+.4f] · 꼴(하한 > 0) %s" % (lab, vl.mean() - vg.mean(), lo, hi, "✅" if lo > 0 else "❌"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
