#!/usr/bin/env python3
"""mini31_analyze.py — 해부 31 판정기: 최소 모형의 역전 실패가 종류인가 범위인가

사전등록: 볼트 `사전등록/해부31-사전등록-2026-10-06.md` (§2 예측 · §3 사다리 · §4 가드)

🔑 **동결된 `mini27_analyze.py` 를 고치지 않는다 — import 해서 쓴다.**
   Δ 추정(`delta_table`) · 재표집 색인(`boot_idx`) · 구간(`ci`) · 교차(`cens_flip`) ·
   구조 가드(`geom_guards`) 를 그 모듈의 것을 **그대로** 부른다. 정의를 다시 쓰지 않으므로
   해부 27 과 Δ 가 한 글자도 달라지지 않는다. 새로 쓴 것은 **μ 사다리와 §2 판정뿐**이다.

왜 새 스크립트인가: `mini27_analyze.py` 는 μ 사다리(`MUS_ALL`)와 판정 시드(`JUDGE_SEEDS`)를
코드에 박고 산출물이 그와 같은지 단언한다. 해부 31 은 사다리가 0…4% 이고 시드가 31001 이라
그 단언에 걸린다. 동결 파일은 건드리지 않는 것이 이 프로젝트 규칙이다.

실행
  python3 mini31_analyze.py dry   DIR            판정 정의 마른 실행 — **Δ 값을 인쇄하지 않는다**
  python3 mini31_analyze.py judge DIR8 DIR32     본판정 (ρ8: loc+glob · ρ32: loc)
"""
import json
import math
import os
import sys
import warnings

import numpy as np

import mini27_analyze as M27        # 동결 모듈 — 읽기만 한다

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ── 해부 31 이 새로 정하는 것 (사전등록 §3) ────────────────────────────────
SCALE31 = 4.0                                   # 원 사다리 × 4.0 → 0…4%
MUS31 = [round(m * SCALE31, 9) for m in M27.MUS_ORIG]
JUDGE_SEEDS31 = list(range(31001, 31031))
PILOT_SEEDS31 = list(range(31901, 31911))
N_BG31 = M27.N_BG                               # 20 — 해부 27 과 같게
GEOMS31 = {8.0: ["loc", "glob"], 32.0: ["loc"]}  # place 는 안 돌린다
P2_BAND = (0.013, 0.035)                        # 사전등록 §2 P2 — 2.14% 의 ±60%
P2_PRED = 0.0214                                # 내가 계산한 예측(원고 값 아님 · §5 고백)


def band_label(x):
    if x is None or not math.isfinite(x):
        return "교차 없음"
    return "%.4f%% %s" % (100 * x, "띠 안 ✅" if P2_BAND[0] <= x <= P2_BAND[1] else "띠 밖 ❌")


def gate31(dirs):
    """해부 31 용 관문. 구조 가드는 M27.geom_guards 를 strict=False 로 부른다
    (팔 수 단언은 M27 의 15점 사다리에 묶여 있어 여기선 따로 센다)."""
    bad, BG, data = [], {}, {}
    for rho, d in dirs.items():
        F = M27.load(d, rho)
        if not F:
            bad.append("ρ%g 산출물 없음: %s" % (rho, d))
            continue
        for s in F:
            if s in PILOT_SEEDS31:
                bad.append("판정 디렉터리에 파일럿 시드 %d" % s)
            if s not in JUDGE_SEEDS31:
                bad.append("후보 밖 시드 %d (해부 31 은 31001~31030)" % s)
        need = GEOMS31[rho]
        qual = []
        for s in sorted(F):
            z = F[s]
            if all(g in z["by_geom"] and z["by_geom"][g]["grow"]["qualified"] for g in need):
                qual.append(s)
            if len(qual) == N_BG31:
                break
        if len(qual) < N_BG31:
            bad.append("ρ%g 자격 배경 %d < %d — 시드를 더 돌릴 것(31030 까지)" % (rho, len(qual), N_BG31))
        BG[rho], data[rho] = qual, F

        for s in qual:
            z = F[s]
            if z["version"] != M27.VERSION:
                bad.append("판 불일치 ρ%g s%d %s (기대 %s)" % (rho, s, z["version"], M27.VERSION))
            for k, v in M27.CONFIG.items():
                if abs(float(z["config"][k]) - float(v)) > 1e-12:
                    bad.append("설정 %s ρ%g s%d %s≠%s" % (k, rho, s, z["config"][k], v))
            # 사다리가 해부 31 것인가 (이 스크립트의 존재 이유)
            got = [round(float(m), 9) for m in z["config"]["mus"]] if "mus" in z["config"] else None
            if got is not None and got != MUS31:
                bad.append("ρ%g s%d μ 사다리가 해부 31 것이 아니다: %s" % (rho, s, got[:4]))
            for g in need:
                M27.geom_guards(z, rho, s, g, bad, strict=False)
                arms = z["by_geom"][g]["arms"]
                exp = len(M27.RULES) * len(MUS31) * (1 + M27.NREF)
                if len(arms) != exp:
                    bad.append("팔 수 ρ%g s%d %s: %d ≠ %d" % (rho, s, g, len(arms), exp))
    return bad, BG, data


def deltas(data, BG, rho, geom, rule):
    """Δ 표 — M27 의 추정기를 그대로 부른다."""
    return M27.delta_table(data, BG, rho, geom, rule, MUS31)


def judge(d8, d32, json_out=None):
    bad, BG, data = gate31({8.0: d8, 32.0: d32})
    print("== 해부 31 관문")
    for b in bad:
        print("   🔴", b)
    if bad:
        print("** 관문 실패 — 판정하지 않는다(사람 결정). 사전등록 §6 실패 조건 4")
        return 3
    print("   ✅ 통과 — 배경 ρ8 %d · ρ32 %d · 사다리 0…%.1f%%"
          % (len(BG[8.0]), len(BG[32.0]), 100 * MUS31[-1]))

    res, curves = {}, {}
    for rho in (8.0, 32.0):
        for geom in GEOMS31[rho]:
            for rule in M27.RULES:
                D, drop, extl = deltas(data, BG, rho, geom, rule)
                n = len(BG[rho])
                row, pts = [], []
                print("\n-- Δ · ρ%g · %s · %s" % (rho, geom, rule))
                for mu in MUS31:
                    idx = M27.boot_idx("D|%g|%s|%s|%g" % (rho, geom, rule, mu), n)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        bm = np.nanmean(D[mu][idx], axis=1)
                    lo, hi = M27.ci(bm)
                    pt = float(np.nanmean(D[mu]))
                    print("   %6.3f%% | %+.3f [%+.3f, %+.3f] | 뺀 %d · 소멸 %d"
                          % (100 * mu, pt, lo, hi, drop[mu], extl[mu]))
                    res["D|%g|%s|%s|%g" % (rho, geom, rule, mu)] = [pt, lo, hi, drop[mu], extl[mu]]
                    row.append((mu, pt, lo, hi, drop[mu]))
                    pts.append(pt)
                    if drop[mu] > M27.MAX_DROP:
                        print("   ** 이 칸 판정 불가(뺀 배경 > %d)" % M27.MAX_DROP)
                curves[(rho, geom, rule)] = row
                res["curve|%g|%s|%s" % (rho, geom, rule)] = pts

    # ── §2 판정 ───────────────────────────────────────────────────────────
    def first_flip(rho, geom, rule):
        """P1 기준 — 상한 < 0 인 첫 μ. 보간하지 않는다(사전등록 §3)."""
        for mu, pt, lo, hi, dr in curves[(rho, geom, rule)]:
            if dr <= M27.MAX_DROP and hi < 0:
                return mu
        return None

    print("\n== 사전등록 §2 판정")
    f_loc = first_flip(8.0, "loc", "roll")
    p1 = f_loc is not None
    print("P1  loc·draw-first·ρ8 가 사다리 안에서 부호를 바꾼다      : %s  (%s)"
          % ("✅" if p1 else "❌", band_label(f_loc)))

    if p1:
        p2 = P2_BAND[0] <= f_loc <= P2_BAND[1]
        print("P2  교차가 %.1f–%.1f%% 안 (예측 %.2f%%)                   : %s"
              % (100 * P2_BAND[0], 100 * P2_BAND[1], 100 * P2_PRED, "✅" if p2 else "❌"))
    else:
        p2 = None
        print("P2  (P1 미통과로 판정 안 함 — 사전등록 §6 실패 조건 1)")

    f_find = first_flip(8.0, "loc", "find")
    p3 = f_find is None
    print("P3  find-first 는 같은 범위에서 교차하지 않는다            : %s  (%s)"
          % ("✅" if p3 else "🔴", band_label(f_find)))

    f_rich = first_flip(32.0, "loc", "roll")
    p4 = f_rich is None
    print("P4  밀도 32(넉넉)는 교차하지 않는다                        : %s  (%s)"
          % ("✅" if p4 else "🔴", band_label(f_rich)))

    f_glob = first_flip(8.0, "glob", "roll")
    if f_glob is None:
        p5, how = True, "교차 없음"
    elif f_loc is None:
        p5, how = None, "loc 이 안 교차해 비교 불가"
    else:
        p5, how = f_glob > f_loc, "glob %.4f%% 대 loc %.4f%%" % (100 * f_glob, 100 * f_loc)
    print("P5  glob 은 교차하지 않거나 loc 보다 늦다                  : %s  (%s)"
          % ({True: "✅", False: "❌", None: "—"}[p5], how))

    if not p3 or not p4:
        print("\n🔴🔴 사전등록 §6 실패 조건 3 — 대조에서도 뒤집혔다. 기제 설명이 무너진다.")
        print("     원고 §4.4 논증을 다시 봐야 한다. 이 결과만으로 Limits 를 고치지 말 것.")

    res["P"] = {"P1": p1, "P2": p2, "P3": p3, "P4": p4, "P5": p5,
                "flip": {"loc_roll_r8": f_loc, "loc_find_r8": f_find,
                         "loc_roll_r32": f_rich, "glob_roll_r8": f_glob},
                "band": list(P2_BAND), "pred": P2_PRED, "mus": MUS31}
    res["BG"] = {str(k): v for k, v in BG.items()}
    if json_out:
        json.dump(res, open(json_out, "w"), ensure_ascii=False)
        print("\n기록:", json_out)
    return 0 if (p1 is not None and p3 and p4) else 1


def dry(d):
    """판정 정의 마른 실행 — Δ 값을 **인쇄하지 않는다**(사전등록 §4).
    표가 만들어지고 P1 판정 코드가 도는지만 본다. 📕 feedback_prereg_definition_dryrun"""
    print("== 마른 실행 — Δ 값은 인쇄하지 않는다")
    bad, BG, data = gate31({8.0: d})
    print("관문 문제 %d 건%s" % (len(bad), "" if not bad else " (마른 실행에서는 치명적이지 않다)"))
    for b in bad[:8]:
        print("   ·", b)
    if not BG.get(8.0):
        print("🔴 배경 0 — 산출물을 먼저 만들 것")
        return 1
    D, drop, extl = deltas(data, BG, 8.0, "loc", "roll")
    shape = {("칸 수"): len(D), ("배경 수"): len(BG[8.0])}
    print("Δ 표 모양:", shape, "· 사다리", ["%.3f%%" % (100 * m) for m in MUS31])
    finite = all(np.isfinite(np.asarray(D[m], dtype=float)).any() for m in MUS31)
    print("각 칸에 유한한 값이 있나:", "✅" if finite else "🔴")
    n = len(BG[8.0])
    mu0 = MUS31[0]
    idx = M27.boot_idx("D|8|loc|roll|%g" % mu0, n)
    print("재표집 색인 모양:", np.shape(idx), "· 기대 (2000, %d)" % n)
    print("P1 판정 코드 경로: 상한 < 0 인 첫 μ 를 찾는다 — 보간 없음 ✅")
    print("P2 띠: %.1f–%.1f%% · 예측 %.2f%%" % (100 * P2_BAND[0], 100 * P2_BAND[1], 100 * P2_PRED))
    print("\n✅ 정의가 입력을 받아 돈다. Δ 값은 보지 않았다.")
    return 0


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "dry":
        return dry(sys.argv[2])
    if len(sys.argv) >= 4 and sys.argv[1] == "judge":
        out = None
        if "--json" in sys.argv:
            out = sys.argv[sys.argv.index("--json") + 1]
        return judge(sys.argv[2], sys.argv[3], out)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
