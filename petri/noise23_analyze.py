#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 23 판정 — 위치 대조 N3 실패는 위치 효과인가, 계측 잡음인가 (새 실행 없음 · 재분석)
사전등록: 해부23-사전등록-2026-10-01.md §4 · §6. 판정선은 본 자료를 보기 전에 고정했다.

  자료(목록 모드 · --wt 밑 코드 · μ 0 · 원래 규칙 — 세 배경 세트 모두 같은 명령):
     S1 첫 세트   inv15B_<밑>/mu0   배경 3~22(20)          밑 코드 여섯     (N3 원본 = 해부 12·13 이웃 기록 — 해부 15 G0 로 목록 = 이웃 확인됨)
     S2 해부15 E  inv15E_<밑>/mu0   배경 26~50 의 20        밑 코드 여섯
     S3 해부20 P  inv20P_<밑>/mu0   배경 76~100 의 20       밑 코드 여덟(acld 는 해부 20 선별로 빠짐)
  한 파일 = 한 (밑 코드, 배경). 팔 = ref 1(밑 코드 · 뒤 난수 A) + neutral 10(밑 코드 · 뒤 난수 B1..B10) + mut(n 삽입 전부 · 뒤 난수 A).
  끼우는 개체는 배경 시드로만 정해져 한 파일의 모든 팔이 같은 자리에 끼운다(dms.js inject).

  s   = ln((m_T + 0.5)/m_0) − ln((n_T − m_T + 0.5)/(n_0 − m_0))                        (해부 13 · 15 · 20 과 같다)
  N3  같은 고리 n 쌍 (a < c 정렬) 의 ΔΔ = 평균_배경 [s(c) − s(a)] · 실패 = |ΔΔ| > 0.15   (해부 13 N3 · 문턱은 바꾸지 않는다)
  잡음 바닥  같은 파일 묶음의 WT 팔 11 개(ref + neutral 10 — 같은 코드 · 같은 끼운 자리 · 뒤 난수만 다름)
             두 팔 (i < j) 에 N3 와 같은 정의를 적용한 값 m_ij = 평균_배경 [s(j) − s(i)] 의 분포
             = '똑같은 코드 둘이 뒤 난수만 달라 얼마나 달라 보이나' · 같은 집계 단위(배경 20 평균)
  q95 (주)   그 분포의 |값| 95% 분위수를 **합친 분산**으로 추정: σ̄² = 평균_배경(배경 안 11 팔 s 의 표본분산, ddof 1)
             SD0 = √(2 σ̄² / n_배경) · q95 = 1.959964 × SD0   (같은 코드 두 팔의 배경 평균 차 ~ N(0, SD0²))
  q95_hi     배경 재추출 2,000 번(σ̄² 를 다시 구함) → q95 의 97.5% 분위수. 난수 = (20261023, 칸 이름) 키
  q95_emp    서술: 55 쌍 |m_ij| 의 경험 95% 분위수(pct) · 경험 초과율(55 쌍 중 |m_ij| > q95 의 몫 — 기대 ≈ 0.05)
             (주 판정에 쓰지 않는 이유 = 사전등록 §4: 팔 11 개에서 나온 서로 얽힌 55 값의 분위수는 아래로 치우치고 흔들린다 — 합성 점검)
  분류(같은 고리 쌍마다)  N 잡음 안: |ΔΔ| ≤ q95 · A 모호: q95 < |ΔΔ| ≤ q95_hi · P 잡음 위: |ΔΔ| > q95_hi
  C0 계측 신뢰  (세트, 밑 코드)마다 해부 20 선별 정의(한 팔이라도 끝/처음 < 0.10 인 배경 수 K) — K ≥ 1 이면 그 (세트, 밑 코드)는 판정 안 함
  P_rob(검토 M3 · 결과 전)  |ΔΔ| > q95_hi · max(1, √R_var) — 쌍 자신의 흔들림(1.96·SD(d)/√n 꼴)도 넘는 P
  종합 V  F = 판정하는 (세트, 밑 코드)의 N3 실패 쌍 · ALL3 = 세 세트 모두에서 판정되는 같은 고리 배열 (밑, a, c) 전부(선택 없음)
     V+ '위치 효과가 잡음 위'   : ALL3 의 어떤 배열이 **세 세트 모두에서 P_rob · ΔΔ 부호 같음**(검토 M2(a) · M3 — 결과 전 변경)
     V? (바닥 위 · 자기 흔들림 위 아님) : V+ 아님 · 세 세트 모두 P 같은 부호인 배열은 있으나 P_rob 로는 안 섬
     V? '실패 쌍 없음'           : 위 둘 아님 · F 가 비었음(검토 S6)
     V0 '잡음으로 설명됨'       : F 에 P 가 하나도 없고 N 이 F 의 절반 이상
     V? '모호'                  : 그 밖
     서술 'V+ 약'               : 정확히 두 세트에서 P · 같은 부호인 배열(그 배열을 F 에 넣은 세트가 어디인지 함께 — 판정 아님)
  서술  N3 '거짓 실패율'(같은 코드 두 팔이 0.15 를 넘을 확률 — 정규 근사 · 55 쌍 경험) · 그 율로 본 기대 실패 수 대 관측
        CRN 비(검토 S1 · mut 코드마다 Var_배경(s_mut − s_ref) / 평균_k Var_배경(s_mut − s_neutral_k) 의 중앙 — 1 미만 = 뒤 난수 A 결합)
        R_var(쌍마다: 배경 간 차 분산 / 2σ̄² — 쌍의 흔들림이 같은 코드 잡음의 몇 배인가)
  G0-a 에 배경 목록 세 겹 일치(규칙 = 하드코딩 §2 표 = 파일 · 각 20) · 읽은 모든 파일 = dms23_data_manifest.sha256(검토 S2 · S3)
실행: python3 noise23_analyze.py [petri 폴더]            → _RESULT_noise23.json (판정)
      PETRI_PILOT=1 python3 noise23_analyze.py [폴더]   → _RESULT_pilot23r.json (파일럿 · 판정 아님 · 검토 개정 뒤)
      python3 noise23_analyze.py --selftest              → 합성 자료로 분류 함수 점검(파일 안 읽음)
"""
import glob
import hashlib
import json
import math
import os
import random
import sys

# ---- 판정선 (사전등록 §6 · 결과 전 고정) ----
TOL_NEG = 0.15          # N3 원래 문턱 — 이 실험은 바꾸지 않는다(사후 문턱 조정 금지)
Q = 0.95                # 잡음 바닥 분위수
N_BOOT = 2000
RNG_SEED = 20261023     # 재추출 난수는 (시드, 칸 이름)으로 키
COLLAPSE = 0.10         # 해부 20 선별 정의 그대로
MIN_SETS_REPL = 3       # V+ — P_rob 가 같은 부호로 나와야 하는 세트 수(검토 M2(a): 세 세트 모두 · 선택 없는 ALL3 배열)
CRN_READ = 0.8          # 서술 읽기: CRN 비 중앙 < 0.8 이고 1 미만인 코드 몫 ≥ 0.75 이면 'A 결합 남음' 단서

# 검토 S2: §2 표의 배경 목록을 하드코딩 — 규칙에서 다시 뽑은 목록 · 실제 파일과 셋이 같아야 한다(각 20)
SEEDS_S1 = list(range(3, 23))
SEEDS_S2 = [26, 27, 28, 30, 31, 32, 33, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 50]
SEEDS_S3 = [76, 77, 79, 80, 81, 83, 84, 85, 86, 87, 88, 89, 92, 93, 94, 95, 96, 98, 99, 100]
MANIFEST = "dms23_data_manifest.sha256"

OLD6 = ["racld", "rascld", "rsacld", "racldx", "acld", "rascled"]
P20 = ["rcald", "crald", "reasccld", "racld", "rascld", "rsacld", "racldx", "rascled"]
FAM13 = {"acld", "rascled"}   # _RESULT_pairs13.json · 나머지 넷은 _RESULT_pairs13_racldfam.json

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
SELFTEST = "--selftest" in sys.argv
BASE = os.path.expanduser(ARGS[0] if ARGS else "~/petri")
PILOT = os.environ.get("PETRI_PILOT") == "1"
bad, notes = [], []
READ = {}   # 읽은 파일(BASE 기준 상대 경로) -> sha256 — 검토 S3 자료 목록 대조


def jload(path):
    """JSON 을 읽으며 sha256 을 기록한다(읽은 모든 입력이 자료 목록에 있어야 한다)"""
    with open(path, "rb") as fh:
        raw = fh.read()
    READ[os.path.relpath(path, BASE)] = hashlib.sha256(raw).hexdigest()
    return json.loads(raw.decode("utf-8"))


# ---------------- 기존 판정기와 글자 그대로 같은 함수 ----------------
def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def s_of(arm):
    r = arm["series"]; n0, m0 = r[0][1], r[0][2]; nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def same_loop_pairs(codes):
    cl = sorted(codes); out = []
    for i, a in enumerate(cl):
        for b in cl[i + 1:]:
            if codes[a] == codes[b]: out.append((a, b))
    return out


# ---------------- 해부 23 의 함수 (판정 정의 — §4) ----------------
def dd_mean(D, seeds, a, c):
    """N3 정의: 평균_배경 [s(c) − s(a)] — D[seed] = {이름: s}"""
    return sum(D[s][c] - D[s][a] for s in seeds) / len(seeds)


def floor_values(W, seeds):
    """잡음 바닥: WT 팔 11 개의 모든 쌍 (i < j) 에 N3 와 같은 정의 — W[seed] = [s_ref, s_n1, …, s_n10]"""
    k = len(W[seeds[0]]); out = []
    for i in range(k):
        for j in range(i + 1, k):
            out.append(sum(W[s][j] - W[s][i] for s in seeds) / len(seeds))
    return out


Z975 = 1.959964        # |N(0,1)| 의 95% 분위수


def pooled_var(W, seeds):
    """배경마다 11 팔 s 의 표본분산(ddof 1) → 배경 평균"""
    v = []
    for s in seeds:
        x = W[s]; m = sum(x) / len(x)
        v.append(sum((a - m) ** 2 for a in x) / (len(x) - 1))
    return sum(v) / len(v)


def floor_q(W, seeds):
    """주 잡음 바닥 q95 — 같은 코드 두 팔의 배경 평균 차 |·| 의 95% 분위수(합친 분산 · 정규 근사)"""
    return Z975 * math.sqrt(2 * pooled_var(W, seeds) / len(seeds))


def floor_q_emp(W, seeds):
    """서술 — 55 쌍의 경험 분위수"""
    return pct([abs(v) for v in floor_values(W, seeds)], Q)


def false_fail_rate(W, seeds):
    """서술 — 같은 코드 두 팔이 N3(0.15)를 '실패' 할 확률: 정규 근사 · 경험(55 쌍)"""
    sd0 = math.sqrt(2 * pooled_var(W, seeds) / len(seeds))
    par = math.erfc(TOL_NEG / (sd0 * math.sqrt(2))) if sd0 > 0 else 0.0
    fv = floor_values(W, seeds)
    return par, sum(1 for v in fv if abs(v) > TOL_NEG) / len(fv)


def floor_q_boot(W, seeds, name):
    rng = random.Random("%d|%s" % (RNG_SEED, name))
    qs = []
    for _ in range(N_BOOT):
        rs = [seeds[rng.randrange(len(seeds))] for _ in seeds]
        qs.append(floor_q(W, rs))
    return pct(qs, 0.025), pct(qs, 0.975)


def var_ratio(D, seeds, a, c, pv):
    """서술 R_var: 쌍의 배경 간 차 분산(ddof 1) / 같은 코드 두 팔 차의 분산(2 σ̄²) — 1 이면 쌍의 흔들림이 뒤 난수 잡음만큼"""
    d = [D[s][c] - D[s][a] for s in seeds]; m = sum(d) / len(d)
    v = sum((x - m) ** 2 for x in d) / (len(d) - 1) if len(d) > 1 else float("nan")
    return v / (2 * pv) if pv > 0 else float("nan")


def classify(x, q, qhi):
    x = abs(x)
    return "N" if x <= q else ("A" if x <= qhi else "P")


def classify_rob(x, qhi, rv):
    """검토 M3(결과 전): P 를 쌍 자신의 흔들림으로 다시 본다 — |ΔΔ| > q95_hi · max(1, √R_var) 이면 P_rob.
    q95_hi · √R_var ≈ 1.96 · SD_배경(d) / √n · (q95_hi / q95) = 쌍 고유 95% 상한. R_var 가 NaN 이면 P_rob 아님."""
    if rv != rv:
        return False
    return abs(x) > qhi * max(1.0, math.sqrt(rv))


def _var(xs):
    m = sum(xs) / len(xs)
    return sum((a - m) ** 2 for a in xs) / (len(xs) - 1) if len(xs) > 1 else float("nan")


def crn_ratio(W, M, seeds):
    """서술(검토 S1 으로 D-CRN 을 바꿈): mut 코드마다 Var_배경(s_mut − s_ref[A]) / 평균_k Var_배경(s_mut − s_neutral_k[B_k]).
    배경 효과는 차에서 빠진다. 뒤 난수 A 가 mut 과 ref 를 묶으면 분자가 작아져 1 미만. 반환 = (중앙, 1 미만인 코드 몫, 코드 수)"""
    vals = []
    for m in M[seeds[0]]:
        num = _var([M[s][m] - W[s][0] for s in seeds])
        den = [_var([M[s][m] - W[s][k] for s in seeds]) for k in range(1, len(W[seeds[0]]))]
        den = [d for d in den if d == d]
        if num == num and den and sum(den) > 0:
            vals.append(num / (sum(den) / len(den)))
    if not vals:
        return float("nan"), float("nan"), 0
    return pct(vals, 0.5), sum(1 for v in vals if v < 1) / len(vals), len(vals)


def verdict(F, CLS):
    """F = [(set, b, a, c, dd)] 실패 쌍 · CLS[(set, b, a, c)] = (분류, dd, P_rob) — 판정 칸의 모든 같은 고리 쌍.
    검토 M2(a)(결과 전 변경): V+ 는 F 와 무관하게(선택 없이) 세 세트 모두에서 판정되는 배열(ALL3)만 본다 —
    세 세트 모두 P_rob · 같은 부호. 반환 = (종합, 재현 목록, 서술 dict)"""
    sets_all = sorted({st for st, _, _, _ in CLS})
    arrs = {}
    for (st, b, a, c), v in CLS.items():
        arrs.setdefault((b, a, c), {})[st] = v
    all3 = sorted(arr for arr, d in arrs.items() if len(d) >= MIN_SETS_REPL and len(sets_all) >= MIN_SETS_REPL and set(d) == set(sets_all))
    repl_rob, repl_p, weak = [], [], []
    for arr in sorted(arrs):
        d = arrs[arr]
        for sign in (1, -1):
            sp = sorted(st for st, (k, dd, rob) in d.items() if k == "P" and dd * sign > 0)
            sr = sorted(st for st, (k, dd, rob) in d.items() if k == "P" and rob and dd * sign > 0)
            if arr in all3 and len(sr) >= MIN_SETS_REPL:
                repl_rob.append((arr, sign, sr))
            elif arr in all3 and len(sp) >= MIN_SETS_REPL:
                repl_p.append((arr, sign, sp))
            elif len(sp) == 2:
                sel = sorted({st for st, b, a, c, _ in F if (b, a, c) == arr})
                weak.append((arr, sign, sp, sel))
    desc = {"n_all3": len(all3), "repl_p_not_rob": [[list(a), s, ss] for a, s, ss in repl_p],
            "weak_2sets": [[list(a), s, ss, sel] for a, s, ss, sel in weak]}
    if repl_rob:
        return "V+ 위치 효과가 잡음 위", repl_rob, desc
    if repl_p:
        return "V? 바닥 위지만 쌍 자신의 흔들림 위는 아님", [], desc
    if not F:
        return "V? 모호 (판정 칸에 N3 실패 쌍 없음)", [], desc
    kinds = [CLS[(st, b, a, c)][0] for st, b, a, c, _ in F]
    if "P" not in kinds and kinds.count("N") * 2 >= len(kinds):
        return "V0 잡음으로 설명됨", [], desc
    return "V? 모호", [], desc


# ---------------- 자기 점검(합성 자료) ----------------
def selftest():
    rng = random.Random(7); seeds = list(range(20)); ok = True
    W, D = {}, {}
    for s in seeds:
        bg = rng.gauss(0, 1.0)
        W[s] = [bg + rng.gauss(0, 0.30) for _ in range(11)]
        D[s] = {"a": bg + rng.gauss(0, 0.30), "c0": bg + rng.gauss(0, 0.30), "c1": bg + 1.0 + rng.gauss(0, 0.30)}
    q = floor_q(W, seeds); lo, hi = floor_q_boot(W, seeds, "selftest")
    k0 = classify(dd_mean(D, seeds, "a", "c0"), q, hi); k1 = classify(dd_mean(D, seeds, "a", "c1"), q, hi)
    # 기대: 같은 분포 잡음(SD 0.30 · 배경 20)의 차 평균 SD = 0.30·√2/√20 = 0.0949 → 참 q95 = 0.186 · 차 0 쌍은 N/A · 차 1.0 쌍은 P
    print("  selftest q95 %.3f [%.3f, %.3f](참 0.186) · 경험 q95 %.3f · 차 0 쌍 %s · 차 1.0 쌍 %s" % (q, lo, hi, floor_q_emp(W, seeds), k0, k1))
    ok &= k1 == "P" and k0 in ("N", "A") and lo <= q <= hi and lo <= 0.186 <= hi
    F = [("S1", "x", "a", "c1", 1.0)]
    S3c = lambda d1, d2, d3, r3=True: {("S1", "x", "a", "c1"): d1, ("S2", "x", "a", "c1"): d2, ("S3", "x", "a", "c1"): d3}
    # 검토 M2 요구 사례: 한 세트에서 실패(그 세트에서 P) + 다른 한 세트에서 P → V+ 이면 안 된다
    v1 = verdict(F, {("S1", "x", "a", "c1"): ("P", 1.0, True), ("S2", "x", "a", "c1"): ("P", 0.9, True), ("S3", "x", "a", "c1"): ("N", 0.01, False)})[0]
    v1b = verdict(F, {("S1", "x", "a", "c1"): ("P", 1.0, True), ("S2", "x", "a", "c1"): ("P", 0.9, True)})[0]   # S3 에 그 칸이 없음(ALL3 아님)
    v5 = verdict(F, S3c(("P", 1.0, True), ("P", 0.9, True), ("P", 0.8, True)))[0]          # 세 세트 P_rob 같은 부호 → V+
    v6 = verdict(F, S3c(("P", 1.0, True), ("P", 0.9, True), ("P", 0.8, False)))[0]         # 하나가 P_rob 아님(M3) → V? 바닥 위 …
    v2 = verdict(F, S3c(("P", 1.0, True), ("P", -0.9, True), ("P", 0.8, True)))[0]         # 부호 반대 → V?
    v3 = verdict(F, {("S1", "x", "a", "c1"): ("N", 1.0, False)})[0]
    v4 = verdict(F, {("S1", "x", "a", "c1"): ("A", 1.0, False)})[0]
    v7 = verdict([], {("S1", "x", "a", "c1"): ("N", 0.0, False)})[0]                       # 검토 S6: F 비었음 → V?
    print("  selftest 종합 규칙: 실패 1 세트 + 다른 1 세트 P → %s · 두 세트뿐(ALL3 아님) → %s · 세 세트 P_rob 같은 부호 → %s · 셋 P 중 하나 P_rob 아님 → %s · 부호 반대 → %s · N → %s · A → %s · F 빔 → %s"
          % (v1, v1b, v5, v6, v2, v3, v4, v7))
    ok &= (not v1.startswith("V+")) and (not v1b.startswith("V+")) and v5.startswith("V+") and v6.startswith("V? 바닥") \
        and v2.startswith("V?") and v3.startswith("V0") and v4.startswith("V?") and v7.startswith("V?")
    # 검토 M3: P_rob 문턱 = q95_hi · max(1, √R_var)
    ok &= classify_rob(0.20, 0.10, 1.0) and not classify_rob(0.20, 0.10, 9.0) and classify_rob(0.31, 0.10, 9.0) and not classify_rob(0.5, 0.1, float("nan"))
    # 검토 S1: CRN 비 — mut 이 ref 와 같은 잡음을 나누면(결합) 1 미만, 독립이면 ≈ 1
    Wc, Mc, Mi = {}, {}, {}
    for s in seeds:
        e = [rng.gauss(0, 0.30) for _ in range(11)]; bg = rng.gauss(0, 1.0)
        Wc[s] = [bg + x for x in e]
        Mc[s] = {"m": bg + e[0] + rng.gauss(0, 0.05)}      # 뒤 난수 A 결합: ref 잡음을 거의 그대로 나눔
        Mi[s] = {"m": bg + rng.gauss(0, 0.30)}             # 결합 없음
    rc, ri = crn_ratio(Wc, Mc, seeds)[0], crn_ratio(Wc, Mi, seeds)[0]
    print("  selftest CRN 비: 결합 %.3f(기대 ≪ 1) · 독립 %.3f(기대 ≈ 1) · P_rob 문턱 점검 %s" % (rc, ri, "통과" if classify_rob(0.31, 0.10, 9.0) else "실패"))
    ok &= rc < 0.2 and 0.5 < ri < 2.0
    # 결정성: 같은 이름 키 → 같은 구간
    ok &= floor_q_boot(W, seeds, "selftest") == (lo, hi)
    # 잡음 바닥 쌍 수 = 11C2
    ok &= len(floor_values(W, seeds)) == 55
    print("  selftest %s" % ("통과 ✅" if ok else "실패 ❌"))
    return ok


if SELFTEST:
    sys.exit(0 if selftest() else 1)


# ---------------- 세트 정의 ----------------
def racld_ok_seeds():
    ok = []
    for s in range(1, 101):
        z = jload(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
            ok.append(s)
    return ok


def pub_judge(sid):
    """판정 세트의 발표 N3 값 위치 — (파일, 키 꼴, 부호). 해부 13 판정기는 음성 쌍을 s(앞) − s(뒤) 로, 해부 15 · 20 은 s(뒤) − s(앞) 로 적었다"""
    return {"S1": lambda b: ("_RESULT_pairs13.json" if b in FAM13 else "_RESULT_pairs13_racldfam.json", "neg|%s|n|%s|%s", -1),
            "S2": lambda b: ("_RESULT_pairs15.json", "E2|%s|%s|%s", 1),
            "S3": lambda b: ("_RESULT_pairs20.json", "neg|%s|%s|%s", 1)}[sid]


def pub_count(pub, keyfmt, b):
    if keyfmt.startswith("E2"): return sum(1 for k in pub if k.startswith("E2|%s|" % b))
    if "|n|" in keyfmt: return sum(1 for k in pub if k.startswith("neg|%s|n|" % b))
    return sum(1 for k in pub if k.startswith("neg|%s|" % b))


if not PILOT:
    # 검토 대비(해부 20 R3 와 같은 뜻): 판정 모드에서는 자료 위치 · 배경 목록을 환경으로 바꿀 수 없다
    for v in ("PETRI_SETS", "PETRI_SEEDS", "PETRI_NBG", "PETRI_ROOT"):
        if os.environ.get(v):
            print("🚨 판정 모드에서 %s 를 쓸 수 없다 — 파일럿이면 PETRI_PILOT=1" % v); sys.exit(1)
    SETS = [
        {"id": "S1", "name": "첫 세트(배경 3~22)", "dir": "inv15B_{b}/mu0", "bases": OLD6, "seeds": SEEDS_S1,
         "pub": pub_judge("S1")},
        {"id": "S2", "name": "해부15 E(배경 26~50)", "dir": "inv15E_{b}/mu0", "bases": OLD6, "seeds": SEEDS_S2,
         "pub": pub_judge("S2")},
        {"id": "S3", "name": "해부20 P(배경 76~100)", "dir": "inv20P_{b}/mu0", "bases": P20, "seeds": SEEDS_S3,
         "pub": pub_judge("S3")},
    ]
    G0B = ("inv12nbr_racld/mu0", SEEDS_S1)
    SCREEN = "_RESULT_screen20.json"
    OUTFN = "_RESULT_noise23"
else:
    SETS = [
        {"id": "S1", "name": "파일럿 첫 세트 꼴(pilot15B · 배경 24 25)", "dir": "pilot15B_{b}/mu0", "bases": ["racld", "rascled"], "seeds": [24, 25], "pub": None},
        {"id": "S2", "name": "파일럿 E 꼴(pilot15E · 배경 51 52)", "dir": "pilot15E_{b}/mu0", "bases": ["racld", "rascled"], "seeds": [51, 52], "pub": None},
        {"id": "S3", "name": "파일럿 20 꼴(pilot20P · 배경 24 25)", "dir": "pilot20P_{b}/mu0", "bases": P20, "seeds": [24, 25],
         "pub": lambda b: ("_RESULT_pilot20.json", "neg|%s|%s|%s", 1)},
    ]
    G0B = ("stage12/pilot12nbr_racld/mu0", [24])
    SCREEN = "_RESULT_pilot20_screen.json"
    OUTFN = "_RESULT_pilot23r"

# 검토 S2: 배경 목록 세 겹 — 규칙(main/ 에서 '멸종 안 함 · 끝 우세 racld' — 해부 13·15·20 발사기 규칙) = 하드코딩 §2 표 · 각 20.
# main/ 은 배경 생성 기록일 뿐 판정 자료가 아니므로 파일럿에서도 같은 점검을 돈다. 파일과의 일치는 아래 G0-a(seen == seeds).
OK = racld_ok_seeds()
RULE = {"S1": OK[:20], "S2": [s for s in OK if s >= 26][:20], "S3": sorted([s for s in OK if s >= 53][::-1][:20])}
for sid, hard in (("S1", SEEDS_S1), ("S2", SEEDS_S2), ("S3", SEEDS_S3)):
    if RULE[sid] != hard or len(hard) != 20 or len(set(hard)) != 20:
        bad.append("배경 목록 %s: 규칙 %s ≠ 하드코딩 %s (또는 20 개 아님)" % (sid, RULE[sid], hard))
if not PILOT:
    for st in SETS:
        if st["seeds"] != {"S1": SEEDS_S1, "S2": SEEDS_S2, "S3": SEEDS_S3}[st["id"]]: bad.append("세트 %s 배경 ≠ 하드코딩" % st["id"])
    if G0B[1] != SEEDS_S1: bad.append("G0-b 배경 ≠ S1")

# 세트끼리 배경이 겹치면 '두 세트에서 재현' 이 같은 자료를 두 번 센 것이 된다(V+ 의 전제)
for i, s1 in enumerate(SETS):
    for s2 in SETS[i + 1:]:
        ov = sorted(set(s1["seeds"]) & set(s2["seeds"]))
        if ov:
            msg = "세트 배경 겹침 %s · %s: %s" % (s1["id"], s2["id"], ov)
            (notes if PILOT else bad).append(msg + (" — 파일럿이라 서술(재현 판정이 독립이 아님)" if PILOT else ""))

JUDGE_SEEDS = set()
if PILOT:
    # 판정 경로의 발표 값 키 꼴 점검 — 발표된 _RESULT_ 파일의 키 · 개수만 본다(본 원자료는 읽지 않는다)
    kc = {}
    for sid, bases in (("S1", OLD6), ("S2", OLD6), ("S3", P20)):
        for b in bases:
            fn, keyfmt, sign = pub_judge(sid)(b)
            pub = jload(os.path.join(BASE, fn))
            prs = same_loop_pairs(inserts(b, "n"))
            miss = [keyfmt % (b, a, c) for a, c in prs if keyfmt % (b, a, c) not in pub]
            npub = pub_count(pub, keyfmt, b)
            kc["%s|%s" % (sid, b)] = {"pairs": len(prs), "pub": npub, "missing": len(miss)}
            if miss or npub != len(prs): bad.append("판정 경로 키 꼴 %s %s: 없음 %d · 쌍 %d ≠ 발표 %d" % (sid, b, len(miss), len(prs), npub))
    print("  판정 경로 키 꼴 점검(발표 파일 키만): %s" % " · ".join("%s %d/%d" % (k, v["pairs"] - v["missing"], v["pub"]) for k, v in kc.items()))
    # 파일럿 배경이 본 판정 배경과 겹치지 않는지(본 자료의 값은 읽지 않는다 — 목록만)
    JUDGE_SEEDS = set(range(3, 23)) | set(range(26, 51)) | set(range(76, 101))
    for st in SETS:
        if set(st["seeds"]) & JUDGE_SEEDS: bad.append("파일럿 배경이 본 판정 배경과 겹친다 %s %s" % (st["id"], sorted(set(st["seeds"]) & JUDGE_SEEDS)))

OUT = {"pilot": PILOT, "tol_neg": TOL_NEG, "q": Q, "n_boot": N_BOOT, "rng_seed": RNG_SEED, "collapse": COLLAPSE,
       "sets": {st["id"]: {"name": st["name"], "dir": st["dir"], "bases": st["bases"], "seeds": st["seeds"]} for st in SETS}}

# 해부 20 선별 통과 목록 = S3 밑 코드
scr = jload(os.path.join(BASE, SCREEN))
if sorted(scr["pass"]) != sorted(P20): bad.append("S3 밑 코드 ≠ %s 의 pass %s" % (SCREEN, scr["pass"]))

# ---------------- 읽기 · 점검(G0) ----------------
DATA = {}   # (set, b) -> {"W": {seed: [11]}, "D": {seed: {code: s}}, "K": int, "minr": float}
for st in SETS:
    for b in st["bases"]:
        codes = list(inserts(b, "n"))
        fs = sorted(glob.glob(os.path.join(BASE, st["dir"].format(b=b), "dms_mat_%s_s*_list.json" % b)))
        W, D, seen, K, minr = {}, {}, [], 0, float("inf")
        gcs = {}
        for f in fs:
            z = jload(f); s = z["seed"]; seen.append(s)
            if abs(z["mu_assay"]) > 1e-12 or z["wt"] != b or z["arms_mode"] != "list" or z["cond"] != "mat": bad.append("설정 %s %s %d" % (st["id"], b, s))
            if z.get("find_first") or z.get("remember_die") or z.get("density") or z.get("grow_mu") is not None or z.get("age0") is not None or z.get("no_mat") or z.get("cosmic_off"):
                bad.append("규칙 · 세계 옵션 %s %s %d" % (st["id"], b, s))
            if z["sim_version"] != "0.3.5": bad.append("sim 판 %s %s %d" % (st["id"], b, s))
            if not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (st["id"], b, s))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (st["id"], b, s))
            if z["grow_extinct"] >= 0 or z["grow_top"][0][0] != "racld": bad.append("배경 %s %s %d" % (st["id"], b, s))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (st["id"], b, s))
            ref = [a for a in z["arms"] if a["kind"] == "ref"]
            neu = [a for a in z["arms"] if a["kind"] == "neutral"]
            muts = [a for a in z["arms"] if a["kind"] == "mut"]
            # WT 팔 형식: ref 1 · neutral 10 · 전부 밑 코드 · 이름 neutral:1..10 순서 · 뒤 난수 11 개 서로 다름 · mut 은 ref 와 같은 뒤 난수(A)
            if len(ref) != 1 or len(neu) != 10: bad.append("WT 팔 수 %s %s %d (ref %d · neutral %d)" % (st["id"], b, s, len(ref), len(neu))); continue
            if any(a["key"] != b for a in ref + neu): bad.append("WT 팔 코드 ≠ 밑 코드 %s %s %d" % (st["id"], b, s))
            if [a["labels"] for a in neu] != [["neutral:%d" % k] for k in range(1, 11)]: bad.append("neutral 이름 · 순서 %s %s %d" % (st["id"], b, s))
            if len({a["post_seed"] for a in ref + neu}) != 11: bad.append("WT 뒤 난수 중복 %s %s %d" % (st["id"], b, s))
            if any(a["post_seed"] != ref[0]["post_seed"] for a in muts): bad.append("mut 뒤 난수 ≠ A %s %s %d" % (st["id"], b, s))
            if len({json.dumps(a["inj"], sort_keys=True) for a in z["arms"]}) != 1: bad.append("끼우기가 팔마다 다르다 %s %s %d" % (st["id"], b, s))
            if [a["key"] for a in muts] != codes or len(z["arms"]) != 11 + len(codes): bad.append("팔 목록 %s %s %d" % (st["id"], b, s))
            W[s] = [s_of(a) for a in ref + neu]
            D[s] = {a["key"]: s_of(a) for a in muts}
            gcs[s] = z["grow_checksum"]
            col = False
            for a in z["arms"]:
                m0, mT = a["series"][0][2], a["series"][-1][2]
                if m0 <= 0: bad.append("처음 개체 0 %s %s %d %s" % (st["id"], b, s, a["key"])); continue
                r = mT / m0; minr = min(minr, r)
                if r < COLLAPSE: col = True
            K += col
        if sorted(seen) != sorted(st["seeds"]): bad.append("배경 목록 %s %s %s ≠ %s" % (st["id"], b, sorted(seen), sorted(st["seeds"])))
        DATA[(st["id"], b)] = {"W": W, "D": D, "K": K, "minr": minr, "gcs": gcs, "codes": inserts(b, "n")}
    # 같은 세트 · 같은 배경이면 밑 코드가 달라도 배경이 같아야 한다
    per_seed = {}
    for b in st["bases"]:
        for s, g in DATA[(st["id"], b)]["gcs"].items(): per_seed.setdefault(s, set()).add(g)
    nd = sum(1 for v in per_seed.values() if len(v) > 1)
    if nd: bad.append("밑 코드 사이 grow_checksum 불일치 %s %d 배경" % (st["id"], nd))

# G0-b: 같은 코드 · 같은 뒤 난수면 궤적이 똑같아야 한다(잡음 바닥의 원천이 뒤 난수뿐임을 확인)
g0b_chk = g0b_same = 0
for f in sorted(glob.glob(os.path.join(BASE, G0B[0], "dms_mat_racld_s*_nbr_racld.json"))):
    z = jload(f)
    if z["seed"] not in G0B[1]: continue
    r = [a for a in z["arms"] if a["kind"] == "ref"]; sf = [a for a in z["arms"] if a["kind"] == "mut" and a["labels"] == ["self"]]
    if len(r) == 1 and len(sf) == 1 and sf[0]["key"] == "racld" and sf[0]["post_seed"] == r[0]["post_seed"]:
        g0b_chk += 1; g0b_same += r[0]["series"] == sf[0]["series"]
if g0b_chk != len(G0B[1]) or g0b_same != g0b_chk: bad.append("G0-b 같은 코드 · 같은 난수 궤적 동일 %d/%d (기대 %d)" % (g0b_same, g0b_chk, len(G0B[1])))
OUT["G0b"] = {"dir": G0B[0], "checked": g0b_chk, "identical": g0b_same}

# G0-c: 다시 셈한 N3 값 = 이미 발표된 판정기 값(같은 정의인지 — 값 차 < 1e-9 · 쌍 수 같음)
g0c = {}
for st in SETS:
    if not st["pub"]:
        notes.append("%s: 대조할 발표 값 없음(파일럿)" % st["id"]); continue
    for b in st["bases"]:
        fn, keyfmt, sign = st["pub"](b)
        pub = jload(os.path.join(BASE, fn))
        dat = DATA[(st["id"], b)]; seeds = sorted(dat["D"])
        prs = same_loop_pairs(dat["codes"]); n = mx = 0
        npub = pub_count(pub, keyfmt, b)
        for a, c in prs:
            k = keyfmt % (b, a, c)
            if k not in pub or not seeds: bad.append("G0-c 발표 값 없음 %s %s" % (st["id"], k)); continue
            mine = dd_mean(dat["D"], seeds, a, c)
            d = abs(mine - sign * pub[k]["dd"]); mx = max(mx, d); n += 1
        if mx > 1e-9: bad.append("G0-c 다시 셈 ≠ 발표 %s %s (최대 차 %.2e)" % (st["id"], b, mx))
        if npub != len(prs): bad.append("G0-c 쌍 수 %s %s %d ≠ 발표 %d" % (st["id"], b, len(prs), npub))
        g0c["%s|%s" % (st["id"], b)] = {"pairs": n, "pub_pairs": npub, "max_abs_diff": mx, "file": fn}
OUT["G0c"] = g0c

# G0-a 자료 목록(검토 S3): 이 판정기가 읽은 모든 파일(main/ · 세 세트 원자료 · G0-b · 발표 _RESULT_) = 동결 때 적은 목록 · 해시 같음
READ_DIGEST = hashlib.sha256("".join("%s  %s\n" % (READ[k], k) for k in sorted(READ)).encode()).hexdigest()
OUT["read_files"] = len(READ); OUT["read_digest"] = READ_DIGEST
MAN_SHA = None
if not PILOT:
    mpath = os.path.join(BASE, MANIFEST)
    if not os.path.exists(mpath):
        bad.append("자료 목록 %s 없음" % MANIFEST)
    else:
        raw = open(mpath, "rb").read(); MAN_SHA = hashlib.sha256(raw).hexdigest()
        man = {}
        for line in raw.decode("utf-8").splitlines():
            if not line.strip(): continue
            h, p = line.split(None, 1); man[os.path.normpath(p.strip().lstrip("*"))] = h
        diff = [k for k in READ if man.get(os.path.normpath(k)) != READ[k]]
        extra = [k for k in man if k not in {os.path.normpath(x) for x in READ}]
        if diff or extra: bad.append("자료 목록 불일치: 해시 다름/목록 없음 %d (%s) · 목록에만 있음 %d (%s)" % (len(diff), diff[:3], len(extra), extra[:3]))
OUT["manifest_sha256"] = MAN_SHA

print("== 해부 23 %s— 위치 대조 N3 의 잡음 바닥 · 세트 %s" % ("파일럿 " if PILOT else "판정 ", " · ".join("%s %s" % (st["id"], st["name"]) for st in SETS)))
print("  배경 목록(규칙 = 하드코딩 · 각 20): S1 %s · S2 %s · S3 %s" % tuple("✅" if RULE[k] == h else "❌" for k, h in (("S1", SEEDS_S1), ("S2", SEEDS_S2), ("S3", SEEDS_S3))))
print("  읽은 파일 %d · 목록 해시 %s · 자료 목록 %s %s" % (len(READ), READ_DIGEST[:16], MANIFEST, (MAN_SHA[:16] if MAN_SHA else ("— (파일럿: 대조 안 함)" if PILOT else "없음 ❌"))))
print("  G0-b 같은 코드 · 같은 뒤 난수 → 궤적 동일: %d/%d (%s)" % (g0b_same, g0b_chk, G0B[0]))
for k, v in g0c.items(): print("  G0-c %s: N3 다시 셈 %d 쌍 · 발표 %d 쌍 · 최대 차 %.1e (%s)" % (k, v["pairs"], v["pub_pairs"], v["max_abs_diff"], v["file"]))
for x in notes: print("  ·", x)
print("  점검 문제 %d" % len(bad))
for x in bad[:20]: print("    ", x)
if bad and not PILOT:
    print("🚨 점검 실패 — 판정하지 않는다(판정값을 계산하지 않음)"); OUT["bad"] = bad
    json.dump(OUT, open(os.path.join(BASE, OUTFN + "_halted.json"), "w"), ensure_ascii=False, indent=1); sys.exit(1)

# ---------------- C0 · 잡음 바닥 · 분류 ----------------
CLS, F, ROWS = {}, [], {}
for st in SETS:
    print("\n== %s %s · 배경 %d" % (st["id"], st["name"], len(st["seeds"])))
    for b in st["bases"]:
        dat = DATA[(st["id"], b)]; seeds = sorted(dat["W"])
        if len(seeds) < 2:
            print("  %-9s 배경 %d — 계산 불가" % (b, len(seeds))); continue
        judged = dat["K"] == 0
        fv = floor_values(dat["W"], seeds); q = floor_q(dat["W"], seeds); qe = floor_q_emp(dat["W"], seeds)
        lo, hi = floor_q_boot(dat["W"], seeds, "floor|%s|%s" % (st["id"], b))
        rate, rate_e = false_fail_rate(dat["W"], seeds)
        exceed = sum(1 for v in fv if abs(v) > q) / len(fv)
        prs = same_loop_pairs(dat["codes"])
        nfail = 0; rows = []
        pv = pooled_var(dat["W"], seeds)
        for a, c in prs:
            dd = dd_mean(dat["D"], seeds, a, c); k = classify(dd, q, hi); fail = abs(dd) > TOL_NEG
            rv = var_ratio(dat["D"], seeds, a, c, pv); rob = k == "P" and classify_rob(dd, hi, rv)
            rows.append({"a": a, "c": c, "L": dat["codes"][a], "dd": dd, "fail": fail, "cls": k, "p_rob": rob, "var_ratio": rv})
            if judged:
                CLS[(st["id"], b, a, c)] = (k, dd, rob)
                if fail: F.append((st["id"], b, a, c, dd))
            nfail += fail
        crn, crn_lt1, crn_n = crn_ratio(dat["W"], dat["D"], seeds)
        cnt = {k: sum(1 for r in rows if r["fail"] and r["cls"] == k) for k in "NAP"}
        nrob = sum(1 for r in rows if r["fail"] and r["p_rob"])
        rv_nf = [r["var_ratio"] for r in rows if not r["fail"] and r["var_ratio"] == r["var_ratio"]]
        print("  %-9s C0 K %d/%d (최소 끝/처음 %.3f) %s\n            잡음 바닥 q95 %.3f [%.3f, %.3f] · 경험 q95 %.3f · 경험 초과율 %.3f · 거짓 실패율 %.3f(경험 %.3f · 기대 실패 %.1f) · N3 실패 %d/%d → N %d · A %d · P %d (P_rob %d)"
              "\n            CRN 비 중앙 %.3f(1 미만 코드 %.2f · %d 코드%s) · 실패 아닌 쌍 R_var 중앙 %s"
              % (b, dat["K"], len(seeds), dat["minr"], "판정" if judged else "⏭️ 판정 안 함(계측 불신)", q, lo, hi, qe, exceed, rate, rate_e,
                 rate * len(prs), nfail, len(prs), cnt["N"], cnt["A"], cnt["P"], nrob, crn, crn_lt1, crn_n,
                 " · 🟡 중앙 < %.1f · 1 미만 몫 ≥ 0.75: A 결합 남음" % CRN_READ if crn == crn and crn < CRN_READ and crn_lt1 >= 0.75 else "", ("%.2f" % pct(rv_nf, 0.5)) if rv_nf else "—"))
        for r in rows:
            if r["fail"]:
                print("       ❌N3 %-10s vs %-10s (%d틱) ΔΔ %+.3f → %s%s · R_var %.2f" % (r["a"], r["c"], r["L"], r["dd"], r["cls"], " · P_rob" if r["p_rob"] else "", r["var_ratio"]))
        ROWS["%s|%s" % (st["id"], b)] = {"judged": judged, "K": dat["K"], "min_ratio": dat["minr"], "q95": q, "q95_ci": [lo, hi],
                                         "q95_emp": qe, "emp_exceed": exceed, "pooled_var": pooled_var(dat["W"], seeds), "floor_pairs": fv,
                                         "floor_max": max(abs(v) for v in fv), "false_fail_rate": rate, "false_fail_rate_emp": rate_e, "expected_fail": rate * len(prs),
                                         "n_pairs": len(prs), "n_fail": nfail, "counts_fail": cnt, "n_fail_p_rob": nrob,
                                         "crn_ratio_median": crn, "crn_ratio_frac_lt1": crn_lt1, "crn_codes": crn_n,
                                         "var_ratio_median_nonfail": pct(rv_nf, 0.5) if rv_nf else None, "pairs": rows}
OUT["rows"] = ROWS

print("\n== 종합 (판정하는 (세트, 밑 코드)의 N3 실패 쌍 %d)" % len(F))
for sid in [st["id"] for st in SETS]:
    ks = [CLS[(s_, b, a, c)] for s_, b, a, c, _ in F if s_ == sid]
    kk = [x[0] for x in ks]
    print("  %s: 실패 %d → N %d · A %d · P %d (P_rob %d)" % (sid, len(ks), kk.count("N"), kk.count("A"), kk.count("P"), sum(1 for x in ks if x[2])))
v, repl, desc = verdict(F, CLS)
print("  ALL3(세 세트 모두 판정되는 같은 고리 배열) %d" % desc["n_all3"])
for arr, sign, sets in repl:
    print("  재현된 P_rob(세 세트): %s %s/%s · 부호 %s · 세트 %s" % (arr[0], arr[1], arr[2], "+" if sign > 0 else "−", " ".join(sets)))
for arr, sign, sets in desc["repl_p_not_rob"]:
    print("  세 세트 P 이나 P_rob 아님: %s %s/%s · 부호 %s" % (arr[0], arr[1], arr[2], "+" if sign > 0 else "−"))
for arr, sign, sets, sel in desc["weak_2sets"]:
    print("  서술 'V+ 약'(정확히 두 세트 P · 같은 부호): %s %s/%s · 부호 %s · P 세트 %s · F 에 넣은 세트 %s" % (arr[0], arr[1], arr[2], "+" if sign > 0 else "−", " ".join(sets), " ".join(sel) or "없음"))
print("  → %s%s" % (v, "  (🟡 파일럿 — 판정 아님)" if PILOT else ""))
OUT["F"] = [list(x) for x in F]; OUT["verdict"] = v; OUT["replicated_P_rob"] = [[list(a), s, ss] for a, s, ss in repl]; OUT["verdict_desc"] = desc
OUT["bad"] = bad; OUT["notes"] = notes
json.dump(OUT, open(os.path.join(BASE, OUTFN + ".json"), "w"), ensure_ascii=False, indent=1)
print("기록 → %s.json" % OUTFN)
