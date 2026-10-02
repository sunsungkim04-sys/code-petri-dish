#!/usr/bin/env python3
"""mini26_analyze.py — 해부 26 판정기 (사전등록 §4 · 결과 전 고정)

  python3 mini26_analyze.py pilot DIR              파일럿: 가드 · 시간 · 배경 조작 확인만. Δ · 계통 R 은 계산하지 않는다.
  python3 mini26_analyze.py judge DIR --repro D2 [--json OUT]  판정: G0 · G1 → P1~P4 · 서술.
      (--repro 는 필수 — 없으면 G0 가 멈춘다 · --json 기본값 = 이 파일 옆 _RESULT_mini26.json)

검토 개정(10-02 · 결과 전 · 판정선 숫자 불변 — 사전등록 §11):
  M1 '판정 불가' 칸(뺀 배경 > 2)이 P1 · P2 · P3 판정을 실제로 막는다 · 곡선에 NaN 이면 cens_flip 이 예외 →
     점추정은 판정 불가 · 재추출 행은 NaN(어느 사건에도 안 들어가고 분모에는 남는다 — 보수적).
  M3 P4 는 자격 배경 k ≥ MIN_BG_P4(15) 일 때만 판정 · 아니면 판정 불가.
  S1 S5 §5.5 의 세 양을 서술로(① T/R − L 구조상 0 · ② = P4a · ③ 불임 자식 비율 대 ln R 기울기 · μ 0.3%).
  S3 가임 개체만 센 Δ 표(서술 · 감도). S4 P(c_rf 유효) 인쇄. S6 n_inject · skipped 인쇄 + 배경 안 팔 사이 동일성(G0).
  N3 μ0 계수기 키 집합 동일성(G0). G0 에 --repro 필수 · G1 칸 뺀 배경 > 2 이면 멈춤.

지표 꼴은 petri/inv21_analyze.py 에서 글자 그대로 옮겼다:
  s = ln((m_T + 0.5)/m_0) − ln((n_T − m_T + 0.5)/(n_0 − m_0))   (m = 끼운 계통 · 혈통으로 · 불임 자손 포함, n = 산 개체 전부)
  Δ = s(mut) − 기준 팔 s 의 평균
  교차 = 평균 Δ 곡선(μ 순서)의 첫 + → ≤0 부호 변화의 선형 보간 (없으면 None)
재추출 = 배경 2,000 · 95% 백분위 · 난수는 ('haebu26', 칸 이름) 으로 키.
"""
import hashlib
import json
import math
import os
import glob
import sys

import numpy as np

N_BOOT = 2000
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
RULES = ["roll", "find"]
RHO_SCARCE, RHO_RICH = 8.0, 32.0
JUDGE_SEEDS = list(range(26001, 26031))      # 후보(시드 순) — 앞 20 자격 배경을 쓴다
PILOT_SEEDS = list(range(26901, 26911))
N_BG = 20
NREF = 5
CONFIG = {"K": 1000, "death": 1.0 / 450, "grow": 6000, "assay": 3000, "every": 250, "n0": 100,
          "frac": 0.1, "nref": NREF, "min_pop": 50}
VERSION = "mini26 v1.0.0"
CROSS_FACTOR = 2.0                            # P2 · P3 '크게 밀린다' = 교차가 두 배 이상(또는 범위 밖)
S_BAND = (math.log(0.8), math.log(1.25))      # P4b — 해부14 §4 A1b 후보선 1.25 를 대칭으로
MIN_POS = 100                                 # P4 — 배경마다 두 계통의 완료 위치 ≥ 100
MIN_BG_P4 = 15                                # P4 — 자격 배경 k ≥ 15(20 중) 아니면 판정 불가 (검토 M3 · 결과 전)
MAX_DROP = 2                                  # 칸에서 뺀 배경 > 2 → 그 칸 판정 불가 (§1-3 · 원래 규칙 · 검토 M1 로 실제 적용)
SLOPE_MU = 0.003                              # S1 ③ — S5 §5.5 의 μ 0.3%
SLOPE_BAND = (0.68 - 0.20, 0.68 + 0.20)       # S1 ③ — 원고의 접시 값 0.68 ± 0.20 (서술 · 판정 아님)
REF_KINDS = ["ref%d" % i for i in range(NREF)]


def rng_of(*parts):
    key = "|".join(str(p) for p in parts).encode()
    h = int.from_bytes(hashlib.sha256(key).digest()[:16], "little")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(h)))


def boot_idx(cell, n):
    return rng_of("haebu26", cell).integers(0, n, size=(N_BOOT, n))


def ci(vals, lo=2.5, hi=97.5):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:                            # 빈 재추출(서술 칸에서만 생길 수 있음) — 예외 대신 NaN
        return float("nan"), float("nan")
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))


def s_of(series):
    n0, m0 = series[0][1], series[0][2]
    nT, mT = series[-1][1], series[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def flip_of(curve):
    for (a, pa), (b, pb) in zip(zip(MUS, curve), zip(MUS[1:], curve[1:])):
        if pa > 0 >= pb:
            return a + (b - a) * pa / (pa - pb)
    return None


class Undecidable(Exception):
    pass


def cens_flip(curve):
    """검열 교차: Δ(0) ≤ 0 → 0 · 1% 까지 + → +inf · 그 밖 보간값.
    곡선에 NaN/inf 가 있으면 조용히 +inf 가 되던 경로(검토 M1b)를 막는다 → Undecidable."""
    if any(not math.isfinite(float(x)) for x in curve):
        raise Undecidable("곡선에 NaN — %s" % ["%.3g" % x for x in curve])
    if curve[0] <= 0:
        return 0.0
    f = flip_of(curve)
    return math.inf if f is None else f


def load(d, rho):
    F = {}
    for f in glob.glob(os.path.join(d, "mini26_r%g_s*.json" % rho)):
        z = json.load(open(f))
        F[z["config"]["seed"]] = z
    return F


def arm_map(z):
    return {(a["rule"], a["mu"], a["kind"]): a for a in z["arms"]}


# ───────────────────────── 파일럿 ─────────────────────────
def pilot(d):
    print("== 해부 26 파일럿 — 가드 · 시간 · 배경 조작 확인만 (Δ · 계통 R 은 계산하지 않는다)")
    bad = []
    for rho in (RHO_SCARCE, RHO_RICH):
        F = load(d, rho)
        for s, z in sorted(F.items()):
            if s not in PILOT_SEEDS:
                bad.append("파일럿 디렉터리에 파일럿 아닌 시드 %d" % s)
            g = z["grow"]
            if not g["cons_ok"]:
                bad.append("성장 보존 r%g s%d %s" % (rho, s, g["cons_bad"]))
            am = arm_map(z)
            for a in z["arms"]:
                if not a["cons_ok"]:
                    bad.append("팔 보존 r%g s%d %s %s %s" % (rho, s, a["rule"], a["mu"], a["kind"]))
            ident = 0
            for (rule, mu, kind), a in am.items():
                if rule == "roll" and mu == 0.0 and ("find", 0.0, kind) in am:
                    b = am[("find", 0.0, kind)]
                    same = a["series"] == b["series"]
                    ident += 1
                    if not same:
                        bad.append("μ0 규칙 동일성 깨짐 r%g s%d %s" % (rho, s, kind))
            gc = g["counters"].get("0_4", {})
            Rbg = gc["draws"] / gc["positions"] if gc.get("positions") else float("nan")
            stall = gc["stalls"] / gc["attempts"] if gc.get("attempts") else float("nan")
            ext = sum(1 for a in z["arms"] if a["extinct_tick"] >= 0)
            skp = sum(a["skipped"] for a in z["arms"])
            if len(set((a["n_inject"], a["skipped"]) for a in z["arms"])) > 1:
                bad.append("끼우기 수 · 건너뜀이 팔마다 다름 r%g s%d" % (rho, s))
            print("  r%-3g s%d  성장 n_end %4d · 점유 %4d/%d · 자격 %s · 보존 검사 %d · 배경 R(상주) %.2f · 멈춤 비 %.3f"
                  " · 팔 %d · μ0 동일성 대조 %d · 전멸 팔 %d · 끼우기 건너뜀 %d · %.0f초"
                  % (rho, s, g["n_end"], g["occ_end"], z["config"]["K"], g["qualified"], g["cons_checks"],
                     Rbg, stall, len(z["arms"]), ident, ext, skp, z["seconds"]))
    print("-- 문제 %d" % len(bad))
    for b in bad:
        print("   ", b)
    return 0 if not bad else 1


# ───────────────────────── 판정 ─────────────────────────
def gate(d, repro):
    bad = []
    BG = {}
    data = {}
    for rho in (RHO_SCARCE, RHO_RICH):
        F = load(d, rho)
        for s in F:
            if s in PILOT_SEEDS:
                bad.append("판정 디렉터리에 파일럿 시드 %d" % s)
            if s not in JUDGE_SEEDS:
                bad.append("후보 밖 시드 %d" % s)
        qual = []
        for s in JUDGE_SEEDS:
            if s not in F:
                break                          # 시드 순으로 끊김 없이 있어야 '앞 20' 이 정의된다
            z = F[s]
            if z["grow"]["qualified"]:
                qual.append(s)
            if len(qual) == N_BG:
                break
        if len(qual) < N_BG:
            bad.append("r%g 자격 배경 %d < %d — 후보를 시드 순으로 더 돌릴 것(26030 까지)" % (rho, len(qual), N_BG))
        BG[rho] = qual
        data[rho] = F
        for s in qual:
            z = F[s]
            if z["version"] != VERSION:
                bad.append("버전 r%g s%d %s" % (rho, s, z["version"]))
            for k, v in CONFIG.items():
                if abs(float(z["config"][k]) - float(v)) > 1e-12:
                    bad.append("설정 %s r%g s%d" % (k, rho, s))
            if [round(m, 6) for m in z["mus"]] != [round(m, 6) for m in MUS] or z["rules"] != RULES:
                bad.append("μ · 규칙 목록 r%g s%d" % (rho, s))
            g = z["grow"]
            if not g["cons_ok"] or g["cons_checks"] < CONFIG["grow"] // CONFIG["every"] + 2:
                bad.append("성장 보존 r%g s%d" % (rho, s))
            am = arm_map(z)
            if len(am) != len(RULES) * len(MUS) * (1 + NREF) or len(z["arms"]) != len(am):
                bad.append("팔 수 %d r%g s%d" % (len(z["arms"]), rho, s))
            for a in z["arms"]:
                need = (CONFIG["assay"] // CONFIG["every"] + 2) if a["extinct_tick"] < 0 else 2
                if not a["cons_ok"] or a["cons_checks"] < need:
                    bad.append("팔 보존 r%g s%d %s %g %s" % (rho, s, a["rule"], a["mu"], a["kind"]))
            for kind in ["mut"] + REF_KINDS:
                a, b = am.get(("roll", 0.0, kind)), am.get(("find", 0.0, kind))
                if a is None or b is None or a["series"] != b["series"]:
                    bad.append("μ0 규칙 동일성 r%g s%d %s" % (rho, s, kind))
                elif set(a["counters"]) != set(b["counters"]):
                    bad.append("μ0 계수기 키 집합 동일성 r%g s%d %s" % (rho, s, kind))
                elif any(a["counters"][k]["positions"] != b["counters"][k]["positions"] for k in a["counters"]):
                    bad.append("μ0 위치 수 동일성 r%g s%d %s" % (rho, s, kind))
            inj = set((a["n_inject"], a["skipped"]) for a in z["arms"])
            if len(inj) != 1:
                bad.append("끼우기 수 · 건너뜀이 팔마다 다름 r%g s%d %s" % (rho, s, sorted(inj)))
    if not repro:
        bad.append("재현 대조 없음 — judge 는 --repro 필수(§1-5 · G0)")
    else:
        n = 0
        for f in sorted(glob.glob(os.path.join(repro, "mini26_r*_s*.json"))):
            z2 = json.load(open(f))
            z1p = os.path.join(d, os.path.basename(f))
            if not os.path.exists(z1p):
                bad.append("재현 대조 원본 없음 %s" % os.path.basename(f))
                continue
            z1 = json.load(open(z1p))
            z1.pop("seconds", None)
            z2.pop("seconds", None)
            if json.dumps(z1, sort_keys=True) != json.dumps(z2, sort_keys=True):
                bad.append("재현 불일치 %s" % os.path.basename(f))
            n += 1
        if n < 2:
            bad.append("재현 대조 파일 %d < 2" % n)
    return bad, BG, data


def s_fert(series):
    """서술(S3): 가임 개체만 — m = 끼운 계통 가임(mf) · n = 가임 전부(nf)"""
    n0, m0 = series[0][4], series[0][3]
    nT, mT = series[-1][4], series[-1][3]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def delta_table(data, BG, rho, rule, sfun=None):
    """D[mu] = 배경 순서의 Δ (전멸 팔이 있으면 nan)"""
    sfun = sfun or s_of
    D = {mu: [] for mu in MUS}
    drop = {mu: 0 for mu in MUS}
    ext_lin = {mu: 0 for mu in MUS}
    for s in BG[rho]:
        am = arm_map(data[rho][s])
        for mu in MUS:
            arms = [am[(rule, mu, k)] for k in ["mut"] + REF_KINDS]
            if any(a["extinct_tick"] >= 0 or a["series"][-1][1] == 0 for a in arms):
                D[mu].append(float("nan"))
                drop[mu] += 1
                continue
            if arms[0]["series"][-1][2] == 0:
                ext_lin[mu] += 1
            base = sum(sfun(a["series"]) for a in arms[1:]) / NREF
            D[mu].append(sfun(arms[0]["series"]) - base)
    return {mu: np.array(v) for mu, v in D.items()}, drop, ext_lin


def neutral_table(data, BG, rho, rule, mu):
    out = []
    for s in BG[rho]:
        am = arm_map(data[rho][s])
        arms = [am[(rule, mu, k)] for k in REF_KINDS]
        if any(a["extinct_tick"] >= 0 or a["series"][-1][1] == 0 for a in arms):
            out.append(float("nan"))
            continue
        out.append(s_of(arms[0]["series"]) - sum(s_of(a["series"]) for a in arms[1:]) / (NREF - 1))
    return np.array(out)


def R_table(data, BG, rho, rule, mu, key):
    """끼운 팔(mut)에서 계통 key('1_2' 짧은 · '0_4' 상주)의 R · T/R · I 성분 — 배경 순서"""
    R, P, TR, E = [], [], [], []
    for s in BG[rho]:
        a = arm_map(data[rho][s])[(rule, mu, "mut")]
        c = a["counters"].get(key)
        if not c or c["positions"] == 0:
            R.append(float("nan")); P.append(0); TR.append(float("nan")); E.append(float("nan"))
            continue
        R.append(c["draws"] / c["positions"])
        P.append(c["positions"])
        TR.append(c["copy_ticks"] / c["attempts"] if c["attempts"] else float("nan"))
        E.append(c["err_realized"] / c["positions"])
    return np.array(R), np.array(P), np.array(TR), np.array(E)


def curve_boot(D, idx):
    """재추출마다 평균 Δ 곡선 → 검열 교차 (곡선에 NaN 이 있는 행은 NaN — 검토 M1b)"""
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            M = np.stack([np.nanmean(D[mu][idx], axis=1) for mu in MUS], axis=1)
    out = []
    for row in M:
        try:
            out.append(cens_flip(list(row)))
        except Undecidable:
            out.append(float("nan"))
    return np.array(out)


def qcens(v, q):
    """검열값(0 · inf)이 섞인 재추출의 분위 — 보간하지 않는다(inf − inf 방지) · NaN 행은 빼고(비율은 따로 인쇄)"""
    s = np.asarray(v, float)
    s = np.sort(s[~np.isnan(s)])
    if len(s) == 0:
        return float("nan")
    return float(s[min(len(s) - 1, max(0, int(math.ceil(q * len(s))) - 1))])


def fmt_c(c):
    if c is None or (isinstance(c, float) and math.isnan(c)):
        return "판정 불가(NaN)"
    return "+inf(1%까지 +)" if c == math.inf else ("0(Δ(0)≤0)" if c == 0 else "%.3f%%" % (100 * c))


def judge(d, repro, json_out=None):
    bad, BG, data = gate(d, repro)
    print("== 해부 26 — 두 가정 최소 모형 · 판정 (%s)" % VERSION)
    for rho in (RHO_SCARCE, RHO_RICH):
        print("  ρ %g 배경 %d: %s" % (rho, len(BG[rho]), " ".join(map(str, BG[rho]))))
    print("-- G0 점검 문제 %d" % len(bad))
    for b in bad[:40]:
        print("   ", b)
    if bad:
        print("** G0 실패 — 판정하지 않는다")
        return 2
    n = N_BG
    # ── G1 계측기 음성 대조 ──
    g1_bad = []
    for mu in (0.0, 0.01):
        v = neutral_table(data, BG, RHO_SCARCE, "roll", mu)
        g1_drop = int(np.isnan(v).sum())
        if g1_drop > MAX_DROP:
            print("-- G1 μ %.2f%% 칸 뺀 배경 %d > %d — 계측기 대조 판정 불가" % (100 * mu, g1_drop, MAX_DROP))
            g1_bad.append(mu)
            continue
        idx = boot_idx("G1|%g" % mu, n)
        bm = np.nanmean(v[idx], axis=1)
        lo, hi = ci(bm, 0.5, 99.5)
        ok = lo <= 0 <= hi
        print("-- G1 기준 팔끼리 Δ(ref0 대 나머지) roll · ρ8 · μ %.2f%%: %+.3f [99%% %+.3f, %+.3f] %s"
              % (100 * mu, np.nanmean(v), lo, hi, "✅" if ok else "❌"))
        if not ok:
            g1_bad.append(mu)
    if g1_bad:
        print("** G1 실패 — 계측기가 같은 코드끼리 0 을 못 낸다. 판정하지 않는다(사람 결정)")
        return 3

    res = {"BG": BG}
    # ── 서술(S6): 끼우기 수 · 건너뜀 ──
    for rho in (RHO_SCARCE, RHO_RICH):
        ni = [data[rho][s]["arms"][0]["n_inject"] for s in BG[rho]]
        sk = [data[rho][s]["arms"][0]["skipped"] for s in BG[rho]]
        print("-- 서술: ρ %g 끼운 개체 수(배경당) %d–%d · 끼우기 건너뜀 합 %d (배경 안 96 팔 사이 동일성은 G0 에서 확인)"
              % (rho, min(ni), max(ni), sum(sk)))
    # ── Δ 표 ──
    T = {}
    UNDEC = set()
    for rho in (RHO_SCARCE, RHO_RICH):
        for rule in RULES:
            D, drop, extl = delta_table(data, BG, rho, rule)
            T[(rho, rule)] = D
            print("\n-- Δ 표 · ρ %g · %s  (팔 전멸로 뺀 배경 · 끼운 계통 소멸 배경)" % (rho, rule))
            for mu in MUS:
                idx = boot_idx("D|%g|%s|%g" % (rho, rule, mu), n)
                bm = np.nanmean(D[mu][idx], axis=1)
                lo, hi = ci(bm)
                print("   %5.2f%% | %+.3f [%+.3f, %+.3f] | 뺀 %d · 소멸 %d"
                      % (100 * mu, np.nanmean(D[mu]), lo, hi, drop[mu], extl[mu]))
                res["D|%g|%s|%g" % (rho, rule, mu)] = [float(np.nanmean(D[mu])), lo, hi, drop[mu], extl[mu]]
                if drop[mu] > MAX_DROP:
                    print("   ** 이 칸 판정 불가(뺀 배경 > %d)" % MAX_DROP)
                    res["undecidable|%g|%s|%g" % (rho, rule, mu)] = True
                    UNDEC.add((rho, rule, mu))

    def mean_curve(D):
        return [float(np.nanmean(D[mu])) for mu in MUS]

    def undec_cells(conds, mus=MUS):
        return [(r, u, m) for (r, u) in conds for m in mus if (r, u, m) in UNDEC]

    def point_cross(Dx):
        try:
            return cens_flip(mean_curve(Dx))
        except Undecidable:
            return float("nan")

    # ── P1 ──
    D = T[(RHO_SCARCE, "roll")]
    r0, r1 = res["D|8|roll|0"], res["D|8|roll|0.01"]
    p1_und = undec_cells([(RHO_SCARCE, "roll")], (0.0, 0.01)) or not (np.isfinite(r0[1]) and np.isfinite(r1[2]))
    p1 = (not p1_und) and r0[1] > 0 and r1[2] < 0
    c_rf = point_cross(D)
    idx = boot_idx("P1|cross", n)
    cb = curve_boot(D, idx)
    p1_mark = "판정 불가(뺀 배경 > %d 또는 NaN)" % MAX_DROP if p1_und else ("✅" if p1 else "❌")
    print("\n-- P1 roll-first · 희소: Δ(0) 하한 %+.3f > 0 · Δ(1%%) 상한 %+.3f < 0 → %s"
          % (r0[1], r1[2], p1_mark))
    if p1_und:
        pass
    elif not p1:
        if r0[1] <= 0:
            print("   갈래: 속도 이점(Δ(0) 하한 > 0)부터 서지 않는다")
        else:
            print("   갈래: 1% 안에서 반전이 구간상 서지 않는다")
    print("   교차(서술 · 0.18%% 와 수치 일치는 판정하지 않는다): %s · 재추출 [%s, %s] · inf 비율 %.3f · 0 비율 %.3f · NaN 비율 %.3f"
          % (fmt_c(c_rf), fmt_c(qcens(cb, 0.025)), fmt_c(qcens(cb, 0.975)),
             float(np.mean(cb == math.inf)), float(np.mean(cb == 0)), float(np.mean(np.isnan(cb)))))
    res["P1"] = "undecidable" if p1_und else p1
    res["cross_roll_scarce"] = c_rf

    # ── P2 · P3 ──
    def compare(name, D_other, paired):
        idx_a = boot_idx("%s|a" % name, n)
        idx_b = idx_a if paired else boot_idx("%s|b" % name, n)
        ca = curve_boot(D, idx_a)
        cbb = curve_boot(D_other, idx_b)
        valid = np.isfinite(ca) & (ca > 0)
        pushed = valid & ((cbb == math.inf) | (cbb >= CROSS_FACTOR * ca))
        kept = valid & np.isfinite(cbb) & (cbb < CROSS_FACTOR * ca)
        p_push, p_keep = float(pushed.mean()), float(kept.mean())
        c_pt = point_cross(D_other)
        verdict = "✅ 사라지거나 크게 밀린다" if p_push >= 0.975 else ("❌ 교차가 남는다(두 배 안)" if p_keep >= 0.975 else "🟡 구별 안 됨")
        p_valid = float(valid.mean())
        nan_rows = int(np.isnan(ca).sum() + np.isnan(cbb).sum())
        return c_pt, p_push, p_keep, verdict, float(np.mean(cbb == math.inf)), p_valid, nan_rows

    for name, key, paired, label in (("P2", (RHO_SCARCE, "find"), True, "find-first · 희소 (배경 짝 재추출)"),
                                     ("P3", (RHO_RICH, "roll"), False, "roll-first · 풍부 (독립 재추출)")):
        c_pt, pp, pk, v, pinf, pval, nnan = compare(name, T[key], paired)
        und = undec_cells([(RHO_SCARCE, "roll"), key])
        if und or p1_und or (isinstance(c_rf, float) and math.isnan(c_rf)) or (isinstance(c_pt, float) and math.isnan(c_pt)):
            v = "판정 불가(칸 %s 뺀 배경 > %d · 또는 P1 판정 불가 · 또는 NaN) — 서술만: " % (
                ",".join("%g/%s/%g%%" % (r, u, 100 * m) for r, u, m in und) or "-", MAX_DROP) + v
        elif not p1:
            v = "판정 보류(P1 ❌) — 서술만: " + v
        print("\n-- %s %s: 교차 %s 대 roll-희소 %s · P(밀림 ≥ %.0f배 또는 없음) %.3f · P(두 배 안에 남음) %.3f · 재추출 inf 비율 %.3f → %s"
              % (name, label, fmt_c(c_pt), fmt_c(c_rf), CROSS_FACTOR, pp, pk, pinf, v))
        print("   P(c_rf 유효 = 유한 > 0) %.3f · NaN 재추출 행 %d (분모는 전체 %d · 유효 비율이 낮으면 🟡 는 c_rf 정밀도 탓 — §4 해석 규칙)"
              % (pval, nnan, N_BOOT))
        res[name] = {"cross": c_pt, "p_push": pp, "p_keep": pk, "verdict": v, "p_crf_valid": pval}

    # ── P4 ──
    R2, P2_, TR2, _ = R_table(data, BG, RHO_SCARCE, "roll", 0.0, "1_2")
    R4, P4_, TR4, _ = R_table(data, BG, RHO_SCARCE, "roll", 0.0, "0_4")
    ok = (P2_ >= MIN_POS) & (P4_ >= MIN_POS) & np.isfinite(R2) & np.isfinite(R4) & (R4 > 1) & (R2 > 1)
    k = int(ok.sum())
    if k < MIN_BG_P4:
        print("\n-- P4 roll-first · 희소 · μ 0 · 끼운 팔 안 (배경 %d/%d) — 자격 배경 < %d → 판정 불가" % (k, n, MIN_BG_P4))
        res["P4"] = {"k": k, "verdict": "undecidable"}
        p4_txt = "판정 불가(자격 배경 %d < %d)" % (k, MIN_BG_P4)
        return finish(res, p1, p1_und, p4_txt, json_out, data, BG, T)
    lnd = np.log(R2[ok]) - np.log(R4[ok])
    lns = np.log((R2[ok] - 1) * 2) - np.log((R4[ok] - 1) * 4)
    idx = boot_idx("P4", k)
    lo_a, hi_a = ci(lnd[idx].mean(axis=1))
    lo_b, hi_b = ci(lns[idx].mean(axis=1))
    p4a = lo_a > 0
    if S_BAND[0] <= lo_b and hi_b <= S_BAND[1]:
        p4b = "✅"
    elif hi_b < S_BAND[0] or lo_b > S_BAND[1]:
        p4b = "❌"
    else:
        p4b = "🟡"
    print("\n-- P4 roll-first · 희소 · μ 0 · 끼운 팔 안 (배경 %d/%d · 판정에 필요한 최소 %d)" % (k, n, MIN_BG_P4))
    print("   R 중앙 L2 %.2f · L4 %.2f · T/R 중앙 L2 %.3f · L4 %.3f"
          % (np.median(R2[ok]), np.median(R4[ok]), np.nanmedian(TR2[ok]), np.nanmedian(TR4[ok])))
    print("   P4a ln R(2) − ln R(4) = %+.3f [%+.3f, %+.3f] → %s" % (lnd.mean(), lo_a, hi_a, "✅" if p4a else "❌"))
    print("   P4b ln[S(2)/S(4)] (S = (R−1)·L) = %+.3f [%+.3f, %+.3f] · 비 %.3f [%.3f, %.3f] · 띠 [0.8, 1.25] → %s"
          " (서술: 독립 시도 기준 비 0.5)"
          % (lns.mean(), lo_b, hi_b, math.exp(lns.mean()), math.exp(lo_b), math.exp(hi_b), p4b))
    p4 = p4a and p4b == "✅"
    p4_txt = "✅" if p4 else ("❌" if (not p4a or p4b == "❌") else "🟡")
    print("   P4 → %s" % p4_txt)
    res["P4"] = {"k": k, "lnR": [float(lnd.mean()), lo_a, hi_a], "lnS": [float(lns.mean()), lo_b, hi_b], "p4b": p4b,
                 "verdict": p4_txt}
    return finish(res, p1, p1_und, p4_txt, json_out, data, BG, T)


def finish(res, p1, p1_und, p4_txt, json_out, data, BG, T):
    n = N_BG

    # ── 서술: R · I 사다리 ──
    print("\n-- 서술: 끼운 팔 안 계통별 R(추첨/완료 위치) 중앙 · 실현 오류/위치 ÷ μ (I) 중앙")
    for rho in (RHO_SCARCE, RHO_RICH):
        for rule in RULES:
            print("   ρ %g · %s" % (rho, rule))
            for mu in MUS:
                a2, _, _, e2 = R_table(data, BG, rho, rule, mu, "1_2")
                a4, _, _, e4 = R_table(data, BG, rho, rule, mu, "0_4")
                I2 = np.nanmedian(e2) / mu if mu > 0 else float("nan")
                I4 = np.nanmedian(e4) / mu if mu > 0 else float("nan")
                print("     %5.2f%% | R L2 %6.2f · L4 %6.2f | I L2 %6.2f · L4 %6.2f"
                      % (100 * mu, np.nanmedian(a2), np.nanmedian(a4), I2, I4))
    # ── 서술: 배경 조작 확인 ──
    print("\n-- 서술: 배경(성장 끝 · 상주만) 조작 확인 — 규칙: 희소 중앙 R ≥ 2 · 풍부 중앙 R ≤ 1.1 이어야 P3 를 '재료 풍부' 로 읽는다")
    man = {}
    for rho in (RHO_SCARCE, RHO_RICH):
        Rs, ne = [], []
        for s in BG[rho]:
            g = data[rho][s]["grow"]
            c = g["counters"]["0_4"]
            Rs.append(c["draws"] / c["positions"])
            ne.append(g["n_end"])
        man[rho] = float(np.median(Rs))
        print("   ρ %g: 배경 R 중앙 %.2f [%.2f–%.2f] · 개체 수 중앙 %d" % (rho, man[rho], min(Rs), max(Rs), int(np.median(ne))))
    manip = man[RHO_SCARCE] >= 2 and man[RHO_RICH] <= 1.1
    print("   → %s" % ("조작 확인 ✅" if manip else "조작 확인 ❌ — P3 를 '재료 풍부' 의 결과로 읽지 않는다"))
    res["manip"] = {"R_scarce": man[RHO_SCARCE], "R_rich": man[RHO_RICH], "ok": manip}

    # ── 서술(S1): S5 §5.5 의 세 양 ──
    print("\n-- 서술(S1 · 판정 아님): S5 §5.5 '어떤 재구현이든' 의 세 양")
    dev = []
    for rho in (RHO_SCARCE, RHO_RICH):
        for rule in RULES:
            for mu in MUS:
                for key, L in (("1_2", 2), ("0_4", 4)):
                    _, _, tr, _ = R_table(data, BG, rho, rule, mu, key)
                    dev += list(tr[np.isfinite(tr)] - L)
    dev = np.array(dev)
    print("   ① T/R − L: 구조상 0(mini26.py 가 시도마다 copy_ticks += L) · 계기 점검만 — 전 칸 최소 %+.4f · 최대 %+.4f · 점 %d"
          % (dev.min(), dev.max(), len(dev)))
    print("   ② 이웃 고리 ln 추첨/글자 하한 > 0 = P4a 와 같은 양(위 P4a 줄)")
    xs, ys, bgi, nz = [], [], [], 0
    for j, s in enumerate(BG[RHO_SCARCE]):
        a = arm_map(data[RHO_SCARCE][s])[("roll", SLOPE_MU, "mut")]
        for key in ("1_2", "0_4"):
            c = a["counters"].get(key)
            if not c or c["positions"] == 0 or (c["kids_ok"] + c["kids_sterile"]) == 0:
                continue
            fr = c["kids_sterile"] / (c["kids_ok"] + c["kids_sterile"])
            if fr <= 0:
                nz += 1
                continue
            xs.append(math.log(c["draws"] / c["positions"]))
            ys.append(math.log(fr))
            bgi.append(j)
    xs, ys, bgi = np.array(xs), np.array(ys), np.array(bgi)

    def slope(sel):
        x, y = xs[sel], ys[sel]
        if len(x) < 3 or np.ptp(x) == 0:
            return float("nan")
        return float(np.polyfit(x, y, 1)[0])
    if len(xs) >= 3:
        b0 = slope(np.ones(len(xs), bool))
        bi = boot_idx("S1|slope", n)
        bs = []
        for row in bi:
            sel = np.concatenate([np.flatnonzero(bgi == j) for j in row])
            bs.append(slope(sel) if len(sel) else float("nan"))
        bs = np.array(bs)
        lo, hi = ci(bs)
        inside = SLOPE_BAND[0] <= lo and hi <= SLOPE_BAND[1]
        print("   ③ roll · ρ8 · μ %.1f%% · 끼운 팔: ln(불임 자식 비율) 대 ln R 기울기(계통 × 배경 점 %d · 불임 0 이라 뺀 점 %d) = %+.3f [%+.3f, %+.3f]"
              " · 원고 접시 값 0.68 ± 0.20 띠 안: %s (서술)"
              % (100 * SLOPE_MU, len(xs), nz, b0, lo, hi, "예" if inside else "아니오"))
        res["S1_slope"] = [b0, lo, hi, inside]
    else:
        print("   ③ 점 %d < 3 — 계산하지 않음" % len(xs))
    # ── 서술(S3): 가임만 센 Δ ──
    print("\n-- 서술(S3 · 감도 · 판정 아님): 가임 개체만 센 Δ(m = 끼운 계통 가임 · n = 가임 전부) — 불임 자손을 m 에 넣는 본 정의는 반전에 불리한 쪽")
    for rho in (RHO_SCARCE, RHO_RICH):
        for rule in RULES:
            Df, _, _ = delta_table(data, BG, rho, rule, s_fert)
            row = " ".join("%+.3f" % float(np.nanmean(Df[mu])) for mu in MUS)
            cf = float("nan")
            try:
                cf = cens_flip([float(np.nanmean(Df[mu])) for mu in MUS])
            except Undecidable:
                pass
            print("   ρ %g · %s | %s | 교차 %s" % (rho, rule, row, fmt_c(cf)))
            res["S3|%g|%s" % (rho, rule)] = [float(np.nanmean(Df[mu])) for mu in MUS]

    def mk(v):
        return v if isinstance(v, str) else ("✅" if v else "❌")
    p1_s = "판정 불가" if p1_und else mk(p1)
    print("\n결과 문장: 최소 모형에서 P1 %s · P2 %s · P3 %s · P4 %s"
          % (p1_s, res["P2"]["verdict"], res["P3"]["verdict"], p4_txt))
    out = json_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "_RESULT_mini26.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda x: None if (isinstance(x, float) and not math.isfinite(x)) else str(x))
    return 0


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    mode, d = sys.argv[1], sys.argv[2]
    repro = None
    if "--repro" in sys.argv:
        repro = sys.argv[sys.argv.index("--repro") + 1]
    json_out = None
    if "--json" in sys.argv:
        json_out = sys.argv[sys.argv.index("--json") + 1]
    if mode == "pilot":
        return pilot(d)
    if mode == "judge":
        return judge(d, repro, json_out)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
