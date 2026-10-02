#!/usr/bin/env python3
"""mini27_analyze.py — 해부 27 판정기 (사전등록 §4 · 결과 전 고정)

  python3 mini27_analyze.py pilot DIR
      파일럿: 가드 · 시간 · 배경 조작 확인만. Δ · 계통 R 은 계산하지 않는다.
  python3 mini27_analyze.py ident IDENT_DIR MINI26_DIR
      G2a 만: mini27 glob 경로(시드 26001 · 원 사다리)가 mini26 출력과 바이트 같은가.
  python3 mini27_analyze.py judge DIR --repro D2 --ident D3 --mini26 D4 [--json OUT]
      판정: G0 · G2a · G1 → P1 · G2b · P3 · P4 · 해석 열쇠 K1 · K2 · 서술.

mini26_analyze.py(a4d4ef2b…)를 고친 것. 지표 꼴은 그대로:
  s = ln((m_T + 0.5)/m_0) − ln((n_T − m_T + 0.5)/(n_0 − m_0))   (m = 끼운 계통 · 혈통 · 불임 자손 포함, n = 산 개체 전부)
  Δ = s(mut) − 기준 팔 다섯 s 평균
  교차 = 주어진 사다리 위 평균 Δ 곡선의 첫 + → ≤0 부호 변화의 선형 보간 · 검열(Δ(0) ≤ 0 → 0 · 끝까지 + → +inf)
재추출 = 배경 2,000 · 95% 백분위(G1 만 99%) · 난수 ('haebu27', 칸 이름).
v1.1 (독립 검토 반영 · 동결 전 · 파일럿 판정선 불변): M1 열쇠 셋 갈래(반전 · 반전 없음 · 불확정) · M3 실행 로그 해시·판 묶기 ·
S1 D(나) 문장 낮춤 · S2 'E/D 구분 불가' · §5 G2b 🟡 꼬리 · 오류 분할 단서 줄.
"""
import hashlib
import json
import math
import os
import re
import glob
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore", message="Mean of empty slice")   # NaN 칸은 판정 불가로 따로 센다

N_BOOT = 2000
MUS_ORIG = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
MUS_EXT = [0.0, 0.001625, 0.00325, 0.004875, 0.0065, 0.008125, 0.0121875, 0.01625]
MUS_ALL = sorted(set(MUS_ORIG) | set(MUS_EXT))
MU_TOP = MUS_EXT[-1]                          # 1.625% = 접시 글자변경 척도 1%
SCALE = 52.0 / 32.0                           # mini μ ÷ SCALE = 접시 글자변경 척도
RULES = ["roll", "find"]
GEOMS = ["glob", "loc", "place"]
RHO_SCARCE, RHO_RICH = 8.0, 32.0
JUDGE_SEEDS = list(range(27001, 27031))
PILOT_SEEDS = list(range(27901, 27911))
IDENT_SEED = 26001
N_BG = 20
NREF = 5
CONFIG = {"K": 1000, "death": 1.0 / 450, "grow": 6000, "assay": 3000, "every": 250, "n0": 100,
          "frac": 0.1, "nref": NREF, "min_pop": 50}
MAT_RADIUS = "1"
GRID = [40, 25]
VERSION = "mini27 v1.0.0"                     # mini27.py 판(출력 JSON 의 version) — 판정기 판은 v1.1
N_ARMS = {"glob": 2 * 15 * (1 + NREF), "loc": 2 * 15 * (1 + NREF), "place": 9 * (1 + NREF)}
CROSS_FACTOR = 2.0
S_BAND = (math.log(0.8), math.log(1.25))      # P4b — 해부26 P4b 와 같은 띠
MIN_POS = 100
MIN_BG_P4 = 15
MAX_DROP = 2
REF_KINDS = ["ref%d" % i for i in range(NREF)]
MINI26_RESULT = "_RESULT_mini26.json"
MINI26_RESULT_SHA = "e863ec859275dce4326ab7e4677e272be33658cd86dc9f3a48c3fe566fc4b43d"   # 해부26 §10 · §12
K_NONE_INF = 0.975                            # M1 — 열쇠 '반전 없음': Δ(1.625%) 하한 > 0 또는 교차 재추출 inf 비율 ≥ 이것
HERE = os.path.dirname(os.path.abspath(__file__))


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def run_log(dirpath):
    """M3 — 발사기 로그(DIR.log) 1~2행: python · numpy 판과 mini27.py sha256 · 끝 줄 ERROR 수"""
    p = dirpath.rstrip("/") + ".log"
    if not os.path.exists(p):
        return None
    L = open(p, encoding="utf-8").read().splitlines()
    m = re.search(r"python (\S+) numpy (\S+)", L[0]) if L else None
    sha = L[1].split()[0] if len(L) > 1 and L[1].split() else ""
    end = re.search(r"ERROR (\d+)\s*$", L[-1]) if L else None
    return {"log": p, "py": m.group(1) if m else None, "np": m.group(2) if m else None, "sha": sha,
            "err": int(end.group(1)) if (end and L[-1].startswith("===== 끝")) else None}


def log_guard(dirs, bad):
    """M3 — 출력 디렉터리마다 로그의 mini27.py 해시 = 지금 판정기 옆 mini27.py · python · numpy 판이 모두 같음 · ERROR 0"""
    want = sha_file(os.path.join(HERE, "mini27.py"))
    seen = {}
    for nm, dd in dirs:
        if not dd:
            continue
        info = run_log(dd)
        if info is None:
            bad.append("M3 실행 로그 없음 %s(%s.log)" % (nm, dd.rstrip("/")))
            continue
        if info["sha"] != want:
            bad.append("M3 %s 로그의 mini27.py 해시 %s… ≠ 지금 %s…" % (nm, info["sha"][:8], want[:8]))
        if info["py"] is None or info["np"] is None:
            bad.append("M3 %s 로그 첫 줄에 python · numpy 판 없음" % nm)
        if info["err"] != 0:
            bad.append("M3 %s 로그 끝 줄 ERROR %s(끝 줄 없음 = None)" % (nm, info["err"]))
        seen[nm] = (info["py"], info["np"])
    if len(set(seen.values())) > 1:
        bad.append("M3 python · numpy 판이 실행마다 다름 %s — 재현 · G2a 를 본 실행과 같은 기계 · 같은 python 에서" % seen)
    return seen


def rng_of(*parts):
    key = "|".join(str(p) for p in parts).encode()
    h = int.from_bytes(hashlib.sha256(key).digest()[:16], "little")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(h)))


def boot_idx(cell, n):
    return rng_of("haebu27", cell).integers(0, n, size=(N_BOOT, n))


def ci(vals, lo=2.5, hi=97.5):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan"), float("nan")
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))


def s_of(series):
    n0, m0 = series[0][1], series[0][2]
    nT, mT = series[-1][1], series[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def s_fert(series):
    n0, m0 = series[0][4], series[0][3]
    nT, mT = series[-1][4], series[-1][3]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def flip_of(curve, mus):
    for (a, pa), (b, pb) in zip(zip(mus, curve), zip(mus[1:], curve[1:])):
        if pa > 0 >= pb:
            return a + (b - a) * pa / (pa - pb)
    return None


class Undecidable(Exception):
    pass


def cens_flip(curve, mus):
    if any(not math.isfinite(float(x)) for x in curve):
        raise Undecidable("곡선에 NaN — %s" % ["%.3g" % x for x in curve])
    if curve[0] <= 0:
        return 0.0
    f = flip_of(curve, mus)
    return math.inf if f is None else f


def load(d, rho):
    F = {}
    for f in glob.glob(os.path.join(d, "mini27_r%g_s*.json" % rho)):
        z = json.load(open(f))
        F[z["config"]["seed"]] = z
    return F


def arm_map(z, geom):
    return {(a["rule"], a["mu"], a["kind"]): a for a in z["by_geom"][geom]["arms"]}


def strip_seconds(z):
    z = json.loads(json.dumps(z))
    z.pop("seconds", None)
    for G in z.get("by_geom", {}).values():
        G.pop("seconds", None)
    return z


# ───────────────────────── 기하별 공통 가드 ─────────────────────────
def geom_guards(z, rho, s, geom, bad, strict):
    G = z["by_geom"][geom]
    g = G["grow"]
    if not g["cons_ok"] or (strict and g["cons_checks"] < CONFIG["grow"] // CONFIG["every"] + 2):
        bad.append("성장 보존 %s r%g s%d %s" % (geom, rho, s, g["cons_bad"]))
    # 처음 배치 건너뜀(init_skipped · 이웃에 글자가 모자란 처음 개체)은 성장 6,000틱 앞의 일이라 멈추지 않고 서술로만 센다
    am = arm_map(z, geom)
    if strict and (len(am) != N_ARMS[geom] or len(G["arms"]) != len(am)):
        bad.append("팔 수 %d %s r%g s%d" % (len(G["arms"]), geom, rho, s))
    for a in G["arms"]:
        need = (CONFIG["assay"] // CONFIG["every"] + 2) if a["extinct_tick"] < 0 else 2
        if not a["cons_ok"] or (strict and a["cons_checks"] < need):
            bad.append("팔 보존 %s r%g s%d %s %g %s %s" % (geom, rho, s, a["rule"], a["mu"], a["kind"], a["cons_bad"]))
    ident = 0
    for kind in ["mut"] + REF_KINDS:
        a, b = am.get(("roll", 0.0, kind)), am.get(("find", 0.0, kind))
        if a is None or b is None:
            if strict:
                bad.append("μ0 규칙 동일성 팔 없음 %s r%g s%d %s" % (geom, rho, s, kind))
            continue
        ident += 1
        if a["series"] != b["series"]:
            bad.append("μ0 규칙 동일성 %s r%g s%d %s" % (geom, rho, s, kind))
        elif set(a["counters"]) != set(b["counters"]):
            bad.append("μ0 계수기 키 집합 동일성 %s r%g s%d %s" % (geom, rho, s, kind))
        elif any(a["counters"][k]["positions"] != b["counters"][k]["positions"] for k in a["counters"]):
            bad.append("μ0 위치 수 동일성 %s r%g s%d %s" % (geom, rho, s, kind))
    inj = set((a["n_inject"], a["skipped"]) for a in G["arms"])
    if len(inj) > 1:
        bad.append("끼우기 수 · 건너뜀이 팔마다 다름 %s r%g s%d %s" % (geom, rho, s, sorted(inj)))
    return ident


# ───────────────────────── 파일럿 ─────────────────────────
def pilot(d):
    print("== 해부 27 파일럿 — 가드 · 시간 · 배경 조작 확인만 (Δ · 계통 R 은 계산하지 않는다)")
    bad = []
    for rho in (RHO_SCARCE, RHO_RICH):
        F = load(d, rho)
        for s, z in sorted(F.items()):
            if s not in PILOT_SEEDS:
                bad.append("파일럿 디렉터리에 파일럿 아닌 시드 %d" % s)
            for geom in z["geoms"]:
                G = z["by_geom"][geom]
                g = G["grow"]
                ident = geom_guards(z, rho, s, geom, bad, strict=True)
                gc = g["counters"].get("0_4", {})
                Rbg = gc["draws"] / gc["positions"] if gc.get("positions") else float("nan")
                stall = gc["stalls"] / gc["attempts"] if gc.get("attempts") else float("nan")
                ext = sum(1 for a in G["arms"] if a["extinct_tick"] >= 0)
                ni = G["arms"][0]["n_inject"] if G["arms"] else -1
                skp = sum(a["skipped"] for a in G["arms"])
                print("  r%-3g s%d %-5s 성장 n_end %4d · 점유 %4d/%d · 자격 %s · 보존 검사 %d · 배경 R(상주) %.2f · 멈춤 비 %.3f"
                      " · 팔 %d · μ0 동일성 대조 %d · 전멸 팔 %d · 끼운 %d · 건너뜀 %d · 처음 배치 건너뜀 %d · %.0f초"
                      % (rho, s, geom, g["n_end"], g["occ_end"], z["config"]["K"], g["qualified"], g["cons_checks"],
                         Rbg, stall, len(G["arms"]), ident, ext, ni, skp, g.get("init_skipped", 0), G["seconds"]))
    logs = log_guard([("파일럿", d)], bad)
    print("-- 실행 로그(M3): %s" % logs)
    print("-- 문제 %d" % len(bad))
    for b in bad:
        print("   ", b)
    return 0 if not bad else 1


# ───────────────────────── G2a ─────────────────────────
def ident_check(ident_dir, mini26_dir):
    """mini27 glob 경로(시드 26001 · 원 사다리 · 기준 5)가 mini26 판정 자료와 바이트 같은가 — 성장 · 팔 전부"""
    bad = []
    n = 0
    for rho in (RHO_SCARCE, RHO_RICH):
        p27 = os.path.join(ident_dir, "mini27_r%g_s%d.json" % (rho, IDENT_SEED))
        p26 = os.path.join(mini26_dir, "mini26_r%g_s%d.json" % (rho, IDENT_SEED))
        if not (os.path.exists(p27) and os.path.exists(p26)):
            bad.append("G2a 파일 없음 r%g" % rho)
            continue
        z7, z6 = json.load(open(p27)), json.load(open(p26))
        G = z7["by_geom"].get("glob")
        if G is None or z7["geoms"] != ["glob"]:
            bad.append("G2a glob 만 있어야 함 r%g" % rho)
            continue
        if [round(m, 9) for m in z7["mus"]] != [round(m, 9) for m in MUS_ORIG]:
            bad.append("G2a μ 목록 r%g" % rho)
        for k in CONFIG:
            if abs(float(z7["config"][k]) - float(z6["config"][k])) > 1e-12:
                bad.append("G2a 설정 %s r%g" % (k, rho))
        if json.dumps(G["grow"], sort_keys=True) != json.dumps(z6["grow"], sort_keys=True):
            bad.append("G2a 성장 다름 r%g" % rho)
        if json.dumps(G["arms"], sort_keys=True) != json.dumps(z6["arms"], sort_keys=True):
            bad.append("G2a 팔 다름 r%g (%d 대 %d)" % (rho, len(G["arms"]), len(z6["arms"])))
        else:
            n += len(G["arms"])
    return bad, n


# ───────────────────────── 판정 ─────────────────────────
def gate(d, repro, ident, mini26):
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
                break
            z = F[s]
            if all(z["by_geom"][g]["grow"]["qualified"] for g in GEOMS if g in z["by_geom"]) and \
                    all(g in z["by_geom"] for g in GEOMS):
                qual.append(s)
            if len(qual) == N_BG:
                break
        if len(qual) < N_BG:
            bad.append("r%g 자격 배경(세 기하 모두) %d < %d — 후보를 시드 순으로 더 돌릴 것(27030 까지)" % (rho, len(qual), N_BG))
        BG[rho] = qual
        data[rho] = F
        for s in qual:
            z = F[s]
            if z["version"] != VERSION:
                bad.append("버전 r%g s%d %s" % (rho, s, z["version"]))
            for k, v in CONFIG.items():
                if abs(float(z["config"][k]) - float(v)) > 1e-12:
                    bad.append("설정 %s r%g s%d" % (k, rho, s))
            if str(z["config"]["mat_radius"]) != MAT_RADIUS or z["grid"] != GRID:
                bad.append("재료 반경 · 격자 r%g s%d" % (rho, s))
            if [round(m, 9) for m in z["mus"]] != [round(m, 9) for m in MUS_ALL] or z["rules"] != RULES \
                    or z["geoms"] != GEOMS:
                bad.append("μ · 규칙 · 기하 목록 r%g s%d" % (rho, s))
            for geom in GEOMS:
                geom_guards(z, rho, s, geom, bad, strict=True)
    if not repro:
        bad.append("재현 대조 없음 — judge 는 --repro 필수")
    else:
        n = 0
        for f in sorted(glob.glob(os.path.join(repro, "mini27_r*_s*.json"))):
            z2 = json.load(open(f))
            z1p = os.path.join(d, os.path.basename(f))
            if not os.path.exists(z1p):
                bad.append("재현 대조 원본 없음 %s" % os.path.basename(f))
                continue
            z1 = json.load(open(z1p))
            if json.dumps(strip_seconds(z1), sort_keys=True) != json.dumps(strip_seconds(z2), sort_keys=True):
                bad.append("재현 불일치 %s" % os.path.basename(f))
            n += 1
        if n < 2:
            bad.append("재현 대조 파일 %d < 2" % n)
    if not (ident and mini26):
        bad.append("G2a 대조 없음 — judge 는 --ident · --mini26 필수")
        n_id = 0
    else:
        b2, n_id = ident_check(ident, mini26)
        bad += b2
    logs = log_guard([("본 실행", d), ("재현", repro), ("G2a", ident)], bad)
    print("  실행 로그(M3): %s" % " · ".join("%s python %s numpy %s" % (k, v[0], v[1]) for k, v in logs.items()))
    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), MINI26_RESULT)
    if not os.path.exists(ref) or hashlib.sha256(open(ref, "rb").read()).hexdigest() != MINI26_RESULT_SHA:
        bad.append("해부26 결과 파일 해시 다름 · 없음 %s" % ref)
        r26 = None
    else:
        r26 = json.load(open(ref))
    return bad, BG, data, n_id, r26


def delta_table(data, BG, rho, geom, rule, mus, sfun=None):
    sfun = sfun or s_of
    D = {mu: [] for mu in mus}
    drop = {mu: 0 for mu in mus}
    ext_lin = {mu: 0 for mu in mus}
    for s in BG[rho]:
        am = arm_map(data[rho][s], geom)
        for mu in mus:
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


def neutral_table(data, BG, rho, geom, rule, mu):
    out = []
    for s in BG[rho]:
        am = arm_map(data[rho][s], geom)
        arms = [am[(rule, mu, k)] for k in REF_KINDS]
        if any(a["extinct_tick"] >= 0 or a["series"][-1][1] == 0 for a in arms):
            out.append(float("nan"))
            continue
        out.append(s_of(arms[0]["series"]) - sum(s_of(a["series"]) for a in arms[1:]) / (NREF - 1))
    return np.array(out)


def R_table(data, BG, rho, geom, rule, mu, key):
    R, P, TR, E = [], [], [], []
    for s in BG[rho]:
        a = arm_map(data[rho][s], geom)[(rule, mu, "mut")]
        c = a["counters"].get(key)
        if not c or c["positions"] == 0:
            R.append(float("nan")); P.append(0); TR.append(float("nan")); E.append(float("nan"))
            continue
        R.append(c["draws"] / c["positions"])
        P.append(c["positions"])
        TR.append(c["copy_ticks"] / c["attempts"] if c["attempts"] else float("nan"))
        E.append(c["err_realized"] / c["positions"])
    return np.array(R), np.array(P), np.array(TR), np.array(E)


def curve_boot(D, idx, mus):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        M = np.stack([np.nanmean(D[mu][idx], axis=1) for mu in mus], axis=1)
    out = []
    for row in M:
        try:
            out.append(cens_flip(list(row), mus))
        except Undecidable:
            out.append(float("nan"))
    return np.array(out)


def qcens(v, q):
    s = np.asarray(v, float)
    s = np.sort(s[~np.isnan(s)])
    if len(s) == 0:
        return float("nan")
    return float(s[min(len(s) - 1, max(0, int(math.ceil(q * len(s))) - 1))])


def fmt_c(c, top=MU_TOP):
    if c is None or (isinstance(c, float) and math.isnan(c)):
        return "판정 불가(NaN)"
    if c == math.inf:
        return "+inf(%.4g%%까지 +)" % (100 * top)
    if c == 0:
        return "0(Δ(0)≤0)"
    return "%.3f%%" % (100 * c)


def p4_stats(data, BG, rho, geom):
    R2, P2_, TR2, _ = R_table(data, BG, rho, geom, "roll", 0.0, "1_2")
    R4, P4_, TR4, _ = R_table(data, BG, rho, geom, "roll", 0.0, "0_4")
    ok = (P2_ >= MIN_POS) & (P4_ >= MIN_POS) & np.isfinite(R2) & np.isfinite(R4) & (R4 > 1) & (R2 > 1)
    lnd = np.log(R2[ok]) - np.log(R4[ok])
    lns = np.log((R2[ok] - 1) * 2) - np.log((R4[ok] - 1) * 4)
    return {"k": int(ok.sum()), "lnd": lnd, "lns": lns, "R2": R2[ok], "R4": R4[ok], "TR2": TR2[ok], "TR4": TR4[ok]}


def judge(d, repro, ident, mini26, json_out=None):
    bad, BG, data, n_id, r26 = gate(d, repro, ident, mini26)
    print("== 해부 27 — 최소 모형 + 국소 결핍 · 판정 (%s)" % VERSION)
    for rho in (RHO_SCARCE, RHO_RICH):
        print("  ρ %g 배경 %d: %s" % (rho, len(BG[rho]), " ".join(map(str, BG[rho]))))
    print("-- G0 · G2a 점검 문제 %d (G2a 바이트 같은 팔 %d)" % (len(bad), n_id))
    for b in bad[:40]:
        print("   ", b)
    if bad:
        print("** G0 · G2a 실패 — 판정하지 않는다")
        return 2
    n = N_BG
    # ── G1 ──
    g1_bad = []
    for mu in (0.0, MU_TOP):
        v = neutral_table(data, BG, RHO_SCARCE, "loc", "roll", mu)
        g1_drop = int(np.isnan(v).sum())
        if g1_drop > MAX_DROP:
            print("-- G1 loc μ %.4g%% 칸 뺀 배경 %d > %d — 판정 불가" % (100 * mu, g1_drop, MAX_DROP))
            g1_bad.append(mu)
            continue
        idx = boot_idx("G1|%g" % mu, n)
        bm = np.nanmean(v[idx], axis=1)
        lo, hi = ci(bm, 0.5, 99.5)
        ok = lo <= 0 <= hi
        print("-- G1 기준 팔끼리 Δ(ref0 대 나머지) loc · roll · ρ8 · μ %.4g%%: %+.3f [99%% %+.3f, %+.3f] %s"
              % (100 * mu, np.nanmean(v), lo, hi, "✅" if ok else "❌"))
        if not ok:
            g1_bad.append(mu)
    if g1_bad:
        print("** G1 실패 — 판정하지 않는다(사람 결정)")
        return 3

    res = {"BG": BG}
    for rho in (RHO_SCARCE, RHO_RICH):
        for geom in GEOMS:
            ni = [data[rho][s]["by_geom"][geom]["arms"][0]["n_inject"] for s in BG[rho]]
            sk = [data[rho][s]["by_geom"][geom]["arms"][0]["skipped"] for s in BG[rho]]
            isk = sum(data[rho][s]["by_geom"][geom]["grow"].get("init_skipped", 0) for s in BG[rho])
            print("-- 서술: ρ %g · %s 끼운 개체 수(배경당) %d–%d · 끼우기 건너뜀 합 %d · 처음 배치 건너뜀 합 %d"
                  % (rho, geom, min(ni), max(ni), sum(sk), isk))
    # ── Δ 표 ──
    T = {}
    UNDEC = set()
    plan = [(g, r) for g in ("loc", "glob") for r in RULES] + [("place", "roll")]
    for rho in (RHO_SCARCE, RHO_RICH):
        for geom, rule in plan:
            mus = MUS_EXT if geom == "place" else MUS_ALL
            D, drop, extl = delta_table(data, BG, rho, geom, rule, mus)
            T[(rho, geom, rule)] = D
            print("\n-- Δ 표 · ρ %g · %s · %s  (사다리 O = 원 · E = 확장 · 팔 전멸로 뺀 배경 · 끼운 계통 소멸)" % (rho, geom, rule))
            for mu in mus:
                idx = boot_idx("D|%g|%s|%s|%g" % (rho, geom, rule, mu), n)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    bm = np.nanmean(D[mu][idx], axis=1)
                lo, hi = ci(bm)
                lad = ("O" if mu in MUS_ORIG else " ") + ("E" if mu in MUS_EXT else " ")
                print("   %7.4f%% %s | %+.3f [%+.3f, %+.3f] | 뺀 %d · 소멸 %d"
                      % (100 * mu, lad, np.nanmean(D[mu]), lo, hi, drop[mu], extl[mu]))
                res["D|%g|%s|%s|%g" % (rho, geom, rule, mu)] = [float(np.nanmean(D[mu])), lo, hi, drop[mu], extl[mu]]
                if drop[mu] > MAX_DROP:
                    print("   ** 이 칸 판정 불가(뺀 배경 > %d)" % MAX_DROP)
                    UNDEC.add((rho, geom, rule, mu))

    def mean_curve(D, mus):
        return [float(np.nanmean(D[mu])) for mu in mus]

    def point_cross(D, mus):
        try:
            return cens_flip(mean_curve(D, mus), mus)
        except Undecidable:
            return float("nan")

    def pform(geom, label):
        """P1 꼴: ρ8 · roll · 확장 사다리 — Δ(0) 하한 > 0 ∧ Δ(1.625%) 상한 < 0"""
        r0 = res["D|8|%s|roll|0" % geom]
        r1 = res["D|8|%s|roll|%g" % (geom, MU_TOP)]
        und = any((RHO_SCARCE, geom, "roll", m) in UNDEC for m in (0.0, MU_TOP)) or \
            not (np.isfinite(r0[1]) and np.isfinite(r1[2]))
        ok = (not und) and r0[1] > 0 and r1[2] < 0
        D = T[(RHO_SCARCE, geom, "roll")]
        c = point_cross(D, MUS_EXT)
        cb = curve_boot(D, boot_idx("%s|cross" % label, n), MUS_EXT)
        mark = "판정 불가" if und else ("✅" if ok else "❌")
        lo_ok = "> 0" if r0[1] > 0 else "≤ 0"
        hi_ok = "< 0" if r1[2] < 0 else "≥ 0"
        print("\n-- %s %s · roll · ρ8 · 확장 사다리: Δ(0) 하한 %+.3f (%s) · Δ(%.4g%%) 상한 %+.3f (%s) → %s"
              % (label, geom, r0[1], lo_ok, 100 * MU_TOP, r1[2], hi_ok, mark))
        if not und and not ok and label == "P1":
            print("   갈래: %s" % ("(가) μ 0 우위(Δ(0) 하한 > 0)부터 서지 않는다" if r0[1] <= 0
                                    else "(나) 우위는 서지만 확장 사다리(접시 척도 1%) 안에서 반전이 구간상 서지 않는다"))
        ceq = c / SCALE if (isinstance(c, float) and math.isfinite(c) and c > 0) else c
        print("   교차(서술): %s (접시 글자변경 척도 %s) · 재추출 [%s, %s] · inf 비율 %.3f · 0 비율 %.3f · NaN 비율 %.3f"
              % (fmt_c(c), fmt_c(ceq, MU_TOP / SCALE), fmt_c(qcens(cb, 0.025)), fmt_c(qcens(cb, 0.975)),
                 float(np.mean(cb == math.inf)), float(np.mean(cb == 0)), float(np.mean(np.isnan(cb)))))
        return ("undecidable" if und else ok), c, r0, r1, float(np.mean(cb == math.inf))

    def kform(geom, label):
        """M1 — 열쇠는 셋 갈래: 반전(P1 과 같은 꼴) · 반전 없음(Δ(1.625%) 하한 > 0 ∨ 교차 재추출 inf 비율 ≥ 0.975) · 불확정"""
        v, c, r0, r1, pinf = pform(geom, label)
        if v == "undecidable":
            k = "undecidable"
        elif v:
            k = "reversal"
        elif r1[1] > 0 or pinf >= K_NONE_INF:
            k = "none"
        else:
            k = "indeterminate"
        print("   열쇠 분류(M1): %s · Δ(%.4g%%) 하한 %+.3f · inf 비율 %.3f(기준 ≥ %.3f)"
              % (KNAME[k], 100 * MU_TOP, r1[1], pinf, K_NONE_INF))
        return k, c

    # ── P1 ──
    p1, c_rf, r0, r1, _ = pform("loc", "P1")
    res["P1"] = p1
    res["cross_loc_roll_scarce"] = c_rf
    p1_und = p1 == "undecidable"
    # ── 해석 열쇠 K1 · K2 (결과 전 고정 · 같은 꼴) ──
    k1, c_g = kform("glob", "K1")
    k2, c_p = kform("place", "K2")
    res["K1_glob"], res["K2_place"] = k1, k2
    res["cross_glob_roll_scarce_ext"], res["cross_place_roll_scarce"] = c_g, c_p

    # ── G2b — 전역 대조가 해부26 을 재현하나 (원 사다리 · 비중단) ──
    Dg = T[(RHO_SCARCE, "glob", "roll")]
    c_go = point_cross(Dg, MUS_ORIG)
    pg = p4_stats(data, BG, RHO_SCARCE, "glob")
    g2b_lines = []
    if r26 is None or pg["k"] < MIN_BG_P4:
        g2b = "판정 불가"
    else:
        lo26, hi26 = r26["P4"]["lnS"][1], r26["P4"]["lnS"][2]
        idx = boot_idx("G2b", pg["k"])
        lo_g, hi_g = ci(pg["lns"][idx].mean(axis=1))
        overlap = not (hi_g < lo26 or lo_g > hi26)
        noflip = c_go == math.inf
        g2b = "✅" if (overlap and noflip) else "🟡"
        g2b_lines.append("   원 사다리 교차 %s (해부26: +inf) · (R−1)L 비 %.3f [%.3f, %.3f] 대 해부26 [%.3f, %.3f] → 겹침 %s"
                         % (fmt_c(c_go, 0.01), math.exp(pg["lns"].mean()), math.exp(lo_g), math.exp(hi_g),
                            math.exp(lo26), math.exp(hi26), "예" if overlap else "아니오"))
        g2b_lines.append("   Δ(0) glob %+.3f 대 해부26 %+.3f [%+.3f, %+.3f] (서술)"
                         % (res["D|8|glob|roll|0"][0], r26["D|8|roll|0"][0], r26["D|8|roll|0"][1], r26["D|8|roll|0"][2]))
    print("\n-- G2b 전역 대조 = 해부26 재현(새 시드 대역 · 비중단): %s" % g2b)
    for l in g2b_lines:
        print(l)
    if g2b == "🟡":
        print("   🟡 → 국소 대 전역 비교는 이 실험 안의 glob 팔과만 한다(해부26 수치와 섞지 않는다) · 시드 대역 차를 사람 확인")
    res["G2b"] = g2b

    # ── P3 ──
    def compare(name, D_a, D_b, paired, mus):
        idx_a = boot_idx("%s|a" % name, n)
        idx_b = idx_a if paired else boot_idx("%s|b" % name, n)
        ca = curve_boot(D_a, idx_a, mus)
        cbb = curve_boot(D_b, idx_b, mus)
        valid = np.isfinite(ca) & (ca > 0)
        pushed = valid & ((cbb == math.inf) | (cbb >= CROSS_FACTOR * ca))
        kept = valid & np.isfinite(cbb) & (cbb < CROSS_FACTOR * ca)
        p_push, p_keep = float(pushed.mean()), float(kept.mean())
        c_pt = point_cross(D_b, mus)
        verdict = "✅ 사라지거나 크게 밀린다" if p_push >= 0.975 else ("❌ 교차가 남는다(두 배 안)" if p_keep >= 0.975 else "🟡 구별 안 됨")
        return c_pt, p_push, p_keep, verdict, float(np.mean(cbb == math.inf)), float(valid.mean()), \
            int(np.isnan(ca).sum() + np.isnan(cbb).sum())

    D_lr, D_lf = T[(RHO_SCARCE, "loc", "roll")], T[(RHO_SCARCE, "loc", "find")]
    c_pt, pp, pk, v, pinf, pval, nnan = compare("P3", D_lr, D_lf, True, MUS_EXT)
    und = [(r, g, u, m) for (r, g, u, m) in UNDEC if r == RHO_SCARCE and g == "loc" and m in MUS_EXT]
    if und or p1_und or (isinstance(c_rf, float) and math.isnan(c_rf)) or (isinstance(c_pt, float) and math.isnan(c_pt)):
        v = "판정 불가(뺀 배경 · P1 판정 불가 · NaN) — 서술만: " + v
    elif not p1:
        v = "판정 보류(P1 ❌) — 서술만: " + v
    print("\n-- P3 loc · find-first · 희소 (배경 짝 재추출 · 확장 사다리): 교차 %s 대 loc-roll %s · P(밀림) %.3f · P(남음) %.3f · inf 비율 %.3f → %s"
          % (fmt_c(c_pt), fmt_c(c_rf), pp, pk, pinf, v))
    print("   P(c_rf 유효) %.3f · NaN 재추출 행 %d (분모 %d)" % (pval, nnan, N_BOOT))
    res["P3"] = {"cross": c_pt, "p_push": pp, "p_keep": pk, "verdict": v, "p_crf_valid": pval}

    # ── P4 ──
    pl = p4_stats(data, BG, RHO_SCARCE, "loc")
    print("\n-- P4 roll · 희소 · μ 0 · 끼운 팔 안 (자격 배경 loc %d · glob %d / %d · 최소 %d)" % (pl["k"], pg["k"], n, MIN_BG_P4))
    if pl["k"] < MIN_BG_P4 or pg["k"] < MIN_BG_P4:
        p4_txt = "판정 불가(자격 배경 부족)"
        res["P4"] = {"verdict": "undecidable", "k_loc": pl["k"], "k_glob": pg["k"]}
        print("   → %s" % p4_txt)
    else:
        ia, ib = boot_idx("P4|loc", pl["k"]), boot_idx("P4|glob", pg["k"])
        dd = pl["lnd"][ia].mean(axis=1) - pg["lnd"][ib].mean(axis=1)
        lo_a, hi_a = ci(dd)
        lo_b, hi_b = ci(pl["lns"][ia].mean(axis=1))
        p4a = lo_a > 0
        if S_BAND[0] <= lo_b and hi_b <= S_BAND[1]:
            p4b = "✅"
        elif hi_b < S_BAND[0] or lo_b > S_BAND[1]:
            p4b = "❌"
        else:
            p4b = "🟡"
        for nm, st in (("loc", pl), ("glob", pg)):
            print("   %-4s R 중앙 L2 %.2f · L4 %.2f · ln R(2) − ln R(4) %+.3f · (R−1)L 비 %.3f · T/R 중앙 %.3f · %.3f"
                  % (nm, np.median(st["R2"]), np.median(st["R4"]), st["lnd"].mean(), math.exp(st["lns"].mean()),
                     np.nanmedian(st["TR2"]), np.nanmedian(st["TR4"])))
        print("   P4a [ln R(2) − ln R(4)]_loc − [..]_glob = %+.3f [%+.3f, %+.3f] (독립 재추출) → %s"
              % (pl["lnd"].mean() - pg["lnd"].mean(), lo_a, hi_a, "✅" if p4a else "❌"))
        print("   P4b loc ln[S(2)/S(4)] = %+.3f [%+.3f, %+.3f] · 비 %.3f [%.3f, %.3f] · 띠 [0.8, 1.25] → %s"
              % (pl["lns"].mean(), lo_b, hi_b, math.exp(pl["lns"].mean()), math.exp(lo_b), math.exp(hi_b), p4b))
        p4 = p4a and p4b == "✅"
        p4_txt = "✅" if p4 else ("❌" if (not p4a or p4b == "❌") else "🟡")
        print("   P4 → %s" % p4_txt)
        res["P4"] = {"k_loc": pl["k"], "k_glob": pg["k"], "dlnR": [float(pl["lnd"].mean() - pg["lnd"].mean()), lo_a, hi_a],
                     "lnS_loc": [float(pl["lns"].mean()), lo_b, hi_b], "p4a": p4a, "p4b": p4b, "verdict": p4_txt}
    return finish(res, p1, p1_und, k1, k2, p4_txt, json_out, data, BG, T)


KNAME = {"reversal": "반전", "none": "반전 없음", "indeterminate": "불확정", "undecidable": "판정 불가"}


def interp(p1, k1, k2, r0_lo):
    """§5 해석 갈래(결과 전 고정 · v1.1 검토 M1 · M2 · S1 · S2 반영) — 기계 적용. k1 · k2 ∈ KNAME 의 키"""
    if p1 == "undecidable":
        return "판정 불가 — 해석하지 않는다"
    if p1 is True:
        if k1 == "undecidable" or k2 == "undecidable":
            return "P1 ✅ · 열쇠 판정 불가 — '국소 결핍이 필요한 성분' 을 쓰지 않는다"
        if k1 == "reversal":
            return "갈래 B: 확장 사다리에서는 전역에서도 반전 — 해부26 ❌ 는 μ 척도(오류 분할) 문제였고, 국소성이 필요하다는 근거는 이 실험에서 안 나온다"
        if k2 == "reversal":
            return "갈래 C: 국소 배치만(재료 전역)으로도 반전 — 반전은 공간 구조(배치 · 국소 경쟁)에 기대고, '재료' 국소성을 필요 성분으로 쓰지 않는다"
        if k1 == "none" and k2 == "none":
            return ("갈래 A: 국소 배치 위에서 재료 국소성이 (이 모형에서) 반전에 필요한 성분 — 전역 · 배치만 대조는 확장 사다리 끝 Δ 하한 > 0"
                    "(또는 교차 재추출 inf ≥ 0.975)으로 반전이 없다. 단서: loc 은 glob 대비 (a) 재료 국소성 (b) 국소 밀도 의존 출생 제한"
                    " (c) 혈연 군집 (d) 불임 자식의 이웃 칸 · 글자 점유를 함께 바꾸고, place 가 (b)~(d) 를 공유하므로 loc 대 place 가 (a) 를 가른다."
                    " 재료만 국소 · 배치 전역 조건은 시험하지 않았다")
        return "P1 ✅ · 열쇠 불확정: 국소에서만 구간상 반전 · 대조는 불확정 — '필요한 성분' 은 쓰지 않는다"
    if k1 == "reversal":
        return "갈래 E: 전역(glob)은 확장 사다리 안에서 반전하는데 국소(loc)는 아니다 — 국소 결핍은 반전에 필요하지 않고 오히려 막는다(이 모형에서)"
    if k1 == "undecidable":
        return "P1 ❌ · K1 판정 불가: E/D 구분 불가 — 전역이 반전하는지 알 수 없어 '국소성이 반전을 막는다' 도 '단순 모형으로 재현 안 됨' 도 쓰지 않는다"
    if r0_lo is not None and r0_lo <= 0:
        return "갈래 D(가): 국소 결핍을 더하면 μ 0 우위부터 안 선다 — 단순 모형으로 재현 안 됨"
    return ("갈래 D(나): 단순 모형(두 가정 + 불임 + 국소 결핍)에서 확장 사다리 안 반전이 구간상 서지 않는다(점추정 교차는 서술) — "
            "접시의 다른 성분(명령 실행 · 코드 구성 · 수명 분포 · 확산 · 오류 분할 등)을 봐야 한다")


def finish(res, p1, p1_und, k1, k2, p4_txt, json_out, data, BG, T):
    print("\n-- 서술: 끼운 팔 안 계통별 R 중앙 · I(실현 오류/위치 ÷ μ) 중앙")
    for rho in (RHO_SCARCE, RHO_RICH):
        for geom in GEOMS:
            rules = ["roll"] if geom == "place" else RULES
            for rule in rules:
                print("   ρ %g · %s · %s" % (rho, geom, rule))
                mus = MUS_EXT if geom == "place" else MUS_ALL
                for mu in mus:
                    a2, _, _, e2 = R_table(data, BG, rho, geom, rule, mu, "1_2")
                    a4, _, _, e4 = R_table(data, BG, rho, geom, rule, mu, "0_4")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        I2 = np.nanmedian(e2) / mu if mu > 0 else float("nan")
                        I4 = np.nanmedian(e4) / mu if mu > 0 else float("nan")
                        print("     %7.4f%% | R L2 %6.2f · L4 %6.2f | I L2 %6.2f · L4 %6.2f"
                              % (100 * mu, np.nanmedian(a2), np.nanmedian(a4), I2, I4))
    print("\n-- 서술: 배경(성장 끝 · 상주만) R 중앙 · 개체 수 — 규칙: 희소 R ≥ 2 · 풍부 R ≤ 1.1 (기하마다)")
    man = {}
    for rho in (RHO_SCARCE, RHO_RICH):
        for geom in GEOMS:
            Rs, ne = [], []
            for s in BG[rho]:
                g = data[rho][s]["by_geom"][geom]["grow"]
                c = g["counters"]["0_4"]
                Rs.append(c["draws"] / c["positions"])
                ne.append(g["n_end"])
            man[(rho, geom)] = float(np.median(Rs))
            print("   ρ %g · %-5s: R 중앙 %.2f [%.2f–%.2f] · 개체 수 중앙 %d" % (rho, geom, man[(rho, geom)], min(Rs), max(Rs), int(np.median(ne))))
    manip = all(man[(RHO_SCARCE, g)] >= 2 and man[(RHO_RICH, g)] <= 1.1 for g in GEOMS)
    print("   → %s" % ("조작 확인 ✅" if manip else "조작 확인 ❌ — 풍부(ρ32) 서술을 '재료 풍부' 로 읽지 않는다"))
    res["manip"] = {"%g|%s" % k: v for k, v in man.items()}
    res["manip_ok"] = manip
    print("\n-- 서술(감도): 가임만 센 Δ (확장 사다리)")
    for rho in (RHO_SCARCE, RHO_RICH):
        for geom, rule in (("loc", "roll"), ("loc", "find"), ("glob", "roll"), ("place", "roll")):
            Df, _, _ = delta_table(data, BG, rho, geom, rule, MUS_EXT, s_fert)
            row = " ".join("%+.3f" % float(np.nanmean(Df[mu])) for mu in MUS_EXT)
            try:
                cf = cens_flip([float(np.nanmean(Df[mu])) for mu in MUS_EXT], MUS_EXT)
            except Undecidable:
                cf = float("nan")
            print("   ρ %g · %-5s · %s | %s | 교차 %s" % (rho, geom, rule, row, fmt_c(cf)))
            res["S3|%g|%s|%s" % (rho, geom, rule)] = [float(np.nanmean(Df[mu])) for mu in MUS_EXT]
    print("\n-- 서술: 풍부(ρ32) 교차 점추정(확장 사다리) — 판정 아님")
    for geom, rule in (("loc", "roll"), ("loc", "find"), ("glob", "roll"), ("place", "roll")):
        D = T[(RHO_RICH, geom, rule)]
        try:
            c = cens_flip([float(np.nanmean(D[mu])) for mu in MUS_EXT], MUS_EXT)
        except Undecidable:
            c = float("nan")
        print("   ρ 32 · %-5s · %s: %s" % (geom, rule, fmt_c(c)))

    def mk(v):
        if v in KNAME:
            return KNAME[v]
        return v if isinstance(v, str) else ("✅" if v else "❌")
    r0_lo = res.get("D|8|loc|roll|0", [None, None])[1]
    it = interp(p1, k1, k2, r0_lo)
    print("\n해석 열쇠: K1 glob %s · K2 place %s" % (mk(k1), mk(k2)))
    if res["G2b"] == "🟡":
        it += " · 전역 대조가 해부26 수치를 새 시드 대역에서 재현하지 못했다(사람 확인)"
    print("해석(§5 기계 적용): %s" % it)
    print("단서(§5 고정): 오류 분할은 고치지 않았다 — 확장 사다리는 글자 변경 확률만 맞추고 오류 종류 몫(접시 바꾸기 3/5 · mini 1/3)은 다르다. 결과가 이 차이에 기댈 수 있다")
    print("결과 문장: 국소 결핍 모형에서 P1 %s · G2b %s · P3 %s · P4 %s"
          % ("판정 불가" if p1_und else mk(p1), res["G2b"], res["P3"]["verdict"], p4_txt))
    res["interp"] = it
    out = json_out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "_RESULT_mini27.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda x: None if (isinstance(x, float) and not math.isfinite(x)) else str(x))
    return 0


def arg(name):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else None


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    mode = sys.argv[1]
    if mode == "pilot":
        return pilot(sys.argv[2])
    if mode == "ident":
        bad, n = ident_check(sys.argv[2], sys.argv[3])
        print("== G2a mini27 glob 대 mini26 (시드 %d · 원 사다리): 바이트 같은 팔 %d · 문제 %d" % (IDENT_SEED, n, len(bad)))
        for b in bad:
            print("   ", b)
        return 0 if not bad else 1
    if mode == "judge":
        return judge(sys.argv[2], arg("--repro"), arg("--ident"), arg("--mini26"), arg("--json"))
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
