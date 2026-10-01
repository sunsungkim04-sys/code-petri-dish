#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 25 판정 — 핵심 결과 둘(짝지은 고리 효과 · 교차 μ 사다리)을 둘째 밀도에서 반복. 결과 전에 썼다.
사전등록: 해부25-사전등록-2026-10-01.md §4 (판정선 🔒) · 정의는 해부 13 · 15 E · 20(짝) · 해부 6C · 21(사다리)과 같다.

두 가지 모드
  qual   python3 pairs25_analyze.py qual <성장 폴더> <밀도 …> [--seeds a b …]
         성장 기록(replay.js)만 읽어 '멸종 없음 · 끝 우세 racld' 자격을 센다. Δ · s 를 계산하지도 인쇄하지도 않는다.
         밀도 둘을 주면 사전등록 §2 의 밀도 규칙(16 의 자격 비율 ≥ 0.6 → 16 · 아니면 12 의 비율 ≥ 0.6 → 12 · 둘 다 아니면 멈춤)을 적용해 인쇄한다.
  판정   python3 pairs25_analyze.py [petri 폴더]
         환경: PETRI_ROOT(inv25) · PETRI_GROW(main25) · PETRI_PILOT=1 이면 파일럿(PETRI_DENS · PETRI_SEEDS 허용 · 판정 아님)

  s  = ln((m_T + 0.5)/m_0) − ln((n_T − m_T + 0.5)/(n_0 − m_0))            (짝 · 해부 13 과 같다)
  ΔΔ = s(고리 L+1 코드) − s(고리 L 코드) · 배경 짝 · 짝 = 이웃 고리 등급의 모든 (짧은, 긴) `n` 삽입 쌍
  Δ  = s(팔) − 기준(WT 팔 11 개 평균) · s 는 ln((m_T+0.5)/m_0) − ln((r_T+0.5)/r_0)   (사다리 · 해부 6C · 21 과 같다)
  교차 = 씨앗 rascld 평균 Δ 곡선의 첫 부호 변화(+ → ≤ 0) 선형 보간 · 이 실험은 **검열**을 더한다:
         Δ(0) ≤ 0 이면 교차 = 0 · 끝까지 + 면 교차 = +∞(범위 위) — 방향 판정이 '정의 안 됨' 으로 빠지지 않게
  재추출: 배경 재추출 2,000 · 95% · 난수는 ('haebu25', 칸 이름)으로 키

  🔒 판정선(사전등록 §4)
    G0  파일 · 설정 · 배경 체크섬(dms/invbud 재성장 = 성장 기록) MATCH · 충실도 · 재료 보존 · 팔 목록 · 같은 시드 파일 사이 grow_checksum 같음 ·
        배경 목록 = 자격 규칙의 앞 20 · 기준 자료 자기 검산(정본 _RESULT_ 와 1e-9) — 하나라도 어긋나면 판정하지 않는다
    P1  주 밑 코드 다섯(racld rascld rsacld racldx rascled · 짝 54)의 합친 틱당 ΔΔ(배경마다 54 짝 평균 → 재추출) 상한 < 0
    P2  54 짝 중 상한 < 0 인 짝 수 ≥ 49 (90%)
    L1  씨앗 Δ(0) 하한 > 0 · Δ(1%) 상한 < 0   (해부 21 R1 과 같다)
    L2  교차(밀도 D) − 교차(밀도 8 기준 60 배경: 해부 6C 20 + 해부 21 40) · 두 묶음 독립 재추출 짝 차 d
        P(d > 0) ≥ 0.975 → ✅ 예측대로 위로 · P(d < 0) ≥ 0.975 → ❌ 반대(아래로) · 그 밖 → 🟡 구별 안 됨
        분모(검토 개정 R6 · 명시): P(d>0) · P(d<0) 의 분모는 언제나 N_BOOT = 2,000. 동률(d = 0 · 예: 둘 다 검열 0)과 ∞ − ∞(둘 다 범위 위)는
        어느 쪽에도 세지 않는다(분모에서 빼지 않는다 — 해부 21 의 '유효 수 분모' 와 달리 🟡 쪽으로 보수적). 뺀 수 · 동률 수를 인쇄한다.
        🟡 일 때 결과 문장은 "구별되지 않는다(차 d 의 95% 구간 [a, b]%p)" — 구간을 반드시 함께 쓰고 '같은 자리' 로 쓰지 않는다(R2).
        서술(R2): d 상한 < +0.037%p(= §3 띠 하단 0.21% − 0.173%)이면 "§3 수치 예측 하단이 배제된다".
    P4(부 예측) 합친 틱당 ΔΔ(D) − 합친 틱당 ΔΔ(밀도 8 기준 40 배경: 해부 15 E 20 + 해부 20 P 20 · 같은 다섯 밑 코드)
        하한 > 0 → 줄었다(예측) · 상한 < 0 → 커졌다 · 그 밖 구별 안 됨
    서술  acld · 해부 20 후보 셋의 W1 · N1 · N3 · 무너진 팔 · N2 · N3(주 다섯) · 사다리 괄호 · 기울기 · 곡률 · 순도 · 조작 확인(R)
          무너진 팔(R11): 해부 20(pairs20_analyze.py:156)과 글자 그대로 — **mut 팔만**, 끝 m < 0.10 × 처음 m 인 팔 수의 합(배경 수도 함께 인쇄)
          조작 확인 규칙(R4 · 서술 규칙 · 파일럿 뒤 · 본 데이터 전에 고정): 배경마다 racld 의 n 삽입 여섯 코드의 ln R 평균 →
          밀도 D 와 밀도 8(해부 15 Eib 20) 각각 배경 재추출 2,000(독립) → 차 d_R(D − 8)의 95% 구간으로 셋:
            상한 < ln 0.5  → "§3 전제만큼 줄었다"(§3 크기 근거의 비 0.27~0.5 중 큰 쪽 · 결과 전 숫자) — L2 를 '부족이 덜한 세계' 로 읽어도 된다
            상한 < 0       → "줄었지만 §3 전제(비 ≤ 0.5)에 못 미친다"(비 e^d_R 를 함께 적는다) — L2 를 '부족이 덜한 세계' 로 읽지 않는다
            그 밖          → "줄지 않았다" — 사전등록 §5 마지막 줄: L2 를 '부족이 덜한 세계' 로 읽지 않는다
    G0 추가(R1): 주 다섯의 54 짝이 배경 20 마다 전부 있어야 한다 — P 를 계산하기 전에 세고, 하나라도 빠지면 멈춘다.
"""
import glob
import json
import math
import os
import random
import sys

# ---- 판정선 상수 (사전등록 §4 · 결과 전 고정) ----
D_FIXED = "16"                 # §2a 밀도 규칙 적용(10-01 파일럿 자격 16: 9/12 · 12: 8/12 → 16) — Δ 보기 전
PRIMARY = ["racld", "rascld", "rsacld", "racldx", "rascled"]
DESC = ["acld", "rcald", "crald", "reasccld"]
PBASES = PRIMARY + DESC
N_PRIMARY_PAIRS = 54
P2_MIN = 49
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
LCODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
P_DIR = 0.975
TOL_NEG, PT_LO, PT_HI = 0.15, -1.0, -0.2       # 해부 13 그대로(서술용)
COLLAPSE = 0.10                                 # 해부 20 §3 무너짐 정의(서술용)
NBG = 20
CANDS_MAIN = list(range(2501, 2551))
PILOT_SEEDS = set(range(2591, 2603))
PRED_BAND = (0.0021, 0.0065)                    # §3 수치 예측(서술 · 판정 아님)
DENS_RULE_FRAC = 0.6


def rng(name):
    return random.Random("haebu25|%s" % name)


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def copies_per_loop(code):
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    return code[start:li].count("c")


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def pairs_of(codes):
    out = []
    Ls = sorted(set(codes.values()))
    for L1, L2 in zip(Ls, Ls[1:]):
        if L2 - L1 != 1: continue
        for a in codes:
            for b in codes:
                if codes[a] == L1 and codes[b] == L2: out.append((a, b))
    return out


def same_loop_pairs(codes):
    cl = sorted(codes); out = []
    for i, a in enumerate(cl):
        for b in cl[i + 1:]:
            if codes[a] == codes[b]: out.append((a, b))
    return out


def s_pair(arm):        # 짝 — 해부 13 · 20 과 같다
    r = arm["series"]; n0, m0 = r[0][1], r[0][2]; nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


s_ladder = s_pair       # 사다리 — 해부 6C · 21 의 s_of 와 같은 식(r = n − m)


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    if xs[f] == xs[c]:
        return xs[f]
    if math.isinf(xs[f]) or math.isinf(xs[c]):
        return xs[c] if (k - f) >= 0.5 else xs[f]
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot_mean(v, name):
    r = rng(name)
    m = sum(v) / len(v)
    bs = [sum(v[r.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975), bs


def fmt(m, lo, hi):
    return "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def fmt_mu(x):
    if x is None: return "-"
    if math.isinf(x): return "> 1% (∞)"
    return "%.3f%%" % (100 * x)


# =====================================================================  qual
def qual_mode(argv):
    if len(argv) < 2:
        print("사용: pairs25_analyze.py qual <성장 폴더> <밀도 …> [--seeds a b …]"); sys.exit(2)
    root = argv[0]
    rest = argv[1:]
    seeds = CANDS_MAIN
    if "--seeds" in rest:
        i = rest.index("--seeds"); seeds = [int(x) for x in rest[i + 1:]]; rest = rest[:i]
    dens = rest
    out = {"root": root, "seeds": seeds, "dens": {}}
    print("== 해부 25 자격(Δ 없음) — 성장 %s · 시드 %d 개(%d~%d)" % (root, len(seeds), min(seeds), max(seeds)))
    for d in dens:
        rows, miss = [], []
        for s in seeds:
            f = os.path.join(root, "mat_d%s_mu0p01_s%05d.json" % (d, s))
            if not os.path.exists(f):
                miss.append(s); continue
            z = json.load(open(f))
            fc = z["final_counts"]; pop = sum(c for _, c in fc)
            top = fc[0][0] if fc else None
            sh = {k: c / pop for k, c in fc[:5]} if pop else {}
            rows.append({"seed": s, "extinct": z["extinct_at"], "top": top, "pop": pop,
                         "racld_share": (dict(fc).get("racld", 0) / pop) if pop else 0.0,
                         "top_share": (fc[0][1] / pop) if pop else 0.0,
                         "top3": [[k, c] for k, c in fc[:3]], "cons_viol": z["conservation_violations"],
                         "ticks": z["ticks_run"], "planned": z["ticks_planned"], "mu": z["opts"]["mu"]})
        q = [r for r in rows if r["extinct"] < 0 and r["top"] == "racld"]
        ext = [r for r in rows if r["extinct"] >= 0]
        other = {}
        for r in rows:
            if r["extinct"] < 0 and r["top"] != "racld": other[r["top"]] = other.get(r["top"], 0) + 1
        frac = len(q) / len(seeds)
        med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
        print("\n-- 밀도 %s: 기록 %d/%d · 멸종 %d · 끝 우세 racld %d → 자격 비율 %.3f (%d/%d)"
              % (d, len(rows), len(seeds), len(ext), len(q), frac, len(q), len(seeds)))
        print("   racld 아닌 끝 우세: %s" % (" · ".join("%s %d" % kv for kv in sorted(other.items(), key=lambda kv: -kv[1])) or "없음"))
        print("   자격 배경의 racld 몫 중앙 %.3f · 개체 수 중앙 %s · 재료 보존 위반 %d" % (med([r["racld_share"] for r in q]), med([r["pop"] for r in rows if r["pop"]]), sum(r["cons_viol"] for r in rows)))
        if miss: print("   🚨 기록 없음 %s" % miss)
        for r in rows:
            print("     s%d  %s  pop %5d  top3 %s" % (r["seed"], ("멸종 %d" % r["extinct"]) if r["extinct"] >= 0 else ("자격" if r in q else "다름"),
                                                    r["pop"], " ".join("%s:%d" % (k, c) for k, c in r["top3"])))
        out["dens"][d] = {"n": len(seeds), "records": len(rows), "missing": miss, "extinct": len(ext), "qualified": [r["seed"] for r in q],
                          "frac": frac, "other_tops": other, "rows": rows}
    if set(dens) == {"12", "16"}:
        n = len(seeds); need = math.ceil(DENS_RULE_FRAC * n - 1e-9)
        f16, f12 = len(out["dens"]["16"]["qualified"]), len(out["dens"]["12"]["qualified"])
        pick = "16" if f16 >= need else ("12" if f12 >= need else None)
        print("\n== 밀도 규칙(사전등록 §2 · 결과 전 고정): 자격 ≥ %d/%d 이면 그 밀도 · 16 먼저" % (need, n))
        print("   밀도 16 %d/%d · 밀도 12 %d/%d → %s" % (f16, n, f12, n, ("밀도 %s" % pick) if pick else "🚨 둘 다 미달 — 멈춤 · 사람 결정"))
        out["rule"] = {"need": need, "q16": f16, "q12": f12, "pick": pick}
    fo = os.environ.get("PETRI_OUT")
    if fo:
        json.dump(out, open(fo, "w"), ensure_ascii=False, indent=1); print("기록 →", fo)


# =====================================================================  판정
BASE = None
bad, notes = [], []


def expected_seeds(GROW, D, PILOT):
    if not PILOT and os.environ.get("PETRI_SEEDS"):
        print("🚨 판정 모드에서 PETRI_SEEDS 를 쓸 수 없다 — 파일럿이면 PETRI_PILOT=1"); sys.exit(1)
    cands = [int(x) for x in os.environ["PETRI_SEEDS"].split()] if os.environ.get("PETRI_SEEDS") else CANDS_MAIN
    picked = []
    for s in cands:
        f = os.path.join(BASE, GROW, "mat_d%s_mu0p01_s%05d.json" % (D, s))
        if not os.path.exists(f):
            bad.append("성장 기록 없음 s%d" % s); continue
        z = json.load(open(f))
        if z["ticks_planned"] != 50000 or abs(z["opts"]["mu"] - 0.01) > 1e-12 or float(z["opts"]["density"]) != float(D) or z["cond"] != "mat":
            bad.append("성장 설정 s%d" % s)
        if z["conservation_violations"]:
            bad.append("성장 재료 보존 s%d" % s)
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
            picked.append(s)
    if os.environ.get("PETRI_SEEDS"):
        if len(picked) != len(cands): bad.append("파일럿 시드 중 자격 아닌 것이 있다")
        return sorted(picked)
    if len(picked) < NBG:
        bad.append("자격 배경 %d < %d — 판정 불가" % (len(picked), NBG))
    return sorted(picked[:NBG])


def check_dms(z, D, mu, wt, codes_expected, tag):
    s = z["seed"]
    if abs(z["mu_assay"] - mu) > 1e-12 or z["wt"] != wt or z["arms_mode"] != "list" or z["cond"] != "mat": bad.append("설정 %s s%d" % (tag, s))
    if float(z.get("density", 8)) != float(D): bad.append("밀도 %s s%d" % (tag, s))
    if any(z.get(k) for k in ("find_first", "remember_die", "no_mat", "cosmic_off", "kids_on")) or z.get("grow_mu") is not None or z.get("age0") is not None or z.get("ancestor"):
        bad.append("규칙 · 세계 옵션 %s s%d" % (tag, s))
    if "exact_series" not in z: bad.append("exact 없음 %s s%d" % (tag, s))
    if not z["checksum_match"]: bad.append("배경 체크섬(재성장 ≠ 성장 기록) %s s%d" % (tag, s))
    if not z["fidelity"]["same"]: bad.append("충실도 %s s%d" % (tag, s))
    if z["grow_extinct"] >= 0 or not z["grow_top"] or z["grow_top"][0][0] != "racld": bad.append("배경 자격 %s s%d" % (tag, s))
    if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s s%d" % (tag, s))
    muts = [a for a in z["arms"] if a["kind"] == "mut"]
    keys = [a["key"] for a in muts]
    # 짝은 해부 20 과 같이 순서까지(목록 모드 팔 순서 = n 삽입 순서) · 사다리는 집합으로(dms.js 팔 순서는 돌연변이 목록 순서)
    same = (keys == list(codes_expected)) if tag.startswith("P_") else (sorted(keys) == sorted(codes_expected))
    if not same or len(z["arms"]) != 11 + len(codes_expected): bad.append("팔 목록 %s s%d" % (tag, s))


def load_pairs_dir(pattern, b):
    """{seed: {code: s}} · 무너진 팔 배경 집합 · grow_checksum"""
    d, coll, gcs = {}, {}, {}
    for f in sorted(glob.glob(pattern)):
        z = json.load(open(f)); s = z["seed"]
        muts = [a for a in z["arms"] if a["kind"] == "mut"]
        d[s] = {a["key"]: s_pair(a) for a in muts}
        coll[s] = sum(1 for a in muts if a["series"][-1][2] < COLLAPSE * a["series"][0][2])   # R11: 해부 20 과 같이 mut 팔만
        gcs[s] = z["grow_checksum"]
    return d, coll, gcs


def primary_pair_problems(S, SEEDS):
    """검토 개정 R1 — 주 다섯의 짝(54)이 배경마다 전부 있어야 한다. P 를 계산하기 전(G0)에 센다."""
    probs = []
    n = sum(len(pairs_of(inserts(b, "n"))) for b in PRIMARY)
    if n != N_PRIMARY_PAIRS:
        probs.append("주 짝 정의 수 %d ≠ %d" % (n, N_PRIMARY_PAIRS))
    for b in PRIMARY:
        for a, c in pairs_of(inserts(b, "n")):
            miss = [s for s in SEEDS if not (s in S.get(b, {}) and a in S[b][s] and c in S[b][s])]
            if miss:
                probs.append("주 짝 %s|%s|%s 배경 빠짐 %d (예 %s)" % (b, a, c, len(miss), miss[:3]))
    if not SEEDS:
        probs.append("배경 0")
    return probs


def halt(OUT):
    print("🚨 점검 실패 — 판정하지 않는다")
    OUT["halted"] = bad[:50]
    json.dump(clean(OUT), open(os.path.join(BASE, "_RESULT_pairs25_halted.json"), "w"), ensure_ascii=False, indent=1)   # 판정 파일 이름을 비워 둔다
    sys.exit(1)


def per_tick_by_bg(S, bases):
    """배경마다 주어진 밑 코드들의 모든 짝 ΔΔ/ΔL 평균"""
    per = {}
    for b in bases:
        codes = inserts(b, "n")
        for a, c in pairs_of(codes):
            dL = codes[c] - codes[a]
            for s in S[b]:
                if a in S[b][s] and c in S[b][s]:
                    per.setdefault(s, []).append((S[b][s][c] - S[b][s][a]) / dL)
    return {s: sum(v) / len(v) for s, v in per.items()}, {s: len(v) for s, v in per.items()}


def ladder_deltas(files_by_mu_seed, bgs):
    delta = {c: {mu: [] for mu in MUS} for c in LCODES}
    pure = {c: {mu: [] for mu in MUS} for c in LCODES + ["ref"]}
    for mu in MUS:
        for s in bgs:
            z = files_by_mu_seed[(mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_ladder(a) for a in wt) / len(wt)
            for a in z["arms"]:
                last = a["series"][-1]
                share = last[3] / last[2] if last[2] else float("nan")
                if a["kind"] == "mut":
                    delta[a["key"]][mu].append(s_ladder(a) - base); pure[a["key"]][mu].append(share)
                elif a["kind"] == "ref":
                    pure["ref"][mu].append(share)
    return delta, pure


def flip_plain(curve):          # 해부 21 과 같다(자기 검산용)
    for (a, pa), (b, pb) in zip(zip(MUS, curve), zip(MUS[1:], curve[1:])):
        if pa > 0 >= pb:
            return a + (b - a) * pa / (pa - pb)
    return None


def flip_cens(curve):           # 이 실험의 검열 교차(§4 L2)
    if curve[0] <= 0:
        return 0.0
    f = flip_plain(curve)
    return math.inf if f is None else f


def fit(xs, ys, deg):
    n = deg + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]; b[i], b[p] = b[p], b[i]
        for r in range(n):
            if r != i:
                f = A[r][i] / A[i][i]
                A[r] = [a - f * c for a, c in zip(A[r], A[i])]; b[r] -= f * b[i]
    return [b[i] / A[i][i] for i in range(n)]


class LSet:
    def __init__(self, name, delta):
        self.name, self.delta = name, delta
        self.n = len(delta["rascld"][0.0])
        r = rng("ladder|" + name)
        self.idx = [[r.randrange(self.n) for _ in range(self.n)] for _ in range(N_BOOT)]

    def mean(self, c, mu, ix=None):
        v = self.delta[c][mu]; ix = range(self.n) if ix is None else ix
        return sum(v[i] for i in ix) / len(ix)

    def bmean(self, c, mu):
        bs = [self.mean(c, mu, ix) for ix in self.idx]
        return self.mean(c, mu), pct(bs, 0.025), pct(bs, 0.975)

    def flips(self):
        pt = flip_cens([self.mean("rascld", mu) for mu in MUS])
        return pt, [flip_cens([self.mean("rascld", mu, ix) for mu in MUS]) for ix in self.idx]

    def curve(self, c, deg):
        pt = fit(MUS, [self.mean(c, mu) for mu in MUS], deg)
        pts = [fit(MUS, [self.mean(c, mu, ix) for mu in MUS], deg) for ix in self.idx]
        return [(pt[k], pct([p[k] for p in pts], 0.025), pct([p[k] for p in pts], 0.975)) for k in range(deg + 1)], pts


def ib_stats(pattern, codes):
    """invbud 팔 → R(시도 / 성공) — 해부 15 E3 와 같은 셈"""
    KEYS = ["copy_ok", "copy_stall", "copy_del", "copy_phase"]
    IB = {c: {} for c in codes}; gcs = {}; nfile = 0
    for f in sorted(glob.glob(pattern)):
        z = json.load(open(f)); s = z["seed"]; nfile += 1; gcs[s] = z["grow_checksum"]
        if list(z["codes"]) != list(codes): bad.append("IB 코드 목록 s%d" % s)
        for a in z["arms"]:
            if a["kind"] != "mut": continue
            d = {k: sum(bn[k][1] for bn in a["bins"]) for k in KEYS}
            if a["stopped"] >= 0 and a["series"][-1][2] == 0: continue
            if d["copy_ok"] == 0: continue
            IB[a["key"]][s] = {"R": (d["copy_ok"] + d["copy_stall"] + d["copy_del"]) / d["copy_ok"], **d}
    return IB, gcs, nfile


def ib_summary(IB, codes, name):
    pooledR = {}
    for c in codes:
        ds = IB[c].values(); ok_ = sum(d["copy_ok"] for d in ds); st = sum(d["copy_stall"] + d["copy_del"] for d in ds)
        if ok_: pooledR[c] = (ok_ + st) / ok_
    Rs = sorted(pooledR.values())
    medR = Rs[len(Rs) // 2] if Rs else float("nan")
    c4 = [c for c in codes if codes[c] == 4]; c5 = [c for c in codes if codes[c] == 5]
    seeds = sorted(set.intersection(*[set(IB[c]) for c in codes])) if all(IB[c] for c in codes) else []
    v = [sum(math.log(IB[c][s]["R"]) for c in c4) / len(c4) - sum(math.log(IB[c][s]["R"]) for c in c5) / len(c5) for s in seeds]
    a1 = boot_mean(v, "ibA1|" + name)[:3] if len(v) >= 2 else (float("nan"),) * 3
    return medR, pooledR, a1, len(seeds)


def clean(x):
    """JSON 에 Infinity/NaN 을 쓰지 않는다(다른 도구가 조용히 못 읽는 함정) — 문자열로"""
    if isinstance(x, float) and (math.isinf(x) or math.isnan(x)):
        return "inf" if x > 0 else ("-inf" if math.isinf(x) else "nan")
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return x


def judge():
    global BASE
    BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
    PILOT = os.environ.get("PETRI_PILOT") == "1"
    ROOT = os.environ.get("PETRI_ROOT", "pilot25" if PILOT else "inv25")
    GROW = os.environ.get("PETRI_GROW", "pilot25g" if PILOT else "main25")
    D = os.environ.get("PETRI_DENS", D_FIXED) if PILOT else D_FIXED
    if not PILOT and (os.environ.get("PETRI_DENS") or os.environ.get("PETRI_ROOT", "inv25") != "inv25" or os.environ.get("PETRI_GROW", "main25") != "main25"):
        print("🚨 판정 모드에서 PETRI_DENS · PETRI_ROOT · PETRI_GROW 를 바꿀 수 없다"); sys.exit(1)
    if D not in ("12", "16"):
        print("🚨 밀도가 정해지지 않았다(D_FIXED=%s) — 판정하지 않는다" % D); sys.exit(1)
    OUT = {"prereg": "해부25-사전등록-2026-10-01.md", "root": ROOT, "grow": GROW, "density": int(D), "pilot": PILOT}
    SEEDS = expected_seeds(GROW, D, PILOT)
    if not PILOT and set(SEEDS) & PILOT_SEEDS: bad.append("본 배경이 파일럿 시드와 겹친다")
    OUT["seeds"] = SEEDS
    print("== 해부 25 %s— 밀도 %s · 루트 %s · 성장 %s · 배경 %d %s" % ("파일럿 " if PILOT else "판정 ", D, ROOT, GROW, len(SEEDS), SEEDS))

    # ---- 짝(P) 파일
    S, COLL, GC = {}, {}, {}
    for b in PBASES:
        codes = inserts(b, "n")
        pat = os.path.join(BASE, "%sP_%s" % (ROOT, b), "mu0", "dms_mat_%s_s*_list_d%s.json" % (b, D))
        for f in sorted(glob.glob(pat)):
            check_dms(json.load(open(f)), D, 0.0, b, codes, "P_" + b)
        S[b], COLL[b], g = load_pairs_dir(pat, b)
        if sorted(S[b]) != SEEDS: bad.append("P 배경 목록 %s (%d/%d)" % (b, len(S[b]), len(SEEDS)))
        for s, x in g.items(): GC.setdefault(s, set()).add(x)
    # ---- 사다리(L) 파일
    FL = {}
    for mu in MUS:
        for f in sorted(glob.glob(os.path.join(BASE, ROOT + "L", mtag(mu), "dms_mat_racld_s*_list_d%s.json" % D))):
            z = json.load(open(f)); check_dms(z, D, mu, "racld", LCODES, "L_%s" % mu)
            FL[(mu, z["seed"])] = z; GC.setdefault(z["seed"], set()).add(z["grow_checksum"])
    miss = [(mu, s) for mu in MUS for s in SEEDS if (mu, s) not in FL]
    extra = sorted({k[1] for k in FL} - set(SEEDS))
    if miss: bad.append("L 파일 빠짐 %d (예 %s)" % (len(miss), miss[:3]))
    if extra: bad.append("L 배경 밖 파일 %s" % extra[:5])
    # ---- 계수(IB) 파일
    rcodes = inserts("racld", "n")
    IBD, gib, nib = ib_stats(os.path.join(BASE, ROOT + "ib_racld", "mu0", "ib_inv_s*_mu0_d%s.json" % D), list(rcodes))
    for f in glob.glob(os.path.join(BASE, ROOT + "ib_racld", "mu0", "ib_inv_s*_mu0_d%s.json" % D)):
        z = json.load(open(f))
        if not z["checksum_match"] or float(z.get("density") or 8) != float(D) or z.get("grow_mu") is not None or z.get("ancestor"):
            bad.append("IB 설정 · 체크섬 s%d" % z["seed"])
    if sorted(gib) != SEEDS: bad.append("IB 배경 목록 (%d/%d)" % (len(gib), len(SEEDS)))
    for s, x in gib.items(): GC.setdefault(s, set()).add(x)
    ndiff = sum(1 for s in GC if len(GC[s]) > 1)
    print("  파일: P %d · L %d · IB %d · 같은 시드 파일 사이 grow_checksum 이 다른 배경 %d/%d"
          % (sum(len(S[b]) for b in PBASES), len(FL), nib, ndiff, len(GC)))
    if ndiff: bad.append("같은 시드 파일 사이 배경 불일치 %d" % ndiff)

    # ---- 밀도 8 기준 자료 · 자기 검산
    r21 = json.load(open(os.path.join(BASE, "_RESULT_inv21.json")))
    rinv = json.load(open(os.path.join(BASE, "_RESULT_inv.json")))
    F8 = {}
    for root in ("inv6", "inv21"):
        for mu in MUS:
            for f in glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json")):
                z = json.load(open(f)); F8[(root, mu, z["seed"])] = z
    old_bgs = sorted({s for (r, mu, s) in F8 if r == "inv6"})
    new_bgs = list(r21["backgrounds_used"])
    ok8 = len(old_bgs) == 20 and len(new_bgs) == 40 and all(("inv6", mu, s) in F8 for mu in MUS for s in old_bgs) and all(("inv21", mu, s) in F8 for mu in MUS for s in new_bgs)
    if not ok8:
        bad.append("밀도 8 사다리 기준(inv6 · inv21) 불완전")
    else:
        dO, _ = ladder_deltas({(mu, s): F8[("inv6", mu, s)] for mu in MUS for s in old_bgs}, old_bgs)
        dN, _ = ladder_deltas({(mu, s): F8[("inv21", mu, s)] for mu in MUS for s in new_bgs}, new_bgs)
        dev = 0.0
        for mu in MUS:
            dev = max(dev, abs(sum(dO["rascld"][mu]) / 20 - rinv["delta"]["rascld"][str(mu)][0]))
            for c in LCODES:
                dev = max(dev, abs(sum(dN[c][mu]) / 40 - r21["delta_new"][c][str(mu)][0]), abs(sum(dO[c][mu]) / 20 - r21["delta_old"][c][str(mu)][0]))
        d8 = {c: {mu: dO[c][mu] + dN[c][mu] for mu in MUS} for c in LCODES}
        fpool = flip_plain([sum(d8["rascld"][mu]) / 60 for mu in MUS])
        dev = max(dev, abs(fpool - r21["desc"]["flip_pooled"]))
        print("  기준 사다리(밀도 8 · 해부 6C 20 + 해부 21 40 = 60): Δ · 합친 교차가 정본(_RESULT_inv.json · _RESULT_inv21.json)과 최대 %.1e 차 · 합친 교차 %s"
              % (dev, fmt_mu(fpool)))
        if dev > 1e-9: bad.append("기준 사다리 자기 검산 %.2e" % dev)
    r15 = json.load(open(os.path.join(BASE, "_RESULT_pairs15.json")))
    r20 = json.load(open(os.path.join(BASE, "_RESULT_pairs20.json")))
    S8 = {}
    dev = 0.0
    for b in PRIMARY:
        codes = inserts(b, "n")
        e, _, _ = load_pairs_dir(os.path.join(BASE, "inv15E_%s" % b, "mu0", "dms_mat_%s_s*_list.json" % b), b)
        p, _, _ = load_pairs_dir(os.path.join(BASE, "inv20P_%s" % b, "mu0", "dms_mat_%s_s*_list.json" % b), b)
        if len(e) != 20 or len(p) != 20 or set(e) & set(p): bad.append("밀도 8 짝 기준 %s (15E %d · 20P %d)" % (b, len(e), len(p)))
        for a, c in pairs_of(codes):
            k = "E1|%s|%s|%s" % (b, a, c)
            if k in r15:
                dev = max(dev, abs(sum(e[s][c] - e[s][a] for s in e) / len(e) - r15[k]["dd"]))
            else:
                bad.append("정본 키 없음 %s" % k)
        w, _ = per_tick_by_bg({b: p}, [b])
        dev = max(dev, abs(sum(w.values()) / len(w) - r20["per_base"][b]["W1"]["m"]))
        S8[b] = dict(e); S8[b].update(p)
    print("  기준 짝(밀도 8 · 해부 15 E 20 + 해부 20 P 20 = 40 · 주 다섯): 짝 ΔΔ · W1 이 정본(_RESULT_pairs15.json · _RESULT_pairs20.json)과 최대 %.1e 차" % dev)
    if dev > 1e-9: bad.append("기준 짝 자기 검산 %.2e" % dev)
    IB8, _, nib8 = ib_stats(os.path.join(BASE, "inv15Eib_racld", "mu0", "ib_inv_s*_mu0.json"), list(rcodes))
    med8, _, a18, _ = ib_summary(IB8, rcodes, "d8")
    if abs(a18[0] - r15["E3A1|racld|4|5"]["d"]) > 1e-9: bad.append("기준 계수 자기 검산 %.2e" % abs(a18[0] - r15["E3A1|racld|4|5"]["d"]))
    print("  기준 계수(밀도 8 · 해부 15 Eib racld %d): ln R(4) − ln R(5) %+.3f = 정본 %+.3f" % (nib8, a18[0], r15["E3A1|racld|4|5"]["d"]))

    pp = primary_pair_problems(S, SEEDS)       # R1 — 판정 전에
    bad.extend(pp)
    print("  주 짝: 정의 %d · 배경마다 빠진 짝 문제 %d" % (sum(len(pairs_of(inserts(b, "n"))) for b in PRIMARY), len(pp)))
    print("  점검 문제 %d" % len(bad))
    for x in bad[:20]: print("    ", x)
    if bad and not PILOT:
        halt(OUT)
    if bad and PILOT and (not SEEDS or not FL):
        print("🚨 파일럿 자료 부족 — 멈춤"); sys.exit(1)

    # ================= P 짝
    print("\n== P — 짝지은 고리 효과 · 밀도 %s · μ 0" % D)
    res = {}
    pair_ok, pair_neg, npairs = 0, 0, 0
    for b in PBASES:
        codes = inserts(b, "n"); role = "주" if b in PRIMARY else "서술"
        k = copies_per_loop(b)
        KB = sum(1 for s in COLL[b] if COLL[b][s] > 0)
        print("\n-- `%s` (%s · 고리 %d · 고리당 복사 %d) · 무너진 팔 %d · 그런 배경 %d/%d(서술 · 해부 20 정의 · mut 팔)" % (b, role, loop_ticks(b), k, sum(COLL[b].values()), KB, len(COLL[b])))
        n1, n2, n3, negs = [], [], [], []
        for a, c in pairs_of(codes):
            dL = codes[c] - codes[a]
            v = [S[b][s][c] - S[b][s][a] for s in sorted(S[b]) if a in S[b][s] and c in S[b][s]]
            if len(v) < 2: continue
            m, lo, hi, _ = boot_mean(v, "pair|%s|%s|%s" % (b, a, c))
            per = m / dL; ok1 = hi < 0; ok2 = PT_LO <= per <= PT_HI
            n1.append(ok1); n2.append(ok2)
            if b in PRIMARY:
                npairs += 1; pair_ok += ok1; pair_neg += m < 0
            print("     ΔΔ %-10s(%d) − %-10s(%d) = %s · 배경 %d · 틱당 %+.3f %s" % (c, codes[c], a, codes[a], fmt(m, lo, hi), len(v), per, "✅" if ok1 else "❌ 상한≥0"))
            OUT["pair|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "per_tick": per, "n_bg": len(v), "N1": ok1, "N2": ok2}
        for a, c in same_loop_pairs(codes):
            v = [S[b][s][c] - S[b][s][a] for s in sorted(S[b]) if a in S[b][s] and c in S[b][s]]
            if len(v) < 2: continue
            m, lo, hi, _ = boot_mean(v, "neg|%s|%s|%s" % (b, a, c)); ok3 = abs(m) <= TOL_NEG; n3.append(ok3); negs.append(abs(m))
            OUT["neg|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "N3": ok3}
            if not ok3: print("     [음성] %-10s vs %-10s (%d틱) = %s ❌ 0.15 밖" % (a, c, codes[a], fmt(m, lo, hi)))
        w, _ = per_tick_by_bg({b: S[b]}, [b])
        if len(w) >= 2:
            wm, wlo, whi, _ = boot_mean([w[s] for s in sorted(w)], "W1|%s" % b)
        else:
            wm = wlo = whi = float("nan")
        print("   W1 합친 틱당 %s %s · N1 %d/%d · N2 %d/%d · N3 %d/%d(최대 |ΔΔ| %.3f)%s"
              % (fmt(wm, wlo, whi), "✅" if whi < 0 else "❌", sum(n1), len(n1), sum(n2), len(n2), sum(n3), len(n3), max(negs) if negs else float("nan"),
                 ("" if k == 1 or wm != wm else " · W4 글자당 %+.3f" % (wm * k))))
        res[b] = {"role": role, "W1": [wm, wlo, whi], "N1": [sum(n1), len(n1)], "N2": [sum(n2), len(n2)], "N3": [sum(n3), len(n3)],
                  "neg_max": max(negs) if negs else None, "K_B": KB, "collapsed_arms": sum(COLL[b].values())}
    OUT["per_base"] = res
    wD, nD = per_tick_by_bg(S, PRIMARY)
    p1m, p1lo, p1hi, bsD = boot_mean([wD[s] for s in sorted(wD)], "P1|D")
    P1 = p1hi < 0
    # R1 — 두 번째 그물: 센 짝 수 · 배경마다 짝 수가 정의와 다르면 판정 전에 멈춘다(위 G0 검사가 있으면 닿지 않는 길)
    if npairs != N_PRIMARY_PAIRS: bad.append("주 짝 수 %d ≠ %d" % (npairs, N_PRIMARY_PAIRS))
    if set(nD.values()) != {N_PRIMARY_PAIRS} or sorted(nD) != SEEDS: bad.append("배경마다 주 짝 수 %s · 배경 %d/%d" % (sorted(set(nD.values())), len(nD), len(SEEDS)))
    if bad and not PILOT:
        halt(OUT)
    P2 = pair_ok >= P2_MIN
    print("\n-- P1 주 다섯 합친 틱당 ΔΔ(배경 %d · 배경마다 짝 %s) %s → %s" % (len(wD), sorted(set(nD.values())), fmt(p1m, p1lo, p1hi), "✅" if P1 else "❌"))
    print("-- P2 주 짝 상한 < 0: %d/%d (선 ≥ %d) → %s · 서술: 점추정 음수 %d/%d" % (pair_ok, npairs, P2_MIN, "✅" if P2 else "❌", pair_neg, npairs))
    w8, _ = per_tick_by_bg(S8, PRIMARY)
    p8m, p8lo, p8hi, bs8 = boot_mean([w8[s] for s in sorted(w8)], "P1|d8")
    dd = sorted(a - b for a, b in zip(bsD, bs8))
    p4lo, p4hi = pct(dd, 0.025), pct(dd, 0.975)
    P4 = "줄었다(예측대로)" if p4lo > 0 else ("커졌다(예측 반대)" if p4hi < 0 else "구별 안 됨")
    print("-- P4(부 예측) 크기: 밀도 %s %+.3f · 밀도 8 기준(40 배경) %s · 차 %+.3f [%+.3f, %+.3f] → %s" % (D, p1m, fmt(p8m, p8lo, p8hi), p1m - p8m, p4lo, p4hi, P4))
    OUT.update({"P1": {"m": p1m, "ci": [p1lo, p1hi], "pass": P1}, "P2": {"ok": pair_ok, "n": npairs, "neg_point": pair_neg, "pass": P2},
                "P4": {"d8": [p8m, p8lo, p8hi], "diff": [p1m - p8m, p4lo, p4hi], "read": P4}})

    # ================= L 사다리
    print("\n== L — 교차 μ 사다리 · 밀도 %s · 씨앗 rascld 대 상주 racld" % D)
    dL_, pL = ladder_deltas(FL, SEEDS)
    SD = LSet("D", dL_)
    S8L = LSet("d8pool", d8) if ok8 else None
    print("   %6s | %-26s | %-26s || rsacld · racldx (밀도 %s)" % ("μ", "씨앗 Δ 밀도 %s (n=%d)" % (D, SD.n), "씨앗 Δ 밀도 8 (n=%d)" % (S8L.n if S8L else 0), D))
    for mu in MUS:
        a = SD.bmean("rascld", mu); b8 = S8L.bmean("rascld", mu) if S8L else (float("nan"),) * 3
        print("   %5.2f%% | %s | %s || %+.3f · %+.3f" % (mu * 100, fmt(*a), fmt(*b8), SD.mean("rsacld", mu), SD.mean("racldx", mu)))
    d0, d1 = SD.bmean("rascld", 0.0), SD.bmean("rascld", 0.01)
    L1 = d0[1] > 0 and d1[2] < 0
    fD, bD = SD.flips()
    f8, b8s = S8L.flips() if S8L else (None, [])
    diffs = [x - y for x, y in zip(bD, b8s) if not (math.isinf(x) and math.isinf(y))]
    if not diffs:
        print("🚨 교차 비교 불가(기준 없음)"); sys.exit(1)
    pu = sum(1 for x in diffs if x > 0) / N_BOOT; pd = sum(1 for x in diffs if x < 0) / N_BOOT   # 분모 = N_BOOT(R6)
    n_infinf = N_BOOT - len(diffs); n_tie = sum(1 for x in diffs if x == 0)
    L2 = "✅ 예측대로 위로" if pu >= P_DIR else ("❌ 반대(아래로)" if pd >= P_DIR else "🟡 구별 안 됨")
    cens_hi = sum(1 for x in bD if math.isinf(x)) / N_BOOT; cens_lo = sum(1 for x in bD if x == 0.0) / N_BOOT
    print("\n-- L1 반전: Δ(0) %s · Δ(1%%) %s → %s" % (fmt(*d0), fmt(*d1), "✅" if L1 else "❌"))
    print("-- L2 교차: 밀도 %s %s [%s, %s] (재추출 범위 위 ∞ %.1f%% · 0 아래 %.1f%%) · 밀도 8 %s [%s, %s]"
          % (D, fmt_mu(fD), fmt_mu(pct(bD, 0.025)), fmt_mu(pct(bD, 0.975)), 100 * cens_hi, 100 * cens_lo, fmt_mu(f8), fmt_mu(pct(b8s, 0.025)), fmt_mu(pct(b8s, 0.975))))
    dpt = (fD - f8) if (fD is not None and f8 is not None) else None
    print("   차(밀도 %s − 8) 점 %s · 재추출 [%s, %s] · P(위로) %.3f · P(아래로) %.3f → **%s**"
          % (D, "∞" if dpt is not None and math.isinf(dpt) else ("%+.3f%%p" % (100 * dpt) if dpt is not None else "-"),
             "∞" if math.isinf(pct(diffs, 0.025)) else "%+.3f" % (100 * pct(diffs, 0.025)), "∞" if math.isinf(pct(diffs, 0.975)) else "%+.3f" % (100 * pct(diffs, 0.975)), pu, pd, L2))
    print("   분모 %d(R6) · 어느 쪽에도 안 센 재추출: ∞ − ∞ %d · 동률(d = 0) %d" % (N_BOOT, n_infinf, n_tie))
    dhi = pct(diffs, 0.975)
    excl_low = (not math.isinf(dhi)) and dhi < (PRED_BAND[0] - f8 if f8 is not None else -math.inf)
    print("   서술(R2): 차 d 의 95%% 구간 [%s, %s]%%p · §3 띠 하단(0.21%%)이 배제되나(d 상한 < 0.21%% − 밀도 8 교차) — %s"
          % ("∞" if math.isinf(pct(diffs, 0.025)) else "%+.3f" % (100 * pct(diffs, 0.025)), "∞" if math.isinf(dhi) else "%+.3f" % (100 * dhi), "예" if excl_low else "아니오"))
    inband = fD is not None and PRED_BAND[0] <= fD <= PRED_BAND[1]
    print("   서술: §3 수치 예측 띠 [0.21%%, 0.65%%] 안인가 — %s" % ("예" if inband else "아니오"))
    # 서술: 괄호 · 기울기 · 곡률 · 순도
    e1, e3 = SD.bmean("rascld", 0.001), SD.bmean("rascld", 0.003)
    lines = {}
    for c in LCODES:
        (a0, b1), pts = SD.curve(c, 1); lines[c] = (b1, [p[1] for p in pts])
    slope8 = S8L.curve("rascld", 1)[0][1] if S8L else (float("nan"),) * 3
    (q0, q1, q2), _ = SD.curve("rascld", 2)
    print("-- 서술: 괄호 Δ(0.1%%) %s · Δ(0.3%%) %s" % (fmt(*e1), fmt(*e3)))
    print("   기울기 씨앗 %+.1f [%+.1f, %+.1f] (밀도 8 %+.1f) · rsacld %+.1f [%+.1f, %+.1f] · racldx %+.1f [%+.1f, %+.1f] · 곡률 c %+.0f [%+.0f, %+.0f]"
          % (*lines["rascld"][0], slope8[0], *lines["rsacld"][0], *lines["racldx"][0], *q2))
    purity = {}
    for mu in MUS:
        row = {}
        for c in ["ref"] + LCODES:
            v = [x for x in pL[c][mu] if not math.isnan(x)]; row[c] = sum(v) / len(v) if v else float("nan")
        purity[str(mu)] = row
    order = all(purity[str(mu)]["ref"] > purity[str(mu)]["racldx"] > purity[str(mu)]["rsacld"] > purity[str(mu)]["rascld"] for mu in MUS if mu > 0)
    print("   순도 WT(ref) > racldx > rsacld > 씨앗 (μ > 0 전부): %s" % ("같음" if order else "다름"))
    OUT.update({"L1": {"d0": list(d0), "d1pct": list(d1), "pass": L1},
                "L2": {"flip_D": fD, "flip_D_ci": [pct(bD, 0.025), pct(bD, 0.975)], "flip_8": f8, "flip_8_ci": [pct(b8s, 0.025), pct(b8s, 0.975)],
                       "diff": dpt, "diff_ci": [pct(diffs, 0.025), pct(diffs, 0.975)], "p_up": pu, "p_down": pd, "cens_hi": cens_hi, "cens_lo": cens_lo, "read": L2,
                       "in_pred_band": inband, "denom": N_BOOT, "n_inf_inf": n_infinf, "n_tie": n_tie, "pred_low_excluded": excl_low},
                "ladder_delta_D": {c: {str(mu): list(SD.bmean(c, mu)) for mu in MUS} for c in LCODES},
                "ladder_desc": {"bracket": [list(e1), list(e3)], "slopes": {c: list(lines[c][0]) for c in LCODES}, "slope_rascld_d8": list(slope8),
                                "curv": list(q2), "purity": purity, "purity_order_same": order}})

    # ================= 조작 확인(서술)
    medD, pooledD, a1D, nsD = ib_summary(IBD, rcodes, "D")
    print("\n== 조작 확인(서술 · invbud · racld 의 n 삽입 여섯 · 고리 4 둘 · 5 넷) — 시도/성공 R 중앙: 밀도 %s %.2f · 밀도 8 %.2f · ln R(4) − ln R(5): 밀도 %s %s · 밀도 8 %s"
          % (D, medD, med8, D, fmt(*a1D), fmt(*a18)))
    # R4 — 조작 확인 규칙(서술 규칙 · 파일럿 뒤 · 본 데이터 전 고정)
    def lnR_by_bg(IB):
        cl = list(rcodes)
        ss = sorted(set.intersection(*[set(IB[c]) for c in cl])) if all(IB[c] for c in cl) else []
        return [sum(math.log(IB[c][s]["R"]) for c in cl) / len(cl) for s in ss]
    vD, v8 = lnR_by_bg(IBD), lnR_by_bg(IB8)
    if len(vD) >= 2 and len(v8) >= 2:
        mRD, _, _, bRD = boot_mean(vD, "manipR|D"); mR8, _, _, bR8 = boot_mean(v8, "manipR|d8")
        dR = [x - y for x, y in zip(bRD, bR8)]; rlo, rhi = pct(dR, 0.025), pct(dR, 0.975)
        Rread = "§3 전제만큼 줄었다" if rhi < math.log(0.5) else ("줄었지만 §3 전제(비 ≤ 0.5)에 못 미친다" if rhi < 0 else "줄지 않았다")
    else:
        mRD = mR8 = rlo = rhi = float("nan"); Rread = "셀 수 없음"
    print("   조작 확인 규칙(R4): 배경별 ln R 평균(코드 %d) 밀도 %s %+.3f (배경 %d) · 밀도 8 %+.3f (배경 %d) · 차 %+.3f [%+.3f, %+.3f] (비 e^d %.2f) → R 이 **%s**%s"
          % (len(rcodes), D, mRD, len(vD), mR8, len(v8), mRD - mR8, rlo, rhi, math.exp(mRD - mR8) if mRD == mRD else float("nan"), Rread,
             "" if Rread.startswith("§3 전제만큼") else " — §5: L2 를 '부족이 덜한 세계' 로 읽지 않는다"))
    pops = [FL[(0.0, s)]["grow_pop"] for s in SEEDS if (0.0, s) in FL]
    shares = [dict(FL[(0.0, s)]["grow_top"]).get("racld", 0) / FL[(0.0, s)]["grow_pop"] for s in SEEDS if (0.0, s) in FL]
    pops8 = [F8[("inv21", 0.0, s)]["grow_pop"] for s in new_bgs] if ok8 else []
    md = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
    print("   배경 개체 수 중앙: 밀도 %s %s · 밀도 8(해부 21 40) %s · racld 몫 중앙 %.3f" % (D, md(pops), md(pops8), md(shares)))
    OUT["manip"] = {"lnR_rule": {"D": [mRD, len(vD)], "d8": [mR8, len(v8)], "diff": [mRD - mR8, rlo, rhi], "read": Rread},
                    "R_median_D": medD, "R_median_8": med8, "R_pooled_D": pooledD, "A1_D": list(a1D), "A1_8": list(a18), "pop_D": md(pops), "pop_8": md(pops8), "racld_share_D": md(shares)}

    # ================= 종합
    core1 = P1 and P2
    print("\n== 종합 (밀도 %s)" % D)
    print("  ① 짝지은 고리 효과: P1 %s · P2 %s → %s" % ("✅" if P1 else "❌", "✅" if P2 else "❌", "선다" if core1 else "서지 않는다"))
    print("  ② 사다리 반전: L1 %s%s" % ("✅" if L1 else "❌", "" if L1 or not (fD is not None and math.isinf(fD)) else " (교차가 1% 위 — §5 읽기 참조)"))
    print("  교차 방향 예측(위로): %s" % L2)
    print("  부 예측 크기(줄어듦): %s" % P4)
    if PILOT: print("  🟡 파일럿 — 판정 아님(배경 %d)" % len(SEEDS))
    OUT["summary"] = {"core_pairs": core1, "P1": P1, "P2": P2, "L1": L1, "L2": L2, "P4": P4}
    OUT["bad"] = bad; OUT["notes"] = notes
    fn = "_RESULT_pilot25.json" if PILOT else "_RESULT_pairs25.json"
    json.dump(clean(OUT), open(os.path.join(BASE, fn), "w"), ensure_ascii=False, indent=1, allow_nan=False)
    print("기록 →", fn)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "qual":
        qual_mode(sys.argv[2:])
    else:
        judge()
