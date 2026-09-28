#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 11 판정 — 수명 사다리의 시간 예산(A) · 밀도 개입(A · L) · 침입 Δ 의 수명 의존(B) · 고리 길이 사다리(C)
사전등록: 해부11-사전등록-2026-09-20.md §2~§5. 판정선은 결과 전에 고정했다. 부트스트랩은 시드(배경) 단위 2,000회.

  A  budget.js v1.2.0 · μ 0 · 창 = t0 ≥ 10000 인 칸.  헛손질/성공 s = stall/ok 를 **정확히** 셋으로 가른다:
       ln s = ln(헛손질 시간 몫 stall/stepped) + ln(출생당 개체-틱 stepped/births) − ln(출생당 글자 ok/births)
     A1 Δln s(1200 − 300) 하한 > 0 (밀도 8 · 세 코드 모두)
     A2 인구학 몫 = Δln(stepped/births) / Δln s  — 점추정 ≥ 0.8 (세 코드 모두) → '인구학이 주'
     A3 |Δln 헛손질 시간 몫| 의 구간이 ±0.2 안 (세 코드 모두) → '시간 몫은 수명에 둔감'
     A5a 밀도 개입: Δln 시간 몫(DHI − 8) 상한 < 0 (수명 300 · 1200 × 세 코드) · A5b 차의 차 [Δln s(DHI)] − [Δln s(8)] 상한 < 0 (세 코드) → '자리가 더 묶는 접시에선 s 의 수명 의존이 더 작다'
  L  lineage --spec · μ 0.1% · I = R × F (해부 8 K0 와 같은 정의)
     L1 밀도 8 에서 Δln I(1200 − 300) 하한 > 0 (E1 재현 · 새 시드 · 길이를 늘린 접시) · L2 차의 차 [Δln I(DHI)] − [Δln I(8)] 의 구간 — 코드마다 0 을 빼는지와 방향(양쪽으로 연다) · L3 ln I(DHI) − ln I(8) 상한 < 0 (두 수명)
  B  dms --ancestor racld --age0. 수명 300 짝 = inv10anc. 배경은 뿌리마다 따로 부트스트랩(짝 아님)
     배경 기준: 멸종 안 함 · 끝 우세 코드가 racld 또는 그 글자 순열(acrld 등 — 파일럿: 수명 1200 · μ 1% 배경은 순열 구름이다)
     B1 수명 1200 의 μ 0 씨앗 Δ 구간 · B2 기울기 차 b(1200) − b(300) 구간이 0 을 빼면 '수명 의존 있음'(방향 함께) · B3 Δ = 0 이 되는 μ* (서술)
  C  dms · μ 0 · 고리 길이 사다리
     C1 배경별 Δ 를 고리 틱 수 {2 · 3 · 4 · 5}에 회귀한 기울기의 상한 < 0 이고 평균이 단조 → '고리가 짧을수록 유리'
     C2 Δ(racldr) − Δ(racldx) 구간(묶인 글자 하나의 값) · C3 Δ(racldn) − Δ(racldx) 구간(한 틱 + n) — 서술
실행: python3 an11_analyze.py [base]   (환경변수로 뿌리 · 시드 · 밀도를 덮어쓴다 — 마른 실행용)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv11")
ANC = os.environ.get("PETRI_ROOT_ANC", "inv10anc")
ilist = lambda k, d: [int(x) for x in os.environ[k].split(",")] if os.environ.get(k) else d
SEEDS = ilist("PETRI_SEEDS", list(range(1401, 1421)))
LSEEDS = ilist("PETRI_LSEEDS", list(range(1421, 1441)))
AGES = ilist("PETRI_AGES", [300, 600, 1200])
BAGES = ilist("PETRI_BAGES", [600, 1200])
DENS = ilist("PETRI_DENS", [8, 14])
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
PARTS = os.environ.get("PETRI_PARTS", "A,L,B,C").split(",")
CODES = ["racld", "rascld", "rsacld"]
LOOP = {"rascld": 2, "rsacld": 3, "racldx": 4, "nracld": 5}
CCODES = ["rascld", "rsacld", "racldx", "racldn", "nracld", "racldr"]
LET = "nsracldjehx"
NOM = 2 + 30 / 11
WIN0 = 10000
N_BOOT = 2000
RNG = random.Random(20260920)
A_LO, A_HI = AGES[0], AGES[-1]
D_LO, D_HI = DENS[0], DENS[-1]
OUT = {}
bad = []


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def ci(v):
    return [pct(v, .025), pct(v, .975)]


def fmt(m, c, d=3):
    return ("%+." + str(d) + "f [%+." + str(d) + "f, %+." + str(d) + "f]") % (m, c[0], c[1])


def boot(n):
    return [[RNG.randrange(n) for _ in range(n)] for _ in range(N_BOOT)]


def halt(msg):
    print("🚨 점검 실패 — 판정하지 않는다: " + msg)
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)


# ---------------------------------------------------------------- A · 시간 예산
def part_A():
    print("\n== A — 수명 사다리의 시간 예산 (μ 0 · 창 t0 ≥ %d × 수명/300 · 접시 길이 2만 × 수명/300)" % WIN0)
    D = {}
    for code in CODES:
        for d in DENS:
            for age in AGES:
                for s in SEEDS:
                    p = os.path.join(BASE, ROOT + "bud", "bud_mat_d%02d_%s_s%05d_a%d.json" % (d, code, s, age))
                    if not os.path.exists(p):
                        bad.append("없음 " + p); continue
                    z = json.load(open(p))
                    if z["budget_version"] != "1.2.0" or z["sim_version"] != "0.3.5" or z["opts"].get("age0") != age or z["opts"].get("ageVar") != age or abs(z["opts"]["density"] - d) > 1e-9 or z["ticks_planned"] != 20000 * age // 300 or z["bin"] != 500 * age // 300:
                        bad.append("버전/opts/길이 " + p); continue
                    D[(code, d, age, s)] = z
    n_exp = len(CODES) * len(DENS) * len(AGES) * len(SEEDS)
    print("  파일 %d/%d · 점검 문제 %d" % (len(D), n_exp, len(bad)))
    if bad or len(D) != n_exp:
        halt("A 파일")

    def row(z):
        B = [b for b in z["bins"] if b["t0"] >= WIN0 * z["opts"]["age0"] // 300]
        g = lambda k: sum(b[k] for b in B)
        code = z["code"]; used = sorted({LET.index(ch) for ch in code})
        tot = [z["bins"][0]["free0"][k] + sum(1 for ch in code if LET.index(ch) == k) for k in range(len(LET))]
        fr = sum(sum(b["free0"][k] for k in used) for b in B) / (len(B) * sum(tot[k] for k in used))
        return {"ok": g("copy_ok"), "stall": g("copy_stall"), "blocked": g("alloc_blocked"), "stepped": g("stepped"), "births": g("births"),
                "deaths": g("deaths"), "dage": g("death_age_sum"), "pop": sum(b["pop0"] for b in B) / len(B), "free_used": fr, "unknown": g("unknown") + g("alloc_dead")}

    R = {k: row(z) for k, z in D.items()}
    res = {}
    for code in CODES:
        for d in DENS:
            # 창업 실패 · 멸종 시드는 그 코드 × 밀도의 사다리 전체에서 뺀다
            keep = [s for s in SEEDS if all(D[(code, d, a, s)]["extinct_at"] < 0 and R[(code, d, a, s)]["births"] > 0 and R[(code, d, a, s)]["ok"] > 0 for a in AGES)]
            res[(code, d)] = keep
    if min(len(v) for v in res.values()) < max(1, len(SEEDS) // 2):
        halt("A 남은 시드가 절반 미만: %s" % {k: len(v) for k, v in res.items()})

    def agg(code, d, age, seeds, ix=None):
        rows = [R[(code, d, age, s)] for s in seeds]
        if ix is not None:
            rows = [rows[i] for i in ix]
        t = {k: sum(r[k] for r in rows) for k in ("ok", "stall", "blocked", "stepped", "births", "deaths", "dage", "unknown")}
        t["pop"] = sum(r["pop"] for r in rows) / len(rows); t["free_used"] = sum(r["free_used"] for r in rows) / len(rows)
        return {"s": t["stall"] / t["ok"], "share": t["stall"] / t["stepped"], "tpb": t["stepped"] / t["births"], "lpb": t["ok"] / t["births"],
                "blocked": t["blocked"] / t["stepped"], "life": t["dage"] / t["deaths"], "pop": t["pop"], "free_used": t["free_used"],
                "bal": (t["births"] - t["deaths"]) / t["births"], "unknown": t["unknown"] / t["stepped"]}

    print("  남은 시드: " + " · ".join("%s d%d %d" % (c, d, len(v)) for (c, d), v in res.items()))
    print("\n  코드 · 밀도 · 수명 | 개체 수 | 죽은 나이 | 헛손질/성공 s | 헛손질 시간 몫 | 막힘 시간 몫 | 출생당 개체-틱 | 출생당 글자 | 쓰는 글자 자유 몫 | (출생−죽음)/출생")
    cells = {}
    worst_bal = 0; worst_unknown = 0
    for code in CODES:
        for d in DENS:
            for age in AGES:
                a = agg(code, d, age, res[(code, d)]); cells["%s|%d|%d" % (code, d, age)] = a
                worst_bal = max(worst_bal, abs(a["bal"])); worst_unknown = max(worst_unknown, a["unknown"])
                print("  %-6s d%-2d a%-4d | %6.0f | %6.0f | %7.2f | %.3f | %.3f | %7.0f | %5.2f | %.3f | %+.4f" % (code, d, age, a["pop"], a["life"], a["s"], a["share"], a["blocked"], a["tpb"], a["lpb"], a["free_used"], a["bal"]))
    OUT["A_cells"] = cells
    print("  A0 가드: |출생 − 죽음|/출생 최대 %.4f (선 0.02) · 모름 몫 최대 %.5f (선 0.001)" % (worst_bal, worst_unknown))
    if worst_bal >= 0.02 or worst_unknown >= 0.001:
        halt("A0 정상 상태 · 모름 가드")

    def contrast(code, d, a_lo, a_hi):
        seeds = res[(code, d)]; IX = boot(len(seeds))
        def one(ix):
            lo, hi = agg(code, d, a_lo, seeds, ix), agg(code, d, a_hi, seeds, ix)
            ds = math.log(hi["s"] / lo["s"]); dsh = math.log(hi["share"] / lo["share"]); dt = math.log(hi["tpb"] / lo["tpb"]); dl = math.log(hi["lpb"] / lo["lpb"])
            return ds, dsh, dt, dl, dt / ds if ds else float("nan"), math.log(hi["pop"] / lo["pop"]), math.log(hi["free_used"] / lo["free_used"]), math.log(hi["life"] / lo["life"])
        p = one(None); bs = [one(ix) for ix in IX]
        return [(p[j], ci([b[j] for b in bs])) for j in range(8)]

    print("\n  Δln (수명 %d − %d) · 정확한 분해 Δln s = Δln 시간 몫 + Δln 출생당 개체-틱 − Δln 출생당 글자" % (A_HI, A_LO))
    a1 = a2 = a3 = True; a2hi = True
    for d in DENS:
        for code in CODES:
            c = contrast(code, d, A_LO, A_HI)
            OUT["A_contrast|%s|%d" % (code, d)] = c
            print("  %-6s d%-2d | Δln s %s | 시간 몫 %s | 출생당 개체-틱 %s | 출생당 글자 %s | 인구학 몫 %s | 닫힘 %.1e" % (
                code, d, fmt(*c[0]), fmt(*c[1]), fmt(*c[2]), fmt(*c[3]), fmt(*c[4], d=2), abs(c[0][0] - c[1][0] - c[2][0] + c[3][0])))
            print("           | Δln 개체 수 %s | Δln 쓰는 글자 자유 몫 %s | Δln 죽은 나이 %s" % (fmt(*c[5]), fmt(*c[6]), fmt(*c[7])))
            if d == D_LO:
                a1 &= c[0][1][0] > 0; a2 &= c[4][0] >= 0.8; a3 &= (c[1][1][0] > -0.2 and c[1][1][1] < 0.2)
    mono = all(cells["%s|%d|%d" % (code, D_LO, AGES[i])]["s"] < cells["%s|%d|%d" % (code, D_LO, AGES[i + 1])]["s"] for code in CODES for i in range(len(AGES) - 1))
    print("\n  A1 %s — Δln s 하한 > 0 (밀도 %d · 세 코드) · 사다리 단조 %s" % ("✅ 헛손질/성공이 수명과 함께 는다" if a1 else "❌", D_LO, mono))
    print("  A2 %s — 인구학 몫 점추정 ≥ 0.8 (세 코드)" % ("✅ 인구학(출생당 개체-틱)이 주" if a2 else "❌ 인구학 몫이 0.8 에 못 미치는 코드가 있다"))
    print("  A3 %s — |Δln 헛손질 시간 몫| 구간이 ±0.2 안 (세 코드)" % ("✅ 헛손질에 쓰는 시간 몫은 수명에 둔감" if a3 else "❌ 시간 몫도 수명에 따라 움직인다"))
    a5 = True
    print("\n  A5 밀도 개입 Δln 헛손질 시간 몫 (밀도 %d − %d)" % (D_HI, D_LO))
    for age in (A_LO, A_HI):
        for code in CODES:
            seeds = [s for s in res[(code, D_LO)] if s in res[(code, D_HI)]]; IX = boot(len(seeds))
            f = lambda ix: math.log(agg(code, D_HI, age, seeds, ix)["share"] / agg(code, D_LO, age, seeds, ix)["share"])
            g = lambda ix: math.log(agg(code, D_HI, age, seeds, ix)["s"] / agg(code, D_LO, age, seeds, ix)["s"])
            m, c_ = f(None), ci([f(ix) for ix in IX]); m2, c2 = g(None), ci([g(ix) for ix in IX])
            OUT["A5|%s|%d" % (code, age)] = {"share": [m, c_], "s": [m2, c2]}
            a5 &= c_[1] < 0
            print("  %-6s a%-4d | 시간 몫 %s | s %s" % (code, age, fmt(m, c_), fmt(m2, c2)))
    for code in CODES:
        seeds = [s for s in res[(code, D_LO)] if s in res[(code, D_HI)]]; IX = boot(len(seeds))
        f = lambda ix: math.log(agg(code, D_HI, A_HI, seeds, ix)["s"] / agg(code, D_HI, A_LO, seeds, ix)["s"]) - math.log(agg(code, D_LO, A_HI, seeds, ix)["s"] / agg(code, D_LO, A_LO, seeds, ix)["s"])
        m, c_ = f(None), ci([f(ix) for ix in IX]); OUT["A5b|%s" % code] = [m, c_]; a2hi &= c_[1] < 0
        print("  %-6s 차의 차 [Δln s(밀도 %d)] − [Δln s(%d)]: %s" % (code, D_HI, D_LO, fmt(m, c_)))
    print("  A5a %s · A5b %s" % ("✅ 넉넉한 재료는 헛손질 시간 몫을 줄인다" if a5 else "❌ 시간 몫 감소가 모든 칸에서 서지 않는다", "✅ 자리가 더 묶는 접시에선 s 의 수명 의존이 더 작다" if a2hi else "❌ 차의 차가 모든 코드에서 서지 않는다"))
    OUT["A_verdict"] = {"A1": a1, "A1_mono": mono, "A2": a2, "A3": a3, "A5a": a5, "A5b": a2hi}


# ---------------------------------------------------------------- L · 부풀림 직접
def part_L():
    print("\n== L — 부풀림 직접 (lineage --spec · μ 0.1%)")
    D = {}
    for code in CODES:
        for d in DENS:
            for age in (A_LO, A_HI):
                for s in LSEEDS:
                    p = os.path.join(BASE, "%slin_d%d" % (ROOT, d), "lin_%s_mu0p001_s%05d_a%d.json" % (code, s, age))
                    if not os.path.exists(p):
                        bad.append("없음 " + p); continue
                    z = json.load(open(p))
                    if z["sim_version"] != "0.3.5" or z["lineage_version"] != "1.5.0" or z["spec"] is None or z["opts"].get("age0") != age or abs(z["opts"]["density"] - d) > 1e-9 or abs(z["mu"] - 0.001) > 1e-12 or z["ticks_planned"] != 20000 * age // 300 or z["win"][0] != 5000 * age // 300:
                        bad.append("버전/opts " + p); continue
                    D[(code, d, age, s)] = z
    n_exp = len(CODES) * len(DENS) * 2 * len(LSEEDS)
    miss = max((z["missed_births"] / max(z["births"], 1) for z in D.values()), default=1)
    print("  파일 %d/%d · 점검 문제 %d · 놓친 출생 몫 최대 %.3f%%" % (len(D), n_exp, len(bad), 100 * miss))
    if bad or len(D) != n_exp:
        halt("L 파일")

    def I_of(rows):
        ok = sum(r["spec"]["ok"] for r in rows); de = sum(r["spec"]["del"] for r in rows)
        su = sum(sum(r["spec"]["sub"]) for r in rows); ins = sum(sum(r["spec"]["ins"]) for r in rows); st = sum(r["spec"]["stall"] for r in rows)
        adv = ok + su + de
        R = (de / adv) / (0.001 / 3); Fv = (1 + ins / de + su / de) / NOM
        return R * Fv, R, Fv, st / max(1, ok)

    keep = {}
    for code in CODES:
        keep[code] = [s for s in LSEEDS if all(D[(code, d, a, s)]["extinct_at"] < 0 and D[(code, d, a, s)]["spec"]["del"] > 0 for d in DENS for a in (A_LO, A_HI))]
    if min(len(v) for v in keep.values()) < max(1, len(LSEEDS) // 2):
        halt("L 남은 시드가 절반 미만: %s" % {k: len(v) for k, v in keep.items()})
    l1 = l3 = True; l2 = []
    for code in CODES:
        seeds = keep[code]; IX = boot(len(seeds))
        rows = lambda d, a, ix: [D[(code, d, a, seeds[i])] for i in (ix if ix is not None else range(len(seeds)))]
        line = "  %-6s (시드 %d)" % (code, len(seeds))
        for d in DENS:
            for a in (A_LO, A_HI):
                I, R, Fv, so = I_of(rows(d, a, None)); OUT["L_cell|%s|%d|%d" % (code, d, a)] = {"I": I, "R": R, "F": Fv, "stall_per_ok": so}
                line += " | d%d a%d I %.2f (R %.2f · F %.2f · s %.1f)" % (d, a, I, R, Fv, so)
        print(line)
        for d in DENS:
            f = lambda ix: math.log(I_of(rows(d, A_HI, ix))[0] / I_of(rows(d, A_LO, ix))[0])
            m, c_ = f(None), ci([f(ix) for ix in IX]); OUT["L_age|%s|%d" % (code, d)] = [m, c_]
            print("           Δln I(수명 %d − %d) · 밀도 %d: %s" % (A_HI, A_LO, d, fmt(m, c_)))
            if d == D_LO: l1 &= c_[0] > 0
        f = lambda ix: math.log(I_of(rows(D_HI, A_HI, ix))[0] / I_of(rows(D_HI, A_LO, ix))[0]) - math.log(I_of(rows(D_LO, A_HI, ix))[0] / I_of(rows(D_LO, A_LO, ix))[0])
        m, c_ = f(None), ci([f(ix) for ix in IX]); OUT["L_did|%s" % code] = [m, c_]; l2.append("%s %s" % (code, "더 작다" if c_[1] < 0 else ("더 크다" if c_[0] > 0 else "검출 안 됨")))
        print("           차의 차 [Δln I(밀도 %d)] − [Δln I(%d)]: %s" % (D_HI, D_LO, fmt(m, c_)))
        for a in (A_LO, A_HI):
            f = lambda ix: math.log(I_of(rows(D_HI, a, ix))[0] / I_of(rows(D_LO, a, ix))[0])
            m, c_ = f(None), ci([f(ix) for ix in IX]); OUT["L_dens|%s|%d" % (code, a)] = [m, c_]
            print("           ln I(밀도 %d) − ln I(%d) · 수명 %d: %s" % (D_HI, D_LO, a, fmt(m, c_)))
            l3 &= c_[1] < 0
    print("  L1 %s · L2 %s · L3 %s" % ("✅ E1 재현(새 시드)" if l1 else "❌ E1 재현 안 됨", "밀도 %d 에서 수명 의존이 — %s" % (D_HI, " · ".join(l2)),
                                  "✅ 넉넉한 재료가 부풀림을 줄인다" if l3 else "❌ 밀도 효과가 모든 칸에서 서지 않는다"))
    OUT["L_verdict"] = {"L1": l1, "L2": l2, "L3": l3}


# ---------------------------------------------------------------- B · C 공용
def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def load_dms(root, mus, expect_age, expect_anc, n_arms, want_main, perm_ok=False):
    F, dropped = {}, {}; load_dms.exact = 0
    for mu in mus:
        for f in sorted(glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f)); s = z["seed"]
            if abs(z["mu_assay"] - mu) > 1e-12: bad.append("μ %s %s %d" % (root, mu, s))
            if z.get("age0") != expect_age or z.get("ancestor") != expect_anc or z.get("find_first") or z.get("cosmic_off"): bad.append("규칙 표시 %s %s %d" % (root, mu, s))
            if expect_age and (z["ticks"] != 3000 * expect_age // 300 or z["grow"] != 50000 * expect_age // 300): bad.append("길이 배수 %s %s %d" % (root, mu, s))
            if not z["fidelity"]["same"]: bad.append("충실도 %s %s %d" % (root, mu, s))
            if want_main and not z["checksum_match"]: bad.append("배경 체크섬 %s %s %d" % (root, mu, s))
            if z["grow_extinct"] >= 0:
                dropped[s] = "멸종"; continue
            top = z["grow_top"][0][0] if z["grow_top"] else ""
            if top == "racld": load_dms.exact += 1
            if not (top == "racld" or (perm_ok and sorted(top) == sorted("racld"))):
                dropped[s] = "끝 우세 %s" % (z["grow_top"][0][0] if z["grow_top"] else "-"); continue
            if len(z["arms"]) != n_arms: bad.append("팔 수 %s %s %d" % (root, mu, s))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %s %d" % (root, mu, s))
            F[(mu, s)] = z
    bgs = sorted({k[1] for k in F}); bgs = [s for s in bgs if all((mu, s) in F for mu in mus)]
    return F, bgs, dropped


def deltas(F, bgs, mus, codes):
    d = {c: {mu: [] for mu in mus} for c in codes}
    for mu in mus:
        for s in bgs:
            z = F[(mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                if a["kind"] == "mut" and a["key"] in d:
                    d[a["key"]][mu].append(s_of(a) - base)
    return d


def part_B():
    print("\n== B — 침입 Δ 의 수명 의존 (둘째 시조 세계 · 길이를 수명 배수로 늘림)")
    roots = {300: ANC}
    for a in BAGES: roots[a] = "%sage%d" % (ROOT, a)
    DL, BG, IXS = {}, {}, {}
    for a, root in roots.items():
        F, bgs, dropped = load_dms(root, MUS, None if a == 300 else a, "racld", 14, False, perm_ok=True)
        exact = load_dms.exact // len(MUS)
        BG[a] = bgs; DL[a] = deltas(F, bgs, MUS, ["rascld", "rsacld", "racldx"]); IXS[a] = boot(max(1, len(bgs)))
        print("  수명 %d · 배경 %d (끝 우세가 정확히 racld 인 것 %d · 뺀 것 %d: %s)" % (a, len(bgs), exact, len(dropped), ", ".join("%d %s" % kv for kv in sorted(dropped.items())) or "없음"))
    if bad or min(len(v) for v in BG.values()) < max(1, NBG // 2):
        halt("B 배경이 절반 미만이거나 점검 문제")
    mx = sum(MUS) / len(MUS)
    def fit(a, c, ix=None):
        ys = []
        for mu in MUS:
            v = DL[a][c][mu]; ys.append(sum(v[i] for i in ix) / len(ix) if ix is not None else sum(v) / len(v))
        my = sum(ys) / len(ys)
        b = sum((m - mx) * (y - my) for m, y in zip(MUS, ys)) / sum((m - mx) ** 2 for m in MUS)
        return b, my - b * mx, ys
    print("\n        μ | " + " | ".join("수명 %-6d" % a for a in roots))
    for j, mu in enumerate(MUS):
        print("  %7g | " % mu + " | ".join("%+.3f     " % fit(a, "rascld")[2][j] for a in roots))
    res = {}
    for a in roots:
        for c in ("rascld", "rsacld", "racldx"):
            b, i0, ys = fit(a, c); bs = [fit(a, c, ix) for ix in IXS[a]]
            res[(a, c)] = {"slope": [b, ci([x[0] for x in bs])], "mu0": [ys[0], ci([x[2][0] for x in bs])],
                           "mustar": [-i0 / b if b else float("nan"), ci([-x[1] / x[0] for x in bs if x[0]])], "_bs": [x[0] for x in bs]}
        r = res[(a, "rascld")]
        print("  수명 %-5d 씨앗: μ0 Δ %s · 기울기 %s · μ* %.4f [%.4f, %.4f] · rsacld 기울기 %s · racldx %s" % (
            a, fmt(*r["mu0"]), fmt(*r["slope"], d=1), r["mustar"][0], r["mustar"][1][0], r["mustar"][1][1], fmt(*res[(a, "rsacld")]["slope"], d=1), fmt(*res[(a, "racldx")]["slope"], d=1)))
    hi = BAGES[-1]
    b1 = res[(hi, "rascld")]["mu0"]
    diff = [x - y for x, y in zip(res[(hi, "rascld")]["_bs"], res[(300, "rascld")]["_bs"])]
    dm, dc = res[(hi, "rascld")]["slope"][0] - res[(300, "rascld")]["slope"][0], ci(diff)
    ratio = res[(hi, "rascld")]["slope"][0] / res[(300, "rascld")]["slope"][0]
    v1 = "μ 0 우위가 남는다" if b1[1][0] > 0 else ("μ 0 우위가 뒤집혔다" if b1[1][1] < 0 else "μ 0 우위가 검출되지 않는다")
    v2 = ("수명 의존 있음 — 긴 수명에서 불이익이 %s" % ("더 가파르다" if dm < 0 else "더 완만하다")) if (dc[0] > 0 or dc[1] < 0) else "기울기의 수명 의존은 검출되지 않는다"
    print("  B1 수명 %d 의 μ 0 씨앗 Δ %s → %s" % (hi, fmt(*b1), v1))
    print("  B2 기울기 차 b(%d) − b(300) = %s · 비 %.2f → %s" % (hi, fmt(dm, dc, d=1), ratio, v2))
    for k in res: res[k].pop("_bs")
    OUT["B"] = {"%d|%s" % k: v for k, v in res.items()}; OUT["B_verdict"] = {"B1": v1, "B2": v2, "diff": [dm, dc], "ratio": ratio, "bg": {str(a): BG[a] for a in BG}}


def part_C():
    print("\n== C — 고리 길이 사다리 (μ 0 · mat 배경)")
    F, bgs, dropped = load_dms(ROOT + "loop", [0.0], None, None, 11 + len(CCODES), True)
    print("  배경 %d (뺀 것 %d)" % (len(bgs), len(dropped)))
    if bad or len(bgs) < max(1, NBG // 2) or dropped:
        halt("C 배경")
    d = deltas(F, bgs, [0.0], CCODES); IX = boot(len(bgs))
    V = {c: d[c][0.0] for c in CCODES}
    m = lambda c, ix=None: sum(V[c][i] for i in ix) / len(ix) if ix is not None else sum(V[c]) / len(V[c])
    for c in CCODES:
        print("  %-6s 고리 %s | Δ %s" % (c, LOOP.get(c, "4*"), fmt(m(c), ci([m(c, ix) for ix in IX]))))
        OUT["C_delta|" + c] = [m(c), ci([m(c, ix) for ix in IX])]
    xs = sorted(LOOP, key=LOOP.get); lx = [LOOP[c] for c in xs]; mlx = sum(lx) / len(lx)
    def sl(ix=None):
        ys = [m(c, ix) for c in xs]; my = sum(ys) / len(ys)
        return sum((x - mlx) * (y - my) for x, y in zip(lx, ys)) / sum((x - mlx) ** 2 for x in lx)
    b, bc = sl(), ci([sl(ix) for ix in IX])
    mono = all(m(xs[i]) > m(xs[i + 1]) for i in range(len(xs) - 1))
    c1 = bc[1] < 0 and mono
    df = lambda a, b_: (m(a) - m(b_), ci([m(a, ix) - m(b_, ix) for ix in IX]))
    c2, c3 = df("racldr", "racldx"), df("racldn", "racldx")
    print("  C1 Δ 의 고리 틱당 기울기 %s · 평균 단조 %s → %s" % (fmt(b, bc), mono, "✅ 고리가 짧을수록 μ 0 침입에 유리하다" if c1 else "❌ 고리 길이 사다리가 서지 않는다"))
    print("  C2 묶인 글자 하나의 값 Δ(racldr) − Δ(racldx) = %s → %s" % (fmt(*c2), "희소성 몫 있음" if (c2[1][0] > 0 or c2[1][1] < 0) else "검출되지 않는다"))
    print("  C3 한 틱 + n 의 값 Δ(racldn) − Δ(racldx) = %s (서술)" % fmt(*c3))
    OUT["C_verdict"] = {"C1": c1, "slope": [b, bc], "mono": mono, "C2": list(c2), "C3": list(c3), "bg": bgs}


print("== 해부 11 판정 — 뿌리 %s · 부분 %s" % (ROOT, ",".join(PARTS)))
for p in PARTS:
    {"A": part_A, "L": part_L, "B": part_B, "C": part_C}[p]()
tag = os.environ.get("PETRI_OUT", "_RESULT_an11")
json.dump(OUT, open(os.path.join(BASE, tag + ".json"), "w"), ensure_ascii=False, indent=1)
print("\n기록 → %s.json" % tag)
