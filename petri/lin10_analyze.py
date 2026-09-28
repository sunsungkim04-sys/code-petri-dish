# -*- coding: utf-8 -*-
"""해부 10 판정 D · E — 확산 사다리(D) · 수명 사다리(E). 사전등록 해부10-사전등록-2026-09-18.md §3 · §4. 결과 전에 썼다(파일럿은 노트 §6).

  D  lineage.js v1.5.0 --spec · 밀도 8 · μ 1% · diffuse {0.05 · 0.15 · 0.5} × 코드 {racld · rascld · rsacld} × 시드 1201~1220 (inv10df/)
  E  lineage.js v1.5.0 --spec · 밀도 8 · μ {0.1% · 1%} · age0 = ageVar ∈ {300 · 600 · 1200} × 코드 3 × 시드 1301~1320 (inv10age/)
  칸마다 시드 합으로 R(주사위 몫) · F(거름 몫) · I = R × F · u_obs · 순도(끝 WT 몫) 를 내고 시드 재추출 2,000회 95%(같은 시드끼리 짝).
  멸종 접시는 그 시드를 그 코드 사다리 전체에서 뺀다(수 보고).

  D1 주사위: 세 코드 모두 ln R(0.05) − ln R(0.5) 하한 > 0 → '확산이 느릴수록 주사위 몫이 크다'
  D2 빠른 확산: 세 코드 R(0.5) 점추정이 [0.9, 1.5] 안 → '빠른 확산은 밀도 32 수준(1.0~1.3)' · 아니면 서술
  D3 씨앗 초과 공급: e = ln(u 씨앗/u racld) − ln S_len (S_len = u_nom(6)/u_nom(5)) · e(0.05) − e(0.5) 하한 > 0 → '확산이 느릴수록 씨앗의 초과 공급이 크다'
  E1 부풀림: 코드 × μ 여섯 칸 모두 ln I(1200) − ln I(300) 구간이 0 을 품는다 → '부풀림은 수명과 무관' · 하나라도 벗어나면 방향을 적는다
  E2 멸종: μ 1% 의 멸종 · 창업 실패 접시 수가 수명에 따라 준다(300 ≥ 600 ≥ 1200) → 서술
  E3 서술: 순도 ln(씨앗/racld) · 씨앗 초과 공급 e · 세대당 헛손질 — 수명별
  점검: 파일 완결 · sim 0.3.5 · lineage 1.5.0 · opts 일치 · 놓친 출생 서술
실행: python3 lin10_analyze.py [petri 폴더]  →  _RESULT_lin10.json
마른 실행: PETRI_ROOT_DF=stage10/pilot10/df PETRI_ROOT_AGE=stage10/pilot10/age PETRI_SEEDS=9501 PETRI_ASEEDS=9501 python3 lin10_analyze.py
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT_DF = os.environ.get("PETRI_ROOT_DF", "inv10df")
ROOT_AGE = os.environ.get("PETRI_ROOT_AGE", "inv10age")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(1201, 1221))
ASEEDS = [int(x) for x in os.environ["PETRI_ASEEDS"].split(",")] if os.environ.get("PETRI_ASEEDS") else list(range(1301, 1321))
CODES = ["racld", "rascld", "rsacld"]
DFS = [0.05, 0.15, 0.5]
AGES = [300, 600, 1200]
AMUS = [0.001, 0.01]
NOM = 2 + 30 / 11
N_BOOT = 2000
RNG = random.Random(20260923)


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def unom(mu, L):
    return 1 - (1 - mu * (2 / 3 + 10 / 11)) ** L


def dftag(df):
    return "_df" + str(df).replace(".", "p")


def load(path, expect):
    if not os.path.exists(path):
        return None, "없음 " + path
    z = json.load(open(path))
    if z["sim_version"] != "0.3.5" or z["lineage_version"] != "1.5.0" or z["spec"] is None:
        return None, "버전/spec " + path
    for k, v in expect.items():
        if abs(z["opts"].get(k, -1e9) - v) > 1e-9:
            return None, "opts %s %s" % (k, path)
    return z, None


D, bad = {}, []
for code in CODES:
    for df in DFS:
        for s in SEEDS:
            z, e = load(os.path.join(BASE, ROOT_DF, "lin_%s_mu0p01_s%05d%s.json" % (code, s, dftag(df))), {"diffuse": df, "mu": 0.01, "density": 8})
            if e: bad.append(e)
            else: D[("df", code, df, 0.01, s)] = z
    for age in AGES:
        for mu in AMUS:
            for s in ASEEDS:
                z, e = load(os.path.join(BASE, ROOT_AGE, "lin_%s_mu%s_s%05d_a%d.json" % (code, str(mu).replace(".", "p"), s, age)), {"age0": age, "ageVar": age, "mu": mu, "density": 8})
                if e: bad.append(e)
                else: D[("age", code, age, mu, s)] = z
total = len(CODES) * (len(DFS) * len(SEEDS) + len(AGES) * len(AMUS) * len(ASEEDS))
print("== 해부 10 D · E — 확산 사다리 · 수명 사다리 — 판정")
print("  파일 %d/%d · 점검 문제 %d · 놓친 출생 몫 최대 %.3f%%" % (len(D), total, len(bad), 100 * max((z["missed_births"] / max(z["births"], 1) for z in D.values()), default=1)))
if bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)


def stats(rows, code, mu):
    ok = sum(r["spec"]["ok"] for r in rows); de = sum(r["spec"]["del"] for r in rows)
    su = sum(sum(r["spec"]["sub"]) for r in rows); ins = sum(sum(r["spec"]["ins"]) for r in rows)
    st = sum(r["spec"]["stall"] for r in rows)
    wpb = sum(r["wt_parent_births"] for r in rows); wpm = sum(r["wt_parent_mut_births"] for r in rows)
    pur_n = sum(next((c for k, c in r["top_end"] if k == code), 0) for r in rows); pop = sum(r["series"][-1][1] for r in rows)
    L = len(code); adv = ok + su + de
    if de == 0 or adv == 0 or wpb == 0:
        return None
    R = (de / adv) / (mu / 3); Fv = (1 + ins / de + su / de) / NOM
    return {"R": R, "F": Fv, "I": R * Fv, "u": wpm / wpb, "I_lin": (wpm / wpb) / unom(mu, L), "pur": (pur_n + 0.5) / (pop + 0.5), "stall_per_ok": st / max(1, ok)}


def cell(kind, code, x, mu, seeds):
    rows = [D[(kind, code, x, mu, s)] for s in seeds if (kind, code, x, mu, s) in D]
    return rows


# 멸종 처리: 시드를 코드 사다리 전체에서 뺀다
def alive_seeds(kind, code, xs, mus, seeds):
    keep = []
    for s in seeds:
        if all(((kind, code, x, mu, s) in D and D[(kind, code, x, mu, s)]["extinct_at"] < 0) for x in xs for mu in mus):
            keep.append(s)
    return keep


out = {"D": {}, "E": {}}
print("\n-- D 확산 사다리 (μ 1% · 밀도 8)  R 주사위 · F 거름 · I · I_계보 · 순도 · 헛손질/맞게")
sd = {c: alive_seeds("df", c, DFS, [0.01], SEEDS) for c in CODES}
boot_d = {}
for c in CODES:
    print("   %-7s 살아남은 시드 %d/%d" % (c, len(sd[c]), len(SEEDS)))
    for df in DFS:
        rows = cell("df", c, df, 0.01, sd[c])
        st = stats(rows, c, 0.01)
        bs = [stats([rows[RNG.randrange(len(rows))] for _ in rows], c, 0.01) for _ in range(N_BOOT)]
        boot_d[(c, df)] = (st, [b for b in bs if b])
        print("      diffuse %-5g R %6.2f [%6.2f, %6.2f] · F %.2f · I %6.2f · I_계보 %6.2f · 순도 %.3f · 헛손질/맞게 %.1f"
              % (df, st["R"], pct([b["R"] for b in bs if b], .025), pct([b["R"] for b in bs if b], .975), st["F"], st["I"], st["I_lin"], st["pur"], st["stall_per_ok"]))
        out["D"]["%s|%g" % (c, df)] = st
# D1 · D2
d1 = {}
for c in CODES:
    (s1, b1), (s2, b2) = boot_d[(c, 0.05)], boot_d[(c, 0.5)]
    n = min(len(b1), len(b2))
    v = [math.log(b1[i]["R"]) - math.log(b2[i]["R"]) for i in range(n)]
    d1[c] = (math.log(s1["R"]) - math.log(s2["R"]), pct(v, .025), pct(v, .975))
ok1 = all(d1[c][1] > 0 for c in CODES)
print("-- D1 ln R(0.05) − ln R(0.5): " + " · ".join("%s %+.2f [%+.2f, %+.2f]" % (c, *d1[c]) for c in CODES) + " → **%s**" % ("확산이 느릴수록 주사위 몫이 크다" if ok1 else "세 코드 모두에서 검출되지는 않는다"))
r05 = {c: boot_d[(c, 0.5)][0]["R"] for c in CODES}
ok2 = all(0.9 <= r05[c] <= 1.5 for c in CODES)
print("-- D2 R(0.5): " + " · ".join("%s %.2f" % (c, r05[c]) for c in CODES) + " → **%s**" % ("빠른 확산은 밀도 32 수준(1.0~1.3)" if ok2 else "밀도 32 수준에 이르지 않는다(서술)"))
# D3 씨앗 초과 공급 — 같은 시드끼리 짝
common = [s for s in SEEDS if s in sd["racld"] and s in sd["rascld"]]
def e_of(df, seeds):
    a = stats(cell("df", "rascld", df, 0.01, seeds), "rascld", 0.01); b = stats(cell("df", "racld", df, 0.01, seeds), "racld", 0.01)
    return math.log(a["u"] / b["u"]) - math.log(unom(0.01, 6) / unom(0.01, 5))
e05, e5 = e_of(0.05, common), e_of(0.5, common)
bs = []
for _ in range(N_BOOT):
    ss = [common[RNG.randrange(len(common))] for _ in common]
    bs.append(e_of(0.05, ss) - e_of(0.5, ss))
lo3, hi3 = pct(bs, .025), pct(bs, .975)
v3 = "확산이 느릴수록 씨앗의 초과 공급이 크다" if lo3 > 0 else ("확산이 느릴수록 초과 공급이 작다" if hi3 < 0 else "차 검출 안 됨")
print("-- D3 씨앗 초과 공급 e(0.05) %+.3f · e(0.5) %+.3f · 차 %+.3f [%+.3f, %+.3f] (짝 시드 %d) → **%s**" % (e05, e5, e05 - e5, lo3, hi3, len(common), v3))
out["D1"] = {c: list(d1[c]) for c in CODES}; out["D2"] = r05; out["D3"] = [e05, e5, e05 - e5, lo3, hi3, v3]

print("\n-- E 수명 사다리 (밀도 8)  I 부풀림 · 순도 · 멸종")
sa = {(c, mu): alive_seeds("age", c, AGES, [mu], ASEEDS) for c in CODES for mu in AMUS}
ext = {}
for c in CODES:
    for mu in AMUS:
        for age in AGES:
            ext[(c, mu, age)] = sum(1 for s in ASEEDS if ("age", c, age, mu, s) in D and D[("age", c, age, mu, s)]["extinct_at"] >= 0)
boot_e = {}
e1 = {}
for c in CODES:
    for mu in AMUS:
        seeds = sa[(c, mu)]
        line = []
        for age in AGES:
            rows = cell("age", c, age, mu, seeds)
            st = stats(rows, c, mu)
            bs = [stats([rows[RNG.randrange(len(rows))] for _ in rows], c, mu) for _ in range(N_BOOT)]
            boot_e[(c, mu, age)] = (st, [b for b in bs if b])
            line.append("age %4d I %5.2f 순도 %.3f 멸종 %d" % (age, st["I"], st["pur"], ext[(c, mu, age)]))
            out["E"]["%s|%g|%d" % (c, mu, age)] = dict(st, extinct=ext[(c, mu, age)])
        (s1, b1), (s2, b2) = boot_e[(c, mu, 1200)], boot_e[(c, mu, 300)]
        n = min(len(b1), len(b2))
        v = [math.log(b1[i]["I"]) - math.log(b2[i]["I"]) for i in range(n)]
        e1[(c, mu)] = (math.log(s1["I"]) - math.log(s2["I"]), pct(v, .025), pct(v, .975))
        print("   %-7s μ %-5g 시드 %2d | " % (c, mu, len(seeds)) + " · ".join(line))
ok_e1 = all(e1[k][1] <= 0 <= e1[k][2] for k in e1)
print("-- E1 ln I(1200) − ln I(300): " + " · ".join("%s/μ%g %+.2f [%+.2f, %+.2f]" % (k[0], k[1], *e1[k]) for k in e1))
print("   → **%s**" % ("부풀림은 수명과 무관" if ok_e1 else "수명에 따라 달라지는 칸이 있다: " + ", ".join("%s/μ%g %s" % (k[0], k[1], "늘어남" if e1[k][1] > 0 else "줄어듦") for k in e1 if not (e1[k][1] <= 0 <= e1[k][2]))))
print("-- E2 멸종(μ 1%): " + " · ".join("%s %s" % (c, "/".join(str(ext[(c, 0.01, a)]) for a in AGES)) for c in CODES) + "  (300/600/1200)")
for mu in AMUS:
    seeds = [s for s in ASEEDS if s in sa[("racld", mu)] and s in sa[("rascld", mu)]]
    parts = []
    for age in AGES:
        a = stats(cell("age", "rascld", age, mu, seeds), "rascld", mu); b = stats(cell("age", "racld", age, mu, seeds), "racld", mu)
        parts.append("age %d 순도차 ln %+.2f · e %+.3f" % (age, math.log(a["pur"] / b["pur"]), math.log(a["u"] / b["u"]) - math.log(unom(mu, 6) / unom(mu, 5))))
    print("-- E3 서술 μ %g (짝 시드 %d): " % (mu, len(seeds)) + " · ".join(parts))
out["E1"] = {"%s|%g" % k: list(v) for k, v in e1.items()}; out["E2"] = {"%s|%g|%d" % k: v for k, v in ext.items()}
json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_lin10.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_lin10.json")
