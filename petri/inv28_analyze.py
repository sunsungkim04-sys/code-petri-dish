# -*- coding: utf-8 -*-
"""해부 28 판정 — 다시 굴림 확률 p 의 용량–반응. 사전등록 사전등록/해부28-사전등록-2026-10-02.md §4 · §5. 결과 전에 썼다.

  같은 mat 배경 20 · 같은 μ(0 · .001 · .003 · .005 · .01) · 같은 팔(ref + neutral 10 + rascld · rsacld · racldx) — 다섯 p 가 **짝**이다.
    p = 0     inv9mem/        기억하는 주사위(해부 9 기록 · 헛손질에서 다시 굴리지 않음)
    p = 0.25  inv28/p0p25/    헛손질마다 확률 0.25 로 기억을 지움 → 다음 시도가 새로 굴림
    p = 0.5   inv28/p0p5/
    p = 0.75  inv28/p0p75/
    p = 1     inv6/           원래 규칙(해부 6C 기록 · 헛손질마다 새로 굴림)
  Δ = s(팔) − WT 팔 11개 평균 · 기울기 b = μ 다섯 점 최소제곱(해부 9 와 같은 정의)
  배경 재추출 2,000회 95% — 다섯 p 가 같은 인덱스(짝) · 난수는 이름으로 연다

  G0 점검 — 하나라도 어긋나면 판정하지 않는다
    파일 완결성 · μ 일치 · 규칙 표시(redraw_p 값 / remember_die / 없음) · 배경 체크섬 MATCH · 충실도 · 팔 14 · 재료 보존
    회귀: inv28reg/p1 ≡ inv6 · inv28reg/none ≡ inv6 · inv28reg/p0 ≡ inv9mem (첫 배경 · μ 0.3% · 1% · 팔 14개 계열 전부)
          inv28reg/det_a ≡ inv28reg/det_b
    계측기 자기 검산: 양 끝 기울기를 이 코드로 다시 계산해 _RESULT_inv9.json 의 orig|rascld · mem|rascld 와 1e-9 안에서 같음

  판정선(§4 — 결과 전 고정)
    D1 단조     : 씨앗 b 의 점추정이 p 를 따라 엄격히 내려간다 b(0) > b(.25) > b(.5) > b(.75) > b(1)
                  그리고 인접 차 b(p_i) − b(p_{i+1}) 넷 중 셋 이상에서 짝 재추출 95% 하한 > 0
    D2 용량–반응: 재추출마다 다섯 점 b(p) 를 p 에 최소제곱 → 기울기 95% 상한 < 0
    종합        : ✅ = D1 ∧ D2 · ❌ = D2 상한 ≥ 0 · 🟡 = 그 밖(D2 는 서지만 D1 이 안 섬)
    D3 모양(서술 등급) : dev = max_{p∈.25,.5,.75} |b(p) − 직선 보간| / |b(1) − b(0)| — ≤ 0.25 '직선에 가깝다' · 아니면 '굽는다'(위 · 아래)
    D4 교차(서술)      : μ*(p) = 씨앗 Δ 평균 곡선의 첫 부호 변화(선형 보간) · 1% 안에 없으면 '1% 너머'
                         예측(§3): p 가 커질수록 μ* 가 줄거나 같다(점추정)
    D5 μ 0 대조        : 각 p 의 μ 0 씨앗 Δ − p 1 의 μ 0 씨앗 Δ 구간이 0 을 품는다(p 는 오류가 없으면 아무것도 안 바꿔야 한다)
    서술               : rsacld · racldx 기울기 · 끼운 계통이 사라진 팔 수
실행: python3 inv28_analyze.py [petri 폴더]  →  _RESULT_inv28.json · 화면 인쇄를 _RESULT_inv28.txt 로
마른 실행: PETRI_ROOT=dry28 PETRI_PS=0.5 PETRI_NBG=1 PETRI_MUS=0,0.01 python3 inv28_analyze.py   (회귀 · 결정론만 의미가 있다)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv28")
NEW_PS = [float(x) for x in os.environ.get("PETRI_PS", "0.25 0.5 0.75").replace(",", " ").split()]
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.003, 0.005, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
PS = [0.0] + NEW_PS + [1.0]
OUT = os.environ.get("PETRI_OUT", "_RESULT_inv28.json")


def ptag(p):
    return "p" + ("%g" % p).replace(".", "p")


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def root_of(p):
    if p == 0.0:
        return "inv9mem"
    if p == 1.0:
        return "inv6"
    return os.path.join(ROOT, ptag(p))


def pct(xs, q):
    xs = sorted(xs)
    k = (len(xs) - 1) * q
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def flag_ok(z, p):
    if p == 0.0:
        return bool(z.get("remember_die")) and "redraw_p" not in z
    if p == 1.0:
        return not z.get("remember_die") and "redraw_p" not in z and not z.get("find_first")
    return z.get("redraw_p") == p and not z.get("remember_die")


# ------------------------------------------------------------------ G0
F, bad = {}, []
for p in PS:
    for mu in MUS:
        for f in sorted(glob.glob(os.path.join(BASE, root_of(p), mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f))
            F[(p, mu, z["seed"])] = z
            if abs(z["mu_assay"] - mu) > 1e-12: bad.append("μ %g %g %d" % (p, mu, z["seed"]))
            if not flag_ok(z, p): bad.append("규칙 표시 p=%g μ=%g s=%d" % (p, mu, z["seed"]))
            if not z["checksum_match"]: bad.append("배경 체크섬 %g %g %d" % (p, mu, z["seed"]))
            if not z["fidelity"]["same"]: bad.append("충실도 %g %g %d" % (p, mu, z["seed"]))
            if len(z["arms"]) != 14: bad.append("팔 수 %g %g %d" % (p, mu, z["seed"]))
            if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %g %g %d" % (p, mu, z["seed"]))
bgs = sorted({k[2] for k in F if k[0] == NEW_PS[0]})[:NBG]
complete = len(bgs) == NBG and all((p, mu, s) in F for p in PS for mu in MUS for s in bgs)


def arms_sig(z):
    return json.dumps([[a["kind"], a["key"], a["series"]] for a in z["arms"]])


def reg_file(sub, mu, s):
    f = os.path.join(BASE, ROOT + "reg", sub, mtag(mu), "dms_mat_racld_s%05d_list.json" % s)
    return json.load(open(f)) if os.path.exists(f) else None


reg = []
if bgs:
    s0 = bgs[0]
    for mu in (0.003, 0.01):
        for sub, ref in (("p1", "inv6"), ("none", "inv6"), ("p0", "inv9mem")):
            a = reg_file(sub, mu, s0)
            rf = os.path.join(BASE, ref, mtag(mu), "dms_mat_racld_s%05d_list.json" % s0)
            b = json.load(open(rf)) if os.path.exists(rf) else None
            ok = a is not None and b is not None and arms_sig(a) == arms_sig(b)
            reg.append(("%s ≡ %s μ %g" % (sub, ref, mu), ok))
    a, b = reg_file("det_a", 0.003, s0), reg_file("det_b", 0.003, s0)
    reg.append(("det_a ≡ det_b μ 0.003", a is not None and b is not None and arms_sig(a) == arms_sig(b)))
for name, ok in reg:
    if not ok: bad.append("회귀 " + name)

print("== 해부 28 — 다시 굴림 확률 p 의 용량–반응 — 판정")
print("  p %s · 배경 %d · μ %d · 파일 %d/%d · 회귀 %d/%d · 점검 문제 %d" % (
    " · ".join("%g" % p for p in PS), len(bgs), len(MUS), sum(1 for k in F if k[2] in bgs),
    len(PS) * len(MUS) * NBG, sum(ok for _, ok in reg), len(reg), len(bad)))
for name, ok in reg:
    print("    회귀 %-26s %s" % (name, "✅" if ok else "❌"))
if not complete or bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:12]:
        print("    ", b)
    sys.exit(1)

# ------------------------------------------------------------------ Δ · 기울기
delta = {p: {c: {mu: [] for mu in MUS} for c in CODES} for p in PS}
lost = {p: 0 for p in PS}
for p in PS:
    for mu in MUS:
        for s in bgs:
            z = F[(p, mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                if a["kind"] == "mut":
                    delta[p][a["key"]][mu].append(s_of(a) - base)
                if a["series"][-1][2] == 0:
                    lost[p] += 1

RNG = random.Random("haebu28-redraw-dose")
IDX = [[RNG.randrange(NBG) for _ in range(NBG)] for _ in range(N_BOOT)]


def mean(p, c, mu, ix=None):
    v = delta[p][c][mu]
    return sum(v[i] for i in ix) / len(ix) if ix is not None else sum(v) / len(v)


def ols(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def b(p, c, ix=None):
    return ols(MUS, [mean(p, c, mu, ix) for mu in MUS])


def ci(v):
    return pct(v, .025), pct(v, .975)


out = {"ps": PS, "mus": MUS, "bgs": bgs, "lost": {"%g" % p: lost[p] for p in PS}, "delta": {}, "slopes": {}}

# 계측기 자기 검산 — 양 끝 기울기가 해부 9 판정 파일과 같은가
selfcheck = None
r9 = os.path.join(BASE, "_RESULT_inv9.json")
if os.path.exists(r9) and NBG == 20 and MUS == [0.0, 0.001, 0.003, 0.005, 0.01]:
    j9 = json.load(open(r9))["slopes"]
    d1 = abs(b(1.0, "rascld") - j9["orig|rascld"][0])
    d0 = abs(b(0.0, "rascld") - j9["mem|rascld"][0])
    selfcheck = max(d1, d0)
    print("  자기 검산: 양 끝 씨앗 기울기 대 _RESULT_inv9.json 최대 차 %.2e %s" % (selfcheck, "✅" if selfcheck < 1e-9 else "❌"))
    if selfcheck >= 1e-9:
        print("🚨 계측기 자기 검산 실패 — 판정하지 않는다")
        sys.exit(1)
out["selfcheck_maxdiff"] = selfcheck

print("\n-- 씨앗 Δ 표 (배경 평균 [95%])")
print("   %6s | " % "μ" + " | ".join("%-24s" % ("p = %g" % p) for p in PS))
for mu in MUS:
    cells = []
    for p in PS:
        m = mean(p, "rascld", mu)
        lo, hi = ci([mean(p, "rascld", mu, ix) for ix in IDX])
        cells.append("%+.3f [%+.3f, %+.3f]" % (m, lo, hi))
        out["delta"]["%g|rascld|%g" % (p, mu)] = [m, lo, hi]
    print("   %6g | %s" % (mu, " | ".join(cells)))

boot_b = {c: {p: [b(p, c, ix) for ix in IDX] for p in PS} for c in CODES}
print("\n-- 기울기 b (Δ / μ · 다섯 점)")
for c in CODES:
    row = []
    for p in PS:
        pt = b(p, c)
        lo, hi = ci(boot_b[c][p])
        out["slopes"]["%g|%s" % (p, c)] = [pt, lo, hi]
        row.append("p %g %+.1f [%+.1f, %+.1f]" % (p, pt, lo, hi))
    print("   %-7s " % c + " · ".join(row))

# ------------------------------------------------------------------ 판정
bs = [b(p, "rascld") for p in PS]
mono = all(bs[i] > bs[i + 1] for i in range(len(PS) - 1))
adj = []
for i in range(len(PS) - 1):
    diffs = [boot_b["rascld"][PS[i]][k] - boot_b["rascld"][PS[i + 1]][k] for k in range(N_BOOT)]
    lo, hi = ci(diffs)
    adj.append((PS[i], PS[i + 1], bs[i] - bs[i + 1], lo, hi))
n_adj = sum(1 for a in adj if a[3] > 0)
D1 = mono and n_adj >= len(adj) - 1
print("\n-- D1 단조: 점추정 엄격 감소 %s · 인접 차 하한 > 0 %d/%d → %s" % ("예" if mono else "아니오", n_adj, len(adj), "✅" if D1 else "❌"))
for a in adj:
    print("     b(%g) − b(%g) = %+.1f [%+.1f, %+.1f]" % a)

dose = [ols(PS, [boot_b["rascld"][p][k] for p in PS]) for k in range(N_BOOT)]
dose_pt = ols(PS, bs)
dlo, dhi = ci(dose)
D2 = dhi < 0
print("-- D2 용량–반응: b 를 p 에 최소제곱한 기울기 %+.1f [%+.1f, %+.1f] → %s" % (dose_pt, dlo, dhi, "✅" if D2 else "❌"))
verdict = "✅" if (D1 and D2) else ("❌" if not D2 else "🟡")
print("-- 종합 %s" % verdict)

span = bs[-1] - bs[0]
devs = []
for i, p in enumerate(PS[1:-1], start=1):
    lin = bs[0] + p * span
    devs.append((p, bs[i] - lin))
dev = max(abs(d) for _, d in devs) / abs(span)
shape = "직선에 가깝다" if dev <= 0.25 else ("굽는다 — 직선보다 덜 깎임(작은 p 에서 늦게)" if sum(d for _, d in devs) > 0 else "굽는다 — 직선보다 더 깎임(작은 p 에서 일찍)")
print("-- D3 모양: 최대 이탈 %.2f(|b(1) − b(0)| 대비) → %s · " % (dev, shape) + " · ".join("p %g %+.1f" % d for d in devs))


def crossing(p):
    ys = [mean(p, "rascld", mu) for mu in MUS]
    for i in range(len(MUS) - 1):
        if ys[i] > 0 >= ys[i + 1]:
            return MUS[i] + (MUS[i + 1] - MUS[i]) * ys[i] / (ys[i] - ys[i + 1])
    return None


cr = {p: crossing(p) for p in PS}
known = [cr[p] if cr[p] is not None else float("inf") for p in PS]
nonincr = all(known[i] >= known[i + 1] for i in range(len(PS) - 1))
print("-- D4 교차 μ*: " + " · ".join("p %g %s" % (p, ("%.3f%%" % (100 * cr[p])) if cr[p] is not None else "1% 너머") for p in PS)
      + " → p 를 따라 줄거나 같다: %s" % ("예" if nonincr else "아니오"))

d5 = []
for p in PS[:-1]:
    v = [mean(p, "rascld", 0.0, ix) - mean(1.0, "rascld", 0.0, ix) for ix in IDX]
    lo, hi = ci(v)
    d5.append((p, mean(p, "rascld", 0.0) - mean(1.0, "rascld", 0.0), lo, hi, lo <= 0 <= hi))
print("-- D5 μ 0 대조(p − p1): " + " · ".join("p %g %+.3f [%+.3f, %+.3f] %s" % (p, m, lo, hi, "✅" if ok else "❌") for p, m, lo, hi, ok in d5))
print("-- 서술: 끼운 계통이 사라진 팔 " + " · ".join("p %g %d" % (p, lost[p]) for p in PS))

out.update({"D1": {"mono": mono, "adj": adj, "n_adj_pos": n_adj, "pass": D1},
            "D2": {"pt": dose_pt, "lo": dlo, "hi": dhi, "pass": D2}, "verdict": verdict,
            "D3": {"dev": dev, "devs": devs, "shape": shape},
            "D4": {"crossing": {"%g" % p: cr[p] for p in PS}, "nonincreasing": nonincr},
            "D5": d5, "regression": reg})
json.dump(out, open(os.path.join(BASE, OUT), "w"), ensure_ascii=False, indent=1)
print("\n저장: %s" % OUT)
