# -*- coding: utf-8 -*-
"""F4′ — 🟡 사후 재분석(판정 아님). 설계검토-2026-09-17.md J1 · 해부7 §8a 정정 콜아웃.

해부 7 F4 는 씨앗의 침입 Δ 를 **명목 μ** 에 놓고 두 규칙 세계의 기울기를 비교했다. 그런데 재료 먼저 규칙은
세 코드 모두의 실현 오류율을 5~12배 낮추므로(해부 8 C: 원래 R 7~22 → 재료 먼저 R ≈ 1), '각 코드가 제 실현 공급만큼
깎인다' 는 **비특이 귀무**만으로도 명목 μ 기울기가 크게 준다. 여기서는 Δ 를 **실현 공급 u**(WT 부모 출생 중 변종 몫)에
놓고 두 세계를 한 곡선에 겹쳐, 귀무가 예측하는 몫과 그 위에 남는 **씨앗 특이 몫**을 가른다.

  자료  : 침입 Δ — inv6/(원래 · μ 0 · .001 · .002 · .003 · .004 · .005 · .0075 · .01) · inv7ff/(재료 먼저 · μ 0 · .001 · .003 · .005 · .01)
          같은 20 배경(짝) · Δ = s(팔) − WT 팔 11개 평균(dms_analyze.py 의 s)
          실현 공급 — lin7/orig · lin7/ff 단일 코드 접시 u_obs = Σ변종/Σ출생(WT 부모) · μ .001 · .005 · .01 · racld · rascld · rsacld
          racldx 는 lineage 기록이 없어 racld 의 부풀림 I 를 L = 6 에 적용(근사 — x 는 틱을 안 쓰므로 고리는 racld 와 같다)
  u(세계, 코드, μ) = u_nom(μ, L) × I(세계, 코드, μ) · I = u_obs/u_nom · ln I 를 μ 에 선형 보간(μ < .001 은 I(.001) 로 둠)
  귀무 곡선 f : 원래 세계의 (u, Δ) 8점 꺾은선(u = 0 포함) → 예측 Δ_ff_null(μ) = f(u_ff(μ))
  초과      : Δ_ff_obs(μ) − Δ_ff_null(μ) — 배경 짝 재추출 2,000회 95% (u 는 고정 — 시드 30개 합이라 흔들림이 작다)
  기울기    : μ {0 · .001 · .003 · .005 · .01} 다섯 점 최소제곱(해부 7 F4 와 같은 점) — 원래 · 재료 먼저 실측 · 귀무 예측
  R_obs = 1 − b_ff/b_orig · R_null = 1 − b_null/b_orig · 특이 몫 S = (b_null − b_ff)/b_null  [귀무가 남긴 기울기 중 실제로 사라진 몫]
  단위 공급당 깎임 : (Δ(u) − Δ(0))/u 를 u 가 비슷한 두 점에서 — 원래 μ .001 대 재료 먼저 μ .01
  한계: u_obs 는 단일 코드 접시 값(침입 공동체 안의 값이 아니다) · 우주선 돌연변이(코드 무관 · 두 세계 같음)는 u 에 안 넣었다
실행: python3 f4prime.py [petri 폴더]  →  _RESULT_f4prime.txt(표준 출력) · _RESULT_f4prime.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOTS = {"orig": "inv6", "ff": "inv7ff"}
MUS = {"orig": [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01], "ff": [0.0, 0.001, 0.003, 0.005, 0.01]}
MUS5 = [0.0, 0.001, 0.003, 0.005, 0.01]
LMUS = [0.001, 0.005, 0.01]
CODES = ["rascld", "rsacld", "racldx"]
UCODE = {"rascld": "rascld", "rsacld": "rsacld", "racldx": "racld"}
LEN = {"rascld": 6, "rsacld": 6, "racldx": 6}
N_BOOT = 2000
RNG = random.Random(20260920)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def mtag(mu):
    return "mu0" if mu == 0 else "mu" + str(mu).replace(".", "p")


def unom(mu, L):
    return 1 - (1 - mu * (2 / 3 + 10 / 11)) ** L


# ---------- 침입 Δ ----------
F = {}
for rule, root in ROOTS.items():
    for mu in MUS[rule]:
        for f in sorted(glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f))
            assert abs(z["mu_assay"] - mu) < 1e-12 and z["checksum_match"] and len(z["arms"]) == 14
            F[(rule, mu, z["seed"])] = z
bgs = sorted({k[2] for k in F if k[0] == "ff"})
assert all((r, mu, s) in F for r in ROOTS for mu in MUS[r] for s in bgs), "배경이 두 세계 · 모든 μ 에 다 있지 않다"
delta = {r: {c: {mu: [] for mu in MUS[r]} for c in CODES} for r in ROOTS}
for r in ROOTS:
    for mu in MUS[r]:
        for s in bgs:
            z = F[(r, mu, s)]
            wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
            base = sum(s_of(a) for a in wt) / len(wt)
            for a in z["arms"]:
                if a["kind"] == "mut":
                    delta[r][a["key"]][mu].append(s_of(a) - base)

# ---------- 실현 공급 ----------
uobs = {}
for rule in ROOTS:
    for code in ("racld", "rascld", "rsacld"):
        for mu in LMUS:
            b = m = 0
            for f in glob.glob(os.path.join(BASE, "lin7", rule, "lin_%s_mu%s_s*%s.json" % (code, str(mu).replace(".", "p"), "_ff" if rule == "ff" else ""))):
                z = json.load(open(f))
                if z["extinct_at"] >= 0:
                    continue
                b += z["wt_parent_births"]
                m += z["wt_parent_mut_births"]
            uobs[(rule, code, mu)] = m / b


def infl(rule, code, mu):
    L = len(code)
    pts = [(m, math.log(uobs[(rule, code, m)] / unom(m, L))) for m in LMUS]
    if mu <= pts[0][0]:
        return math.exp(pts[0][1])
    for (m0, y0), (m1, y1) in zip(pts, pts[1:]):
        if m0 <= mu <= m1:
            return math.exp(y0 + (y1 - y0) * (mu - m0) / (m1 - m0))
    return math.exp(pts[-1][1])


def u_of(rule, code, mu):
    return unom(mu, LEN[code]) * infl(rule, UCODE[code], mu)


def interp(xs, ys, x):
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0) if x1 > x0 else y0
    return ys[-1] if x > xs[-1] else ys[0]


def slope(mus, ys):
    mx = sum(mus) / len(mus)
    my = sum(ys) / len(ys)
    return sum((m - mx) * (y - my) for m, y in zip(mus, ys)) / sum((m - mx) ** 2 for m in mus)


IDX = [[RNG.randrange(len(bgs)) for _ in bgs] for _ in range(N_BOOT)]
out = {"bgs": bgs, "u": {}, "codes": {}}
print("== F4′ — Δ 를 실현 공급 u 에 놓고 두 규칙 세계를 겹친다 (🟡 사후 · 판정 아님)")
print("  배경 %d · 원래 μ %d · 재료 먼저 μ %d · u_obs: lin7 단일 접시(시드 701~730 합)" % (len(bgs), len(MUS["orig"]), len(MUS["ff"])))
print("\n-- 실현 공급 u = u_nom × I  (I = 부풀림, ln I 를 μ 에 보간)")
for code in CODES:
    row = []
    for rule in ROOTS:
        for mu in MUS[rule]:
            out["u"]["%s|%s|%g" % (rule, code, mu)] = u_of(rule, code, mu)
    print("   %-7s 원래   " % code + " · ".join("μ%g u %.4f(I %.2f)" % (mu, u_of("orig", code, mu), infl("orig", UCODE[code], mu)) for mu in MUS["orig"]))
    print("   %-7s 재료먼저" % "" + " · ".join("μ%g u %.4f(I %.2f)" % (mu, u_of("ff", code, mu), infl("ff", UCODE[code], mu)) for mu in MUS["ff"]))

for code in CODES:
    def calc(ix):
        mean = lambda r, mu: (sum(delta[r][code][mu][i] for i in ix) / len(ix)) if ix is not None else (sum(delta[r][code][mu]) / len(bgs))
        uo = [u_of("orig", code, mu) for mu in MUS["orig"]]
        do = [mean("orig", mu) for mu in MUS["orig"]]
        order = sorted(range(len(uo)), key=lambda i: uo[i])
        xs, ys = [uo[i] for i in order], [do[i] for i in order]
        dff = [mean("ff", mu) for mu in MUS5]
        dnull = [interp(xs, ys, u_of("ff", code, mu)) for mu in MUS5]
        dorig5 = [mean("orig", mu) for mu in MUS5]
        b_orig, b_ff, b_null = slope(MUS5, dorig5), slope(MUS5, dff), slope(MUS5, dnull)
        exc = [o - n for o, n in zip(dff, dnull)]
        # 단위 공급당 깎임: 원래 μ .001 대 재료 먼저 μ .01
        po = (mean("orig", 0.001) - mean("orig", 0.0)) / u_of("orig", code, 0.001)
        pf = (mean("ff", 0.01) - mean("ff", 0.0)) / u_of("ff", code, 0.01)
        return {"dff": dff, "dnull": dnull, "dorig5": dorig5, "exc": exc, "b_orig": b_orig, "b_ff": b_ff, "b_null": b_null,
                "R_obs": 1 - b_ff / b_orig if b_orig else float("nan"), "R_null": 1 - b_null / b_orig if b_orig else float("nan"),
                "S": (b_null - b_ff) / b_null if b_null else float("nan"), "per_u_orig": po, "per_u_ff": pf,
                "per_u_ratio": pf / po if po else float("nan")}
    pt = calc(None)
    bs = [calc(ix) for ix in IDX]
    ci = lambda key, j=None: (pct([(b[key][j] if j is not None else b[key]) for b in bs], .025), pct([(b[key][j] if j is not None else b[key]) for b in bs], .975))
    print("\n-- %s" % code)
    print("   %6s | %8s | %8s | %9s | %8s | %-22s" % ("μ", "u_ff", "Δ_orig", "Δ_ff 실측", "귀무 예측", "초과 = 실측 − 귀무 [95%]"))
    for j, mu in enumerate(MUS5):
        lo, hi = ci("exc", j)
        print("   %6g | %8.4f | %+8.3f | %+9.3f | %+8.3f | %+.3f [%+.3f, %+.3f]" % (mu, u_of("ff", code, mu), pt["dorig5"][j], pt["dff"][j], pt["dnull"][j], pt["exc"][j], lo, hi))
    r = {k: (pt[k],) + ci(k) for k in ("b_orig", "b_ff", "b_null", "R_obs", "R_null", "S", "per_u_ratio")}
    print("   기울기(명목 μ · 다섯 점): 원래 %+.1f [%+.1f, %+.1f] · 재료 먼저 실측 %+.1f [%+.1f, %+.1f] · 귀무 예측 %+.1f [%+.1f, %+.1f]"
          % (r["b_orig"] + r["b_ff"] + r["b_null"]))
    print("   줄어든 몫 R_obs %.2f [%.2f, %.2f] · 귀무 R_null %.2f [%.2f, %.2f] · **특이 몫 S = (귀무가 남긴 기울기 중 사라진 몫) %.2f [%.2f, %.2f]**"
          % (r["R_obs"] + r["R_null"] + r["S"]))
    print("   단위 공급당 깎임 (Δ(u) − Δ(0))/u: 원래(μ .001 · u %.4f) %+.2f · 재료 먼저(μ .01 · u %.4f) %+.2f · 비 %.2f [%.2f, %.2f]"
          % (u_of("orig", code, 0.001), pt["per_u_orig"], u_of("ff", code, 0.01), pt["per_u_ff"], *r["per_u_ratio"]))
    out["codes"][code] = dict(pt, ci={k: list(v[1:]) for k, v in r.items()}, exc_ci=[list(ci("exc", j)) for j in range(len(MUS5))])

print("\n읽기(🟡): 초과가 μ 0.5% · 1% 에서 하한 > 0 이면 '실현 공급이 같아도 재료 먼저 세계의 씨앗이 덜 깎인다 — 씨앗 특이 몫이 있다'.")
print("          S 는 그 몫의 크기(귀무가 남긴 기울기 중 실제로 사라진 몫). racldx 의 초과가 0 근처면 계산 틀이 헛것을 만들지 않는다는 대조.")
json.dump(out, open(os.path.join(BASE, "_RESULT_f4prime.json"), "w"), ensure_ascii=False, indent=1)
print("저장: _RESULT_f4prime.json")
