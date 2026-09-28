# -*- coding: utf-8 -*-
"""해부 5B 판정 — 촘촘한 μ 사다리. 사전등록 해부5-사전등록-2026-09-16.md §3 · §4. 결과 전에 썼다.

  mat 밀도 8 · racld 대 rascld · μ 0·0.1%·0.2%·0.3%·0.4%·0.5%·0.75%·1% · 시드 501~560 · 20,000틱
  층 1 (전 접시)          : D10 = ln(P(1만)+0.5) 짝 차이 · 멸종 접시는 P = 0. **해부 4A 와 견주는 참고값**
  층 2 (둘 다 산 접시만)  : 같은 D10 · n 을 함께 적는다. **판정은 이 층으로 한다**
     — 마른 실행에서 멸종 접시 하나가 층 1 의 평균을 +0.39 → +3.24 로 끌고 가는 것을 확인했다.
       성장 비교는 두 팔이 다 살아 있을 때만 정의되고, 멸종 자체는 아래 Wilson 구간으로 따로 판정한다.
  멸종                    : 팔별 비율 + Wilson 95% · 짝 차이
  교차 · 기울기           : **층 2** 의 D10 으로 — 기울기 하한 > 0 이면 '오류가 많을수록 racld 가 유리'
  🚨 두 층의 점추정 부호가 다르고 둘 중 하나라도 0 을 품지 않으면 그 μ 는 **보류**
실행: python3 mufine_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_mufine.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
CODES = ["racld", "rascld"]
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(501, 561))
SUB = os.environ.get("PETRI_SUB", "mufine")        # 마른 실행에서만 바꾼다
N_BOOT = 2000
RNG = random.Random(20260931)
IDX = [[RNG.randrange(len(SEEDS)) for _ in range(len(SEEDS))] for _ in range(N_BOOT)]


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot(vals, idx=None):
    if not vals:
        return (float("nan"),) * 3
    m = sum(vals) / len(vals)
    if idx is None:
        idx = [[RNG.randrange(len(vals)) for _ in range(len(vals))] for _ in range(N_BOOT)]
    bs = [sum(vals[i] for i in ix) / len(ix) for ix in idx]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r
    return None


def alive_at(z, t):
    return not (z["extinct_at"] >= 0 and z["extinct_at"] <= t)


def lnp(z, t):
    r = at(z, t)
    if not alive_at(z, t) or r is None:
        return math.log(0.5)
    return math.log(r[1] + 0.5)


def verdict(lo, hi):
    return "racld 가 더 잘 자란다" if lo > 0 else ("씨앗이 더 잘 자란다" if hi < 0 else "구별 안 됨")


# ---------- 읽기 · 점검 ----------
M = {}
for f in sorted(glob.glob(os.path.join(BASE, SUB, "mono_mat_*.json"))):
    z = json.load(open(f))
    M[(z["code"], z["mu"], z["seed"])] = z
missing = [(c, mu, s) for c in CODES for mu in MUS for s in SEEDS if (c, mu, s) not in M]
viol = sum(z["conservation_violations"] for z in M.values())
mu0_bad = [k for k, z in M.items() if k[1] == 0.0 and z["extinct_at"] < 0 and z["samples"][-1][3] != z["samples"][-1][1]]
print("== 해부 5B 촘촘한 μ 사다리 — 판정")
print("  파일 %d/%d · 재료 보존 위반 %d · 점검 ③ μ0 에서 출발 코드 아닌 개체가 있는 접시 %d"
      % (len(M), len(CODES) * len(MUS) * len(SEEDS), viol, len(mu0_bad)))
if missing or viol or mu0_bad:
    print("🚨 계측기/완결성 점검 실패 — 판정하지 않는다")
    for m in missing[:10]:
        print("    없음", m)
    sys.exit(1)

out = {"rows": [], "cross": None, "slope": None}
print("\n   %6s | %-28s | %-28s | %-24s | 판정" % ("μ", "층1 D(1만) 전 접시", "층2 D(1만) 둘 다 산 접시", "멸종 racld / 씨앗"))
held = []
for mu in MUS:
    d1 = [lnp(M[("racld", mu, s)], 10000) - lnp(M[("rascld", mu, s)], 10000) for s in SEEDS]
    pairs = [s for s in SEEDS if alive_at(M[("racld", mu, s)], 10000) and alive_at(M[("rascld", mu, s)], 10000)]
    d2 = [lnp(M[("racld", mu, s)], 10000) - lnp(M[("rascld", mu, s)], 10000) for s in pairs]
    e1 = [lnp(M[("racld", mu, s)], 20000) - lnp(M[("rascld", mu, s)], 20000) for s in SEEDS]
    m1, l1, h1 = boot(d1, IDX)
    m2, l2, h2 = boot(d2)
    n1, n2, n3 = boot(e1, IDX)
    ext = {c: sum(1 for s in SEEDS if M[(c, mu, s)]["extinct_at"] >= 0) for c in CODES}
    wa, wb = wilson(ext["racld"], len(SEEDS)), wilson(ext["rascld"], len(SEEDS))
    hold = (m1 * m2 < 0) and not (l1 <= 0 <= h1 and l2 <= 0 <= h2)
    if hold:
        held.append(mu)
    v = "🚨 보류 — 두 층의 방향이 다르다" if hold else verdict(l2, h2)
    print("   %5.2f%% | %+7.3f [%+.3f, %+.3f] | %+7.3f [%+.3f, %+.3f] n=%2d | %2d [%.2f-%.2f] / %2d [%.2f-%.2f] | %s"
          % (mu * 100, m1, l1, h1, m2, l2, h2, len(pairs),
             ext["racld"], wa[0], wa[1], ext["rascld"], wb[0], wb[1], v))
    out["rows"].append({"mu": mu, "d10_all": [m1, l1, h1], "d10_surv": [m2, l2, h2], "n_surv": len(pairs),
                        "d20_all": [n1, n2, n3], "ext": ext, "hold": hold, "verdict": v})

print("\n   %6s | %-28s | 판정(2만틱)" % ("μ", "D(2만틱) 전 접시"))
for r in out["rows"]:
    m, lo, hi = r["d20_all"]
    print("   %5.2f%% | %+7.3f [%+.3f, %+.3f] | %s" % (r["mu"] * 100, m, lo, hi, verdict(lo, hi)))

signs = [r["d10_surv"][0] for r in out["rows"]]          # 판정 층(둘 다 산 접시)
cross = None
for i in range(1, len(signs)):
    if signs[i] * signs[0] < 0:
        cross = MUS[i]
        break
# 기울기도 층 2 로 — 접시 수가 μ 마다 다르므로 μ 안에서 복원추출한다
xm = sum(MUS) / len(MUS)
sxx = sum((x - xm) ** 2 for x in MUS)
per = {}
for mu in MUS:
    pairs = [s for s in SEEDS if alive_at(M[("racld", mu, s)], 10000) and alive_at(M[("rascld", mu, s)], 10000)]
    per[mu] = [lnp(M[("racld", mu, s)], 10000) - lnp(M[("rascld", mu, s)], 10000) for s in pairs]
sl = []
for _ in range(N_BOOT):
    ys = []
    for mu in MUS:
        v = per[mu]
        ys.append(sum(v[RNG.randrange(len(v))] for _ in v) / len(v) if v else 0.0)
    ym = sum(ys) / len(ys)
    sl.append(sum((x - xm) * (y - ym) for x, y in zip(MUS, ys)) / sxx)
ym = sum(signs) / len(signs)
sm = sum((x - xm) * (y - ym) for x, y in zip(MUS, signs)) / sxx
slo, shi = pct(sl, 0.025), pct(sl, 0.975)
out["cross"] = cross
out["slope"] = [sm, slo, shi]
print("\n   ▸ (아래 교차·기울기는 판정 층인 층 2 로 잰다)")
print("   ▸ 부호가 바뀌는 가장 작은 μ: %s" % ("없다 — 0~1% 안에서 뒤집히지 않는다" if cross is None else "%.2f%%" % (cross * 100)))
print("   ▸ D(1만) 의 μ 기울기: %+.2f [%+.2f, %+.2f] → %s"
      % (sm, slo, shi, "오류가 많아질수록 racld 가 유리" if slo > 0
         else ("오류가 많아질수록 씨앗이 유리" if shi < 0 else "기울기 검출 안 됨")))
print("   ▸ 보류된 μ: %s" % ("없다" if not held else ", ".join("%.2f%%" % (m * 100) for m in held)))
print("\n   ▸ 미리 박은 대조 — 해부 2A 의 침입 대결은 μ 0 → 1% 에서 부호가 뒤집혔다(+0.39 → −0.75).")
print("     여기 성장 대결에서는 %s → **%s**"
      % ("교차가 없다" if cross is None else "μ %.2f%% 에서 교차한다" % (cross * 100),
         "그 뒤집힘은 경쟁 안에서만 일어나는 일이다" if cross is None else "성장만으로도 뒤집힌다"))

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_mufine.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_mufine.json"))
