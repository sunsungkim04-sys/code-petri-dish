# -*- coding: utf-8 -*-
"""해부 6C 판정 — 촘촘한 침입 사다리. 사전등록 해부6-사전등록-2026-09-16.md §4 · §5. 결과 전에 썼다.

  배경 mat 20 접시(해부 1 · 2A 와 같은 규칙) · 팔 = ref + neutral 10 + rascld · rsacld · racldx · --exact 1
  s = ln((m_T + 0.5)/m_0) − ln((r_T + 0.5)/r_0)  (dms_analyze.py 와 같다) · 기준 = WT 팔 11개 평균 · Δ = s − 기준
  배경 재추출 부트스트랩 2,000회 95%

  C0 회귀   : μ 0 팔 전부 = dms_main 기록 · μ 1% 팔 전부 = dms_mu1 기록(기록 앞 3칸) — 하나라도 다르면 중단
  C1 교차   : 씨앗 Δ 가 처음 '해로움'(상한 < 0)이 되는 μ · 점추정 부호가 바뀌는 자리
              기존 두 점(+0.388 · −0.753)의 직선 예측 = 0.34%
              → 처음 해로움이 0.1% 이하 : '이른 자리에서 뒤집힌다 — 담는 수(해부 5B)와 같은 자리'
              → 점추정 부호 변화가 [0.2%, 0.5%] : '직선 예측(0.34%) 근처에서 뒤집힌다'
              → 그 밖 : '둘 다 아님'
              곡률은 따로 적는다: Δ = a + bμ + cμ² 의 c — 하한 > 0 '앞쪽이 가파르고 뒤로 갈수록 완만'
              · 상한 < 0 '뒤로 갈수록 가파르다' · 그 밖 '직선과 구별 안 됨'
              (마른 실행에서 고침 — 가운데 결과 이름이 '오류율에 비례' 였는데 판정선은 교차 자리만 본다)
  C2 음성   : racldx 의 μ 기울기 구간이 0 을 품어야 한다(끼어든 x 는 틱을 안 쓴다) — 못 품으면 🚩 '모든 변종이 깎인다'
  C3 다시 시도 : 기울기(rsacld) − 기울기(rascld) 하한 > 0 → '같은 길이에서 느리게 다시 시도하는 코드는 덜 깎인다'
  C4 계통 순도 : 끝 기록에서 계통 안 '끼운 코드 그대로' 몫 — 씨앗 팔 대 WT 팔(ref) 대 rsacld 팔 (서술)
실행: python3 inv_analyze.py [petri 폴더]  →  _RESULT_inv.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv6")
MUS = [float(x) for x in os.environ["PETRI_MUS"].split(",")] if os.environ.get("PETRI_MUS") else [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
NBG = int(os.environ.get("PETRI_NBG", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
RNG = random.Random(20261003)
REF = {0.0: "dms_main", 0.01: "dms_mu1"}


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


F = {}
for mu in MUS:
    for f in sorted(glob.glob(os.path.join(BASE, ROOT, mtag(mu), "dms_mat_racld_s*_list.json"))):
        z = json.load(open(f))
        F[(mu, z["seed"])] = z
bgs = sorted({k[1] for k in F})
print("== 해부 6C 촘촘한 침입 사다리 — 판정")
bad = []
for (mu, s), z in F.items():
    if abs(z["mu_assay"] - mu) > 1e-12:
        bad.append("μ 어긋남 %s s%d" % (mu, s))
    if not z["checksum_match"]:
        bad.append("배경 체크섬 %s s%d" % (mu, s))
    if not z["fidelity"]["same"]:
        bad.append("복제 충실도 %s s%d" % (mu, s))
    if len(z["arms"]) != 14:
        bad.append("팔 수 %d %s s%d" % (len(z["arms"]), mu, s))
    if not all(a["cons_ok"] for a in z["arms"]):
        bad.append("재료 보존 %s s%d" % (mu, s))
complete = len(bgs) == NBG and all((mu, s) in F for mu in MUS for s in bgs)
# C0 회귀
reg_n, reg_same = 0, 0
for mu, folder in REF.items():
    if mu not in MUS:
        continue
    for s in bgs:
        old_f = os.path.join(BASE, folder, "dms_mat_racld_s%05d_all.json" % s)
        if not os.path.exists(old_f):
            bad.append("회귀 기록 없음 %s" % old_f)
            continue
        old = {(a["kind"], a["key"], a["post_seed"]): a["series"] for a in json.load(open(old_f))["arms"]}
        for a in F[(mu, s)]["arms"]:
            reg_n += 1
            o = old.get((a["kind"], a["key"], a["post_seed"]))
            reg_same += o is not None and o == [r[:3] for r in a["series"]]
print("  파일 %d/%d · 배경 %d · 점검 문제 %d · C0 회귀 %d/%d 팔 같음" % (len(F), len(MUS) * NBG, len(bgs), len(bad), reg_same, reg_n))
if not complete or bad or reg_same != reg_n or reg_n == 0:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in bad[:10]:
        print("    ", b)
    sys.exit(1)

IDX = [[RNG.randrange(len(bgs)) for _ in bgs] for _ in range(N_BOOT)]
delta = {c: {} for c in CODES}          # delta[c][mu] = [bg 별 Δ]
pure = {c: {} for c in CODES + ["ref"]}
for mu in MUS:
    for c in CODES + ["ref"]:
        pure[c][mu] = []
    for c in CODES:
        delta[c][mu] = []
    for s in bgs:
        z = F[(mu, s)]
        wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt) / len(wt)
        for a in z["arms"]:
            last = a["series"][-1]
            share = last[3] / last[2] if last[2] else float("nan")
            if a["kind"] == "mut":
                delta[a["key"]][mu].append(s_of(a) - base)
                pure[a["key"]][mu].append(share)
            elif a["kind"] == "ref":
                pure["ref"][mu].append(share)


def bmean(v):
    m = sum(v) / len(v)
    bs = [sum(v[i] for i in ix) / len(ix) for ix in IDX]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fit(xs, ys, deg):
    # 최소제곱 (deg 1 또는 2) — 정규방정식
    n = deg + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]; b[i], b[p] = b[p], b[i]
        for r in range(n):
            if r != i:
                f = A[r][i] / A[i][i]
                A[r] = [a - f * c for a, c in zip(A[r], A[i])]
                b[r] -= f * b[i]
    return [b[i] / A[i][i] for i in range(n)]


def curve(c, deg):
    pts = [fit(MUS, [sum(delta[c][mu][i] for i in ix) / len(ix) for mu in MUS], deg) for ix in IDX]
    pt = fit(MUS, [sum(delta[c][mu]) / len(bgs) for mu in MUS], deg)
    return [(pt[k], pct([p[k] for p in pts], 0.025), pct([p[k] for p in pts], 0.975)) for k in range(deg + 1)], pts


out = {"delta": {}, "C1": {}, "C2": {}, "C3": {}, "C4": {}}
print("\n   %6s | %-26s | %-26s | %-26s" % ("μ", "씨앗 rascld Δ", "rsacld Δ", "racldx Δ(음성)"))
for mu in MUS:
    cells = []
    for c in CODES:
        m, lo, hi = bmean(delta[c][mu])
        out["delta"].setdefault(c, {})[str(mu)] = [m, lo, hi]
        tag = "+" if lo > 0 else ("−" if hi < 0 else "0")
        cells.append("%+.3f [%+.3f,%+.3f] %s" % (m, lo, hi, tag))
    print("   %5.2f%% | %s" % (mu * 100, " | ".join(cells)))

# C1
seed = out["delta"]["rascld"]
first_harm = next((mu for mu in MUS if seed[str(mu)][2] < 0), None)
last_help = max((mu for mu in MUS if seed[str(mu)][1] > 0), default=None)
flip = None
for a, b in zip(MUS, MUS[1:]):
    if seed[str(a)][0] > 0 >= seed[str(b)][0]:
        pa, pb = seed[str(a)][0], seed[str(b)][0]
        flip = a + (b - a) * pa / (pa - pb)
        break
(q0, q1, q2), _ = curve("rascld", 2)
if first_harm is not None and first_harm <= 0.001:
    v1 = "이른 자리에서 뒤집힌다 — 담는 수(해부 5B)와 같은 자리"
elif flip is not None and 0.002 <= flip <= 0.005:
    v1 = "직선 예측(0.34%) 근처에서 뒤집힌다"
else:
    v1 = "둘 다 아님"
out["C1"] = {"first_harmful_mu": first_harm, "last_helpful_mu": last_help, "flip_interp": flip,
             "linear_prediction": 0.0034, "quad": [q0, q1, q2], "verdict": v1}
print("\n-- C1 씨앗의 교차")
print("   마지막 '이로움' μ %s · 처음 '해로움' μ %s · 점추정 부호 변화(보간) %s · 직선 예측 0.34%%"
      % ("-" if last_help is None else "%.2f%%" % (last_help * 100), "-" if first_harm is None else "%.2f%%" % (first_harm * 100),
         "-" if flip is None else "%.2f%%" % (flip * 100)))
print("   곡률 c = %+.0f [%+.0f, %+.0f] (Δ = a + bμ + cμ²) → %s" % (q2[0], q2[1], q2[2],
      "앞쪽이 가파르고 뒤로 갈수록 완만" if q2[1] > 0 else ("뒤로 갈수록 가파르다" if q2[2] < 0 else "직선과 구별 안 됨")))
print("   → **%s**" % v1)

# C2 · C3
lines = {}
for c in CODES:
    (a0, b1), pts = curve(c, 1)
    lines[c] = (b1, [p[1] for p in pts])
xb = lines["racldx"][0]
c2 = xb[1] <= 0 <= xb[2]
out["C2"] = {"slope": list(xb), "passes": c2}
print("\n-- C2 음성 대조 racldx 기울기 %+.1f [%+.1f, %+.1f] (Δ 단위 / μ) → %s"
      % (*xb, "통과 — 모든 변종이 깎이는 게 아니다" if c2 else "🚩 끼어든 x 도 깎인다 — 모든 변종이 깎일 수 있다"))
d = [a - b for a, b in zip(lines["rsacld"][1], lines["rascld"][1])]
dm = lines["rsacld"][0][0] - lines["rascld"][0][0]
dlo, dhi = pct(d, 0.025), pct(d, 0.975)
c3 = dlo > 0
out["C3"] = {"slope_rascld": list(lines["rascld"][0]), "slope_rsacld": list(lines["rsacld"][0]), "diff": [dm, dlo, dhi], "verdict": c3}
print("-- C3 기울기 씨앗 %+.1f [%+.1f, %+.1f] · rsacld %+.1f [%+.1f, %+.1f] · 차 %+.1f [%+.1f, %+.1f] → %s"
      % (*lines["rascld"][0], *lines["rsacld"][0], dm, dlo, dhi,
         "**같은 길이에서 느리게 다시 시도하는 코드는 덜 깎인다**" if c3 else "덜 깎인다고 할 수 없다"))

# C4
print("\n-- C4 끝 기록의 계통 순도(끼운 코드 그대로인 몫, 배경 평균) — 서술")
print("   %6s | %8s %8s %8s %8s" % ("μ", "WT(ref)", "씨앗", "rsacld", "racldx"))
for mu in MUS:
    row = []
    for c in ["ref"] + CODES:
        v = [x for x in pure[c][mu] if not math.isnan(x)]
        row.append(sum(v) / len(v) if v else float("nan"))
    out["C4"][str(mu)] = dict(zip(["ref"] + CODES, row))
    print("   %5.2f%% | %s" % (mu * 100, " ".join("%7.1f%%" % (100 * x) for x in row)))

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_inv.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_inv.json"))
