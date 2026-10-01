# -*- coding: utf-8 -*-
"""해부 21 판정 — §3.2 교차를 새 배경에서 재현. 사전등록 해부21-사전등록-2026-10-01.md §4 · §5. 결과 전에 썼다.

  계측기 · 팔 · Δ 정의는 해부 6C(inv_analyze.py)와 **같다**:
    s = ln((m_T + 0.5)/m_0) − ln((r_T + 0.5)/r_0) · 기준 = WT 팔 11개(ref + neutral 10) 평균 · Δ = s − 기준
    배경 재추출 부트스트랩 2,000회 95% · 교차 = 씨앗 Δ 점추정 평균 곡선의 첫 부호 변화(선형 보간)
  옛 값(해부 6C · inv6/ 배경 시드 3~22)을 같은 코드로 다시 계산해 _RESULT_inv.json 과 대조한다(계측기 자기 검산).

  G0 점검 — 하나라도 어긋나면 판정하지 않는다
    파일 완결성(후보 × μ) · μ 일치 · 복제 충실도 · 팔 14 · 재료 보존 · 새 배경은 main 기록 없음(main_checksum null)
    · 같은 배경의 μ 파일 8개가 grow_checksum · grow_top · grow_extinct 모두 같음(배경이 μ 와 무관하게 결정론적으로 같다)
    · 회귀 R: 시드 3 · μ 0.2% · 0.3% 의 팔 14개 궤적이 inv6/ 기록과 전부 같음 + 본실험 체크섬 MATCH
    · 옛 값 재계산이 _RESULT_inv.json 의 Δ(8 μ) · 교차와 1e-9 안에서 같음
    · 자격 배경(끝 우세 racld · 멸종 없음) ≥ 20 — 시드 순 앞 40 개까지 쓴다. 뺀 배경은 이유와 함께 인쇄한다

  판정선(§4 — 결과 전 고정)
    R1 반전       : 새 배경 씨앗 Δ(0) 하한 > 0 · Δ(1%) 상한 < 0
    R2 괄호       : Δ(0.1%) 하한 > 0 · Δ(0.3%) 상한 < 0   (원고 §3.2 "0.1% 에서 앞서고 0.3% 부터 뒤진다")
    R3 범주       : 해부 6C 의 C1 판정 규칙을 그대로 적용한 범주가 '둘 다 아님' 인가 (원고 "두 예측 사이")
    R4 교차 자리  : 교차(새) − 교차(옛) 의 부트스트랩 95% 구간(두 묶음 독립 재추출)이 0 을 품는다
                    · 재추출에서 교차가 정의되지 않는 몫이 어느 쪽이든 > 5% 면 R4 는 '판정 불가'
    종합          : ✅ 재현 = R1 ∧ R2 ∧ R4
                    ❌ 재현 실패 = ¬R1 ∨ 교차(새) 없음 ∨ (R4 ≠ ✅ ∧ 교차(새) ∉ [0.1%, 0.3%])   (검토 개정: R4 판정 불가도 포함)
                    🟡 부분 재현 = 그 밖(R1 은 서지만 R2 또는 R4 가 안 섬 · 교차는 [0.1%, 0.3%] 안, 또는 괄호 밖이지만 R4 ✅)
    부 판정(원고 §3.2 둘째 문단 — 같은 판정선)
      S1 = C2 racldx 기울기 구간이 0 을 품는다 · S2 = C3 기울기(rsacld) − 기울기(씨앗) 하한 > 0
      S3 = 곡률 c 하한 > 0(앞쪽이 가파르다) · S4 서술 = μ > 0 전부에서 계통 순도 WT > racldx > rsacld > 씨앗
    서술(판정 아님): 새 배경 앞 20 개만의 교차 · 옛 + 새 합친 교차 · 배경별 교차의 중앙값 · 사분위
실행: python3 inv21_analyze.py [petri 폴더]   환경: PETRI_ROOT(inv21) · PETRI_CANDS · PETRI_NMAX(40) · PETRI_NMIN(20) · PETRI_OUT
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv21")
OLD = os.environ.get("PETRI_OLD", "inv6")
CANDS = [int(x) for x in os.environ.get("PETRI_CANDS", " ".join(map(str, range(4201, 4251)))).split()]
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
REG_MUS = [0.002, 0.003]
NMAX = int(os.environ.get("PETRI_NMAX", "40"))
NMIN = int(os.environ.get("PETRI_NMIN", "20"))
CODES = ["rascld", "rsacld", "racldx"]
N_BOOT = 2000
UNDEF_MAX = 0.05
OUT = os.environ.get("PETRI_OUT", "_RESULT_inv21.json")
# 부트스트랩 난수는 묶음마다 이름으로 따로 연다(한 흐름을 나눠 쓰면 앞 단계 변경이 뒤 단계를 흔든다)
RNG_NEW = random.Random("haebu21-new-backgrounds")
RNG_OLD = random.Random("haebu21-old-inv6")
RNG_N20 = random.Random("haebu21-new-first20")
RNG_POOL = random.Random("haebu21-pooled")


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


def load(root):
    F = {}
    for mu in MUS:
        for f in sorted(glob.glob(os.path.join(BASE, root, mtag(mu), "dms_mat_racld_s*_list.json"))):
            z = json.load(open(f))
            F[(mu, z["seed"])] = z
    return F


def check_file(z, mu, s, bad, new):
    if abs(z["mu_assay"] - mu) > 1e-12:
        bad.append("μ 어긋남 %s s%d" % (mu, s))
    if not z["fidelity"]["same"]:
        bad.append("복제 충실도 %s s%d" % (mu, s))
    if z["grow_extinct"] < 0 and len(z["arms"]) != 14:
        bad.append("팔 수 %d %s s%d" % (len(z["arms"]), mu, s))
    if not all(a["cons_ok"] for a in z["arms"]):
        bad.append("재료 보존 %s s%d" % (mu, s))
    if new and z["main_checksum"] is not None:
        bad.append("새 배경인데 main 기록이 있다 s%d" % s)
    if not new and not z["checksum_match"]:
        bad.append("배경 체크섬 %s s%d" % (mu, s))


def deltas(F, bgs):
    """delta[c][mu] = 배경 순서대로 Δ 목록 · pure[c][mu] = 끝 계통 순도"""
    delta = {c: {mu: [] for mu in MUS} for c in CODES}
    pure = {c: {mu: [] for mu in MUS} for c in CODES + ["ref"]}
    for mu in MUS:
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
    return delta, pure


def flip_of(curve):
    """curve = MUS 순서의 평균 Δ. 첫 부호 변화(+ → ≤0)의 선형 보간. 없으면 None"""
    for (a, pa), (b, pb) in zip(zip(MUS, curve), zip(MUS[1:], curve[1:])):
        if pa > 0 >= pb:
            return a + (b - a) * pa / (pa - pb)
    return None


def fit(xs, ys, deg):
    n = deg + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]
        b[i], b[p] = b[p], b[i]
        for r in range(n):
            if r != i:
                f = A[r][i] / A[i][i]
                A[r] = [a - f * c for a, c in zip(A[r], A[i])]
                b[r] -= f * b[i]
    return [b[i] / A[i][i] for i in range(n)]


class Set:
    """한 배경 묶음의 Δ · 부트스트랩 · 교차"""

    def __init__(self, name, delta, pure, rng):
        self.name, self.delta, self.pure = name, delta, pure
        self.n = len(delta["rascld"][0.0])
        self.idx = [[rng.randrange(self.n) for _ in range(self.n)] for _ in range(N_BOOT)]

    def mean(self, c, mu, ix=None):
        v = self.delta[c][mu]
        ix = range(self.n) if ix is None else ix
        return sum(v[i] for i in ix) / len(ix)

    def bmean(self, c, mu):
        bs = [self.mean(c, mu, ix) for ix in self.idx]
        return self.mean(c, mu), pct(bs, 0.025), pct(bs, 0.975)

    def flips(self):
        pt = flip_of([self.mean("rascld", mu) for mu in MUS])
        bs = [flip_of([self.mean("rascld", mu, ix) for mu in MUS]) for ix in self.idx]
        return pt, bs

    def curve(self, c, deg):
        pt = fit(MUS, [self.mean(c, mu) for mu in MUS], deg)
        pts = [fit(MUS, [self.mean(c, mu, ix) for mu in MUS], deg) for ix in self.idx]
        return [(pt[k], pct([p[k] for p in pts], 0.025), pct([p[k] for p in pts], 0.975)) for k in range(deg + 1)], pts

    def per_bg_flips(self):
        return [flip_of([self.delta["rascld"][mu][i] for mu in MUS]) for i in range(self.n)]


def c1_verdict(S):
    """해부 6C C1 규칙 그대로"""
    seed = {mu: S.bmean("rascld", mu) for mu in MUS}
    first_harm = next((mu for mu in MUS if seed[mu][2] < 0), None)
    flip = flip_of([seed[mu][0] for mu in MUS])
    if first_harm is not None and first_harm <= 0.001:
        return "이른 자리에서 뒤집힌다 — 담는 수(해부 5B)와 같은 자리", first_harm, flip
    if flip is not None and 0.002 <= flip <= 0.005:
        return "직선 예측(0.34%) 근처에서 뒤집힌다", first_harm, flip
    return "둘 다 아님", first_harm, flip


def fmt_mu(x):
    return "-" if x is None else "%.3f%%" % (100 * x)


def main():
    out = {"prereg": "해부21-사전등록-2026-10-01.md", "root": ROOT, "candidates": CANDS}
    bad = []
    print("== 해부 21 — §3.2 교차를 새 배경에서 재현 · 판정")
    # ---- 새 배경
    F = load(ROOT)
    missing = [(mu, s) for mu in MUS for s in CANDS if (mu, s) not in F]
    extra = sorted({k[1] for k in F} - set(CANDS))
    if missing:
        bad.append("파일 빠짐 %d개(예: %s)" % (len(missing), missing[:3]))
    if extra:
        bad.append("후보 밖 시드 파일 %s" % extra[:5])
    for (mu, s), z in F.items():
        check_file(z, mu, s, bad, new=True)
    # 같은 배경의 μ 파일이 같은 배경을 키웠나
    bgfacts = {}
    for s in CANDS:
        facts = {(F[(mu, s)]["grow_checksum"], json.dumps(F[(mu, s)]["grow_top"][:1]), F[(mu, s)]["grow_extinct"]) for mu in MUS if (mu, s) in F}
        if len(facts) > 1:
            bad.append("시드 %d: μ 파일들의 배경이 다르다(grow_checksum/top/extinct %d 가지)" % (s, len(facts)))
        z = next((F[(mu, s)] for mu in MUS if (mu, s) in F), None)
        if z is not None:
            bgfacts[s] = (z["grow_extinct"], z["grow_top"][0][0] if z["grow_top"] else None, z["grow_top"][0][1] if z["grow_top"] else 0, z["grow_pop"])
    qual = [s for s in CANDS if s in bgfacts and bgfacts[s][0] < 0 and bgfacts[s][1] == "racld"]
    dropped = [(s, bgfacts.get(s)) for s in CANDS if s not in qual]
    used = qual[:NMAX]
    print("  후보 %d · 파일 %d/%d · 자격(끝 우세 racld · 멸종 없음) %d · 판정에 씀 %d(시드 순 앞 %d 까지)"
          % (len(CANDS), len(F), len(CANDS) * len(MUS), len(qual), len(used), NMAX))
    print("  뺀 배경 %d:" % len(dropped))
    for s, f in dropped:
        print("     s%d  %s" % (s, "파일 없음" if f is None else ("멸종 틱 %d" % f[0] if f[0] >= 0 else "끝 우세 %s (%d/%d)" % (f[1], f[2], f[3]))))
    if len(qual) > NMAX:
        print("  자격이지만 상한 밖이라 안 씀: %s" % qual[NMAX:])
    out["backgrounds_used"] = used
    out["backgrounds_dropped"] = [[s, None if f is None else list(f)] for s, f in dropped]
    # ---- 회귀 R
    reg_n, reg_same = 0, 0
    for mu in REG_MUS:
        fs = glob.glob(os.path.join(BASE, ROOT + "reg", mtag(mu), "dms_mat_racld_s00003_list.json"))
        of = os.path.join(BASE, OLD, mtag(mu), "dms_mat_racld_s00003_list.json")
        if not fs or not os.path.exists(of):
            bad.append("회귀 파일 없음 μ %s" % mu)
            continue
        z, o = json.load(open(fs[0])), json.load(open(of))
        check_file(z, mu, 3, bad, new=False)
        old = {(a["kind"], a["key"], a["post_seed"]): a["series"] for a in o["arms"]}
        for a in z["arms"]:
            reg_n += 1
            reg_same += old.get((a["kind"], a["key"], a["post_seed"])) == a["series"]
        if len(o["arms"]) != len(z["arms"]):
            bad.append("회귀 팔 수 다름 μ %s: %d 대 %d" % (mu, len(z["arms"]), len(o["arms"])))
    print("  회귀 R(시드 3 · μ 0.2 · 0.3%%): 팔 %d/%d 궤적 같음 · 본실험 체크섬 대조 포함" % (reg_same, reg_n))
    if reg_n == 0 or reg_same != reg_n:
        bad.append("회귀 R 실패 %d/%d" % (reg_same, reg_n))
    # ---- 옛 값 재계산(해부 6C)
    FO = load(OLD)
    old_bgs = sorted({k[1] for k in FO})
    old_ok = len(old_bgs) == 20 and all((mu, s) in FO for mu in MUS for s in old_bgs)
    if not old_ok:
        bad.append("옛 기록 inv6 불완전(배경 %d)" % len(old_bgs))
    ref_json = json.load(open(os.path.join(BASE, "_RESULT_inv.json")))
    dO, pO = deltas(FO, old_bgs) if old_ok else ({}, {})
    SO = Set("옛(해부 6C)", dO, pO, RNG_OLD) if old_ok else None
    if SO:
        dev = max(abs(SO.mean("rascld", mu) - ref_json["delta"]["rascld"][str(mu)][0]) for mu in MUS)
        fo_pt = flip_of([SO.mean("rascld", mu) for mu in MUS])
        dev = max(dev, abs(fo_pt - ref_json["C1"]["flip_interp"]))
        print("  옛 값 재계산: 배경 %d · Δ · 교차가 _RESULT_inv.json 과 최대 %.1e 차 (교차 %s)" % (len(old_bgs), dev, fmt_mu(fo_pt)))
        if dev > 1e-9:
            bad.append("옛 값 재계산이 정본과 다르다(%.2e)" % dev)
    print("  점검 문제 %d" % len(bad))
    if bad or len(used) < NMIN:
        print("🚨 점검 실패 또는 자격 배경 부족(%d < %d) — 판정하지 않는다" % (len(used), NMIN) if len(used) < NMIN else "🚨 점검 실패 — 판정하지 않는다")
        for b in bad[:15]:
            print("    ", b)
        out["halted"] = bad[:50] + ([] if len(used) >= NMIN else ["자격 배경 %d < %d" % (len(used), NMIN)])
        json.dump(out, open(os.path.join(BASE, OUT), "w"), ensure_ascii=False, indent=1)
        sys.exit(1)

    # ---- 새 배경 판정
    dN, pN = deltas(F, used)
    SN = Set("새", dN, pN, RNG_NEW)
    out["delta_new"] = {c: {str(mu): list(SN.bmean(c, mu)) for mu in MUS} for c in CODES}
    out["delta_old"] = {c: {str(mu): list(SO.bmean(c, mu)) for mu in MUS} for c in CODES}
    print("\n   %6s | %-28s | %-28s || %-26s" % ("μ", "씨앗 Δ 새(n=%d)" % SN.n, "씨앗 Δ 옛(n=%d)" % SO.n, "rsacld · racldx Δ 새"))
    for mu in MUS:
        a, b = SN.bmean("rascld", mu), SO.bmean("rascld", mu)
        r, x = SN.bmean("rsacld", mu), SN.bmean("racldx", mu)
        print("   %5.2f%% | %+.3f [%+.3f, %+.3f] %s | %+.3f [%+.3f, %+.3f] %s || %+.3f · %+.3f"
              % (mu * 100, *a, "+" if a[1] > 0 else ("−" if a[2] < 0 else "0"), *b, "+" if b[1] > 0 else ("−" if b[2] < 0 else "0"), r[0], x[0]))

    d0, d1 = SN.bmean("rascld", 0.0), SN.bmean("rascld", 0.01)
    R1 = d0[1] > 0 and d1[2] < 0
    e1, e3 = SN.bmean("rascld", 0.001), SN.bmean("rascld", 0.003)
    R2 = e1[1] > 0 and e3[2] < 0
    v3, fh, fl3 = c1_verdict(SN)
    R3 = v3 == "둘 다 아님"
    fN, bN = SN.flips()
    fO, bO = SO.flips()
    undefN = sum(x is None for x in bN) / N_BOOT
    undefO = sum(x is None for x in bO) / N_BOOT
    diffs = [a - b for a, b in zip(bN, bO) if a is not None and b is not None]
    dlo, dhi = (pct(diffs, 0.025), pct(diffs, 0.975)) if diffs else (float("nan"), float("nan"))
    nlo, nhi = (pct([x for x in bN if x is not None], 0.025), pct([x for x in bN if x is not None], 0.975)) if undefN < 1 else (float("nan"),) * 2
    olo, ohi = (pct([x for x in bO if x is not None], 0.025), pct([x for x in bO if x is not None], 0.975)) if undefO < 1 else (float("nan"),) * 2
    if undefN > UNDEF_MAX or undefO > UNDEF_MAX or fN is None:
        R4 = None
    else:
        R4 = dlo <= 0 <= dhi
    print("\n-- R1 반전: Δ(0) %+.3f [%+.3f, %+.3f] · Δ(1%%) %+.3f [%+.3f, %+.3f] → %s" % (*d0, *d1, "✅" if R1 else "❌"))
    print("-- R2 괄호: Δ(0.1%%) %+.3f [%+.3f, %+.3f] · Δ(0.3%%) %+.3f [%+.3f, %+.3f] → %s" % (*e1, *e3, "✅" if R2 else "❌"))
    print("-- R3 범주(해부 6C C1 규칙): 처음 해로움 %s · 교차 %s → '%s' → %s" % (fmt_mu(fh), fmt_mu(fl3), v3, "✅ 같은 범주" if R3 else "❌ 범주가 바뀜"))
    print("-- R4 교차 자리: 새 %s [%s, %s] (정의 안 됨 %.1f%%) · 옛 %s [%s, %s] (정의 안 됨 %.1f%%)"
          % (fmt_mu(fN), fmt_mu(nlo), fmt_mu(nhi), 100 * undefN, fmt_mu(fO), fmt_mu(olo), fmt_mu(ohi), 100 * undefO))
    print("   새 − 옛 = %s [%s, %s] → %s" % ("-" if fN is None else "%+.3f%%p" % (100 * (fN - fO)), "%+.3f" % (100 * dlo), "%+.3f" % (100 * dhi),
          "판정 불가" if R4 is None else ("✅ 0 을 품는다" if R4 else "❌ 옮겨졌다")))
    inside = fN is not None and 0.001 <= fN <= 0.003
    # 검토 개정(10-01 · 결과 전): R4 가 '판정 불가'(None)이고 교차가 괄호 밖이어도 ❌ — 옛 자리와 같다는 근거 없이 괄호 밖이면 실패
    if not R1 or fN is None or (R4 is not True and not inside):
        verdict = "❌ 재현 실패"
    elif R1 and R2 and R4:
        verdict = "✅ 재현"
    else:
        verdict = "🟡 부분 재현"
    print("   → **%s**" % verdict)

    # ---- 부 판정
    lines = {}
    for c in CODES:
        (a0, b1), pts = SN.curve(c, 1)
        lines[c] = (b1, [p[1] for p in pts])
    xb = lines["racldx"][0]
    S1 = xb[1] <= 0 <= xb[2]
    dd = [a - b for a, b in zip(lines["rsacld"][1], lines["rascld"][1])]
    dm = lines["rsacld"][0][0] - lines["rascld"][0][0]
    S2lo, S2hi = pct(dd, 0.025), pct(dd, 0.975)
    S2 = S2lo > 0
    (q0, q1, q2), _ = SN.curve("rascld", 2)
    S3 = q2[1] > 0
    print("\n-- 부 판정 (원고 §3.2 둘째 문단)")
    print("   S1 racldx 기울기 %+.1f [%+.1f, %+.1f] → %s (옛 +0.7)" % (*xb, "✅ 0 을 품는다" if S1 else "❌"))
    print("   S2 기울기 씨앗 %+.1f [%+.1f, %+.1f] · rsacld %+.1f [%+.1f, %+.1f] · 차 %+.1f [%+.1f, %+.1f] → %s (옛 −110 · −28 · +83)"
          % (*lines["rascld"][0], *lines["rsacld"][0], dm, S2lo, S2hi, "✅" if S2 else "❌"))
    print("   S3 곡률 c %+.0f [%+.0f, %+.0f] → %s (옛 +10,198)" % (*q2, "✅ 앞쪽이 가파르다" if S3 else "❌"))
    order_ok = []
    print("   S4 서술 — 끝 계통 순도(배경 평균) WT · 씨앗 · rsacld · racldx")
    purity = {}
    for mu in MUS:
        row = {}
        for c in ["ref"] + CODES:
            v = [x for x in pN[c][mu] if not math.isnan(x)]
            row[c] = sum(v) / len(v) if v else float("nan")
        purity[str(mu)] = row
        ok = row["ref"] > row["racldx"] > row["rsacld"] > row["rascld"]
        if mu > 0:
            order_ok.append(ok)
        print("   %5.2f%% | %5.1f %5.1f %5.1f %5.1f %s" % (mu * 100, *(100 * row[c] for c in ["ref", "rascld", "rsacld", "racldx"]), "" if mu == 0 else ("순서 같음" if ok else "순서 다름")))
    S4 = all(order_ok)

    # ---- 서술
    S20 = Set("새 앞 20", {c: {mu: v[:20] for mu, v in dN[c].items()} for c in CODES}, None, RNG_N20)
    f20, _ = S20.flips()
    pooled = {c: {mu: dN[c][mu] + dO[c][mu] for mu in MUS} for c in CODES}
    SP = Set("옛 + 새", pooled, None, RNG_POOL)
    fP, bP = SP.flips()
    bPd = [x for x in bP if x is not None]
    pbN = [x for x in SN.per_bg_flips()]
    pbO = [x for x in SO.per_bg_flips()]

    def q3(v):
        w = [x for x in v if x is not None]
        return (pct(w, 0.25), pct(w, 0.5), pct(w, 0.75), len(v) - len(w)) if w else (float("nan"),) * 3 + (len(v),)
    print("\n-- 서술(판정 아님)")
    print("   새 앞 20 배경만의 교차 %s · 옛 + 새(n=%d) 교차 %s [%s, %s]" % (fmt_mu(f20), SP.n, fmt_mu(fP), fmt_mu(pct(bPd, 0.025)), fmt_mu(pct(bPd, 0.975))))
    a, b = q3(pbN), q3(pbO)
    print("   배경별 교차 중앙값 [사분위]: 새 %s [%s, %s] (교차 없는 배경 %d) · 옛 %s [%s, %s] (%d)"
          % (fmt_mu(a[1]), fmt_mu(a[0]), fmt_mu(a[2]), a[3], fmt_mu(b[1]), fmt_mu(b[0]), fmt_mu(b[2]), b[3]))

    out.update({
        "n_new": SN.n, "n_old": SO.n,
        "R1": {"d0": list(d0), "d1pct": list(d1), "pass": R1},
        "R2": {"d0p1": list(e1), "d0p3": list(e3), "pass": R2},
        "R3": {"category": v3, "first_harmful": fh, "flip": fl3, "pass": R3},
        "R4": {"flip_new": fN, "flip_new_ci": [nlo, nhi], "undef_new": undefN, "flip_old": fO, "flip_old_ci": [olo, ohi], "undef_old": undefO,
               "diff": None if fN is None else fN - fO, "diff_ci": [dlo, dhi], "pass": R4},
        "verdict": verdict,
        "S1": {"slope_racldx": list(xb), "pass": S1},
        "S2": {"slope_rascld": list(lines["rascld"][0]), "slope_rsacld": list(lines["rsacld"][0]), "diff": [dm, S2lo, S2hi], "pass": S2},
        "S3": {"quad": [list(q0), list(q1), list(q2)], "pass": S3},
        "S4": {"purity": purity, "order_same_all_mu": S4},
        "desc": {"flip_new_first20": f20, "flip_pooled": fP, "flip_pooled_ci": [pct(bPd, 0.025), pct(bPd, 0.975)],
                 "per_bg_flip_new_q": list(a), "per_bg_flip_old_q": list(b)},
    })
    json.dump(out, open(os.path.join(BASE, OUT), "w"), ensure_ascii=False, indent=1)
    print("\n저장: " + OUT)


if __name__ == "__main__":
    main()
