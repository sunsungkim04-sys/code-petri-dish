# -*- coding: utf-8 -*-
"""해부 5A 판정 — 재료 밀도 사다리. 사전등록 해부5-사전등록-2026-09-16.md §2 · §4. 결과 전에 썼다.

  칸 = mat 밀도 4·8·14·20·32 + space · 코드 = acld(4) racld(5) rascld(6) rsacld(6) rascled(7)
  μ 0 · 시드 501~530 · 20,000틱 · 시드 짝 부트스트랩 2,000회 95%

  A1 racld 대 rascld : D10 = ln(P(1만)+0.5) 짝 차이 · D20 = 2만틱
     재료 가설 = D10 이 밀도 4→32 로 단조 감소 · space 에서 0 을 품는다 · D20 은 모든 칸에서 0
  A2 길이 사다리     : 칸마다 ln P(1만) 을 코드 길이에 회귀한 기울기. **두 벌** — 네 길이(4·5·6·7) 와
     acld 뺀 세 길이(5·6·7). acld 는 다른 셋보다 개체가 수십 배 적어 회귀를 혼자 끌 수 있다(마른 실행에서 확인).
     재료 가설 = 빡빡한 칸(4·8)에서 **둘 중 하나라도** 기울기 음수 · 넉넉한 칸(32·space)에서 0
  A3 음성 대조       : rascld 대 rsacld(둘 다 6글자) — '밀도 8 → space' 로 D 가 움직인 폭을 A1 과 견준다(차이의 차이).
     밀도 4 는 개체가 수백이라 구간이 ±2 로 벌어져 대조에서 뺀다(마른 실행에서 확인 · 표에는 그대로 적는다).
  직접 측정          : 끝 자유 재료 비율(<20% 바닥 · >50% 남아돎) · 점유율(>90% 자리 참) · 개체당 출산율
실행: python3 dens_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_dens.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
RUNGS = [("mat_d04", 4, "밀도 4"), ("mat_d08", 8, "밀도 8"), ("mat_d14", 14, "밀도 14"),
         ("mat_d20", 20, "밀도 20"), ("mat_d32", 32, "밀도 32"), ("space", None, "space(재료없음)")]
DENS_RUNGS = [r for r in RUNGS if r[1] is not None]
CODES = ["acld", "racld", "rascld", "rsacld", "rascled"]
LADDER = ["acld", "racld", "rascld", "rascled"]        # 길이 4·5·6·7
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(501, 531))
SUB = os.environ.get("PETRI_SUB", "dens")          # 마른 실행에서만 바꾼다
SUB2 = os.environ.get("PETRI_SUB2", "dens_twice")
DISH = 10732
N_BOOT = 2000
RNG = random.Random(20260930)
IDX = [[RNG.randrange(len(SEEDS)) for _ in range(len(SEEDS))] for _ in range(N_BOOT)]


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot(vals):
    """짝 차이 리스트 → (점추정, 하한, 상한)"""
    m = sum(vals) / len(vals)
    bs = [sum(vals[i] for i in ix) / len(ix) for ix in IDX]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def verdict(lo, hi, pos, neg):
    return pos if lo > 0 else (neg if hi < 0 else "구별 안 됨")


def at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r
    return None


def lnp(z, t):
    r = at(z, t)
    if z["extinct_at"] >= 0 and (r is None or z["extinct_at"] <= t):
        return math.log(0.5)
    if r is None:
        r = z["samples"][-1]
    return math.log(r[1] + 0.5)


# ---------- 읽기 · 점검 ----------
D = {}
for folder, dens, _ in RUNGS:
    for f in sorted(glob.glob(os.path.join(BASE, SUB, folder, "mono_*_mu0_s*.json"))):
        z = json.load(open(f))
        D[(folder, z["code"], z["seed"])] = z

missing = [(fo, c, s) for fo, _, _ in RUNGS for c in CODES for s in SEEDS if (fo, c, s) not in D]
viol = sum(z["conservation_violations"] for z in D.values())
mu0_bad = [k for k, z in D.items() if z["extinct_at"] < 0 and z["samples"][-1][3] != z["samples"][-1][1]]
space_bad = [k for k, z in D.items() if (k[0] == "space") != (z["free_end"] is None)]
twice_bad, twice_n = [], 0
for f in sorted(glob.glob(os.path.join(BASE, SUB2, "*.json"))):
    z = json.load(open(f))
    o = D.get(("mat_d08", z["code"], z["seed"]))
    if o:
        twice_n += 1
        if o["checksum"] != z["checksum"]:
            twice_bad.append(z["code"])

print("== 해부 5A 재료 밀도 사다리 — 판정")
print("  파일 %d/%d · 재료 보존 위반 %d · 점검 ③ μ0 인데 출발 코드 아닌 개체 %d접시 · 점검 ④ space 재료기록 어긋남 %d · 점검 ⑤ 같은 시드 두 번 %d/%d 일치"
      % (len(D), len(RUNGS) * len(CODES) * len(SEEDS), viol, len(mu0_bad), len(space_bad), twice_n - len(twice_bad), twice_n))
print("  점검 ① mono.js v1.1.0 회귀 검사는 발사 전 따로 통과(해부 5 노트 §4)")
if missing or viol or mu0_bad or space_bad or twice_bad or twice_n == 0:
    print("🚨 계측기/완결성 점검 실패 — 판정하지 않는다")
    for m in missing[:10]:
        print("    없음", m)
    sys.exit(1)

out = {"rungs": {}, "A1": {}, "A2": {}, "A3": {}}

# ---------- 칸별 직접 측정 ----------
print("\n-- 무엇이 제약인가 (칸마다 · racld/rascld 평균 · 2만틱)")
print("   %-16s %9s %9s %11s %11s" % ("칸", "끝 개체", "점유%", "자유재료%", "개체당출산/천틱"))
for folder, dens, name in RUNGS:
    pops, occ, fr, br = [], [], [], []
    for c in ("racld", "rascld"):
        for s in SEEDS:
            z = D[(folder, c, s)]
            last = z["samples"][-1]
            a = at(z, 10000)
            pops.append(last[1]); occ.append(100.0 * last[1] / DISH)
            if z["free_end"] is not None and z["mat_total0"]:
                fr.append(100.0 * z["free_end"] / z["mat_total0"])
            mid = [r[1] for r in z["samples"] if 10000 <= r[0] <= 20000]
            if a and last[0] > a[0] and mid and sum(mid) > 0:
                br.append(1000.0 * (last[4] - a[4]) / ((last[0] - a[0]) * (sum(mid) / len(mid))))
    p = sum(pops) / len(pops)
    tag = []
    if fr:
        f = sum(fr) / len(fr)
        tag.append("재료 바닥에 가깝다" if f < 20 else ("총량으로는 남아돈다" if f > 50 else "중간"))
    if sum(occ) / len(occ) > 90:
        tag.append("자리가 찼다")
    print("   %-16s %9.0f %9.1f %11s %11s   %s"
          % (name, p, sum(occ) / len(occ), "-" if not fr else "%.1f" % (sum(fr) / len(fr)),
             "-" if not br else "%.2f" % (sum(br) / len(br)), " · ".join(tag)))
    out["rungs"][name] = {"pop": p, "occupancy": sum(occ) / len(occ),
                          "free_pct": None if not fr else sum(fr) / len(fr),
                          "birth_per_capita_k": None if not br else sum(br) / len(br)}

# ---------- A1 ----------
print("\n-- A1  racld 대 씨앗 rascld  (D = racld − rascld · 95%)")
print("   %-16s %26s %26s %s" % ("칸", "D(1만틱)", "D(2만틱)", "판정(1만틱)"))
A1 = {}
for folder, dens, name in RUNGS:
    d10 = [lnp(D[(folder, "racld", s)], 10000) - lnp(D[(folder, "rascld", s)], 10000) for s in SEEDS]
    d20 = [lnp(D[(folder, "racld", s)], 20000) - lnp(D[(folder, "rascld", s)], 20000) for s in SEEDS]
    m1, l1, h1 = boot(d10)
    m2, l2, h2 = boot(d20)
    A1[name] = {"d": dens, "d10": [m1, l1, h1], "d20": [m2, l2, h2]}
    print("   %-16s %+9.3f [%+.3f, %+.3f] %+9.3f [%+.3f, %+.3f]  %s"
          % (name, m1, l1, h1, m2, l2, h2, verdict(l1, h1, "racld 가 더 잘 자란다", "씨앗이 더 잘 자란다")))
seq = [A1[n]["d10"][0] for _, _, n in DENS_RUNGS]
mono = all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1))
lo4, hi32 = A1["밀도 4"]["d10"][1], A1["밀도 32"]["d10"][2]
sp = A1["space(재료없음)"]["d10"]
d20_pos = [n for _, _, n in RUNGS if A1[n]["d20"][1] > 0]
print("\n   밀도 4→32 의 D(1만) 점추정: %s" % " → ".join("%+.3f" % v for v in seq))
print("   ▸ 단조 감소인가: %s → **%s**" % ("예" if mono else "아니다",
      "재료 가설과 맞다" if mono else "재료 가설과 어긋난다"))
print("   ▸ D(밀도4) 하한 %+.3f vs D(밀도32) 상한 %+.3f → %s"
      % (lo4, hi32, "밀도가 이점을 줄인다" if lo4 > hi32 else "밀도가 이점을 줄인다고 할 수 없다"))
print("   ▸ space 칸 D(1만) %+.3f [%+.3f, %+.3f] → **%s**"
      % (sp[0], sp[1], sp[2], "재료 규칙이 없으면 이점도 없다" if sp[1] <= 0 <= sp[2]
         else ("재료와 무관한 이점이 있다" if sp[1] > 0 else "재료 없으면 씨앗이 낫다")))
print("   ▸ 담는 수(2만틱) — 재료 가설은 D20 > 0(racld 가 더 많이 담는다)을 예측한다: %s"
      % ("**어느 칸에서도 하한 > 0 이 아니다 → 담는 수에서는 racld 의 이점이 없다**" if not d20_pos
         else "하한 > 0 인 칸: " + ", ".join(d20_pos)))
out["A1"] = {"per_rung": A1, "monotone": mono, "lo_d4": lo4, "hi_d32": hi32, "d20_positive_rungs": d20_pos}

# ---------- A2 ----------
def slope_ci(per, codes):
    xs = [len(c) for c in codes]
    xm = sum(xs) / len(xs)
    sxx = sum((x - xm) ** 2 for x in xs)
    means = [sum(per[c]) / len(SEEDS) for c in codes]
    ym = sum(means) / len(means)
    m = sum((x - xm) * (y - ym) for x, y in zip(xs, means)) / sxx
    bs = []
    for ix in IDX:
        ys = [sum(per[c][i] for i in ix) / len(ix) for c in codes]
        ym2 = sum(ys) / len(ys)
        bs.append(sum((x - xm) * (y - ym2) for x, y in zip(xs, ys)) / sxx)
    return means, (m, pct(bs, 0.025), pct(bs, 0.975))


SHORT = ["racld", "rascld", "rascled"]        # acld 를 뺀 세 길이 5·6·7
print("\n-- A2  길이 사다리  (%s · ln P(1만) 을 길이에 회귀 · 기울기 두 벌)"
      % " ".join("%s(%d)" % (c, len(c)) for c in LADDER))
print("   %-16s %s %25s %25s" % ("칸", "".join("%9s" % c for c in LADDER), "기울기 4·5·6·7", "기울기 5·6·7(acld 뺌)"))
A2 = {}
for folder, dens, name in RUNGS:
    per = {c: [lnp(D[(folder, c, s)], 10000) for s in SEEDS] for c in LADDER}
    means, s4 = slope_ci(per, LADDER)
    _, s3 = slope_ci(per, SHORT)
    A2[name] = {"means": means, "slope4": list(s4), "slope3": list(s3)}
    print("   %-16s %s %+8.3f [%+.3f,%+.3f] %+8.3f [%+.3f,%+.3f]"
          % (name, "".join("%9.2f" % v for v in means), s4[0], s4[1], s4[2], s3[0], s3[1], s3[2]))
tight_neg = all(A2[n]["slope4"][2] < 0 or A2[n]["slope3"][2] < 0 for n in ("밀도 4", "밀도 8"))
loose_zero = all(A2[n]["s"][1] <= 0 <= A2[n]["s"][2] for n in ("밀도 32", "space(재료없음)")
                 for A2[n]["s"] in [A2[n]["slope3"]])
print("   ▸ 빡빡한 칸(밀도 4·8)에서 기울기가 음수인가: %s" % ("예" if tight_neg else "아니다"))
print("   ▸ 넉넉한 칸(밀도 32 · space)에서 0 을 품는가(5·6·7): %s — %s"
      % ("예" if loose_zero else "아니다",
         " · ".join("%+.3f [%+.3f, %+.3f]" % tuple(A2[n]["slope3"]) for n in ("밀도 32", "space(재료없음)"))))
print("   ▸ **%s**" % ("짧을수록 유리하고 재료를 주면 사라진다 — 재료 가설과 맞다" if (tight_neg and loose_zero)
                       else ("길이 효과는 있지만 재료를 넉넉히 줘도 그대로다 — 재료 때문이 아니다" if tight_neg
                             else "길이가 짧아서 유리한 것이 아니다")))
out["A2"] = {"per_rung": A2, "tight_negative": tight_neg, "loose_zero": loose_zero}

# ---------- A3 ----------
print("\n-- A3  음성 대조: rascld 대 rsacld (둘 다 6글자)")
print("   %-16s %26s | A1 의 같은 칸" % ("칸", "D(1만)"))
A3, d3 = {}, {}
for folder, dens, name in RUNGS:
    d = [lnp(D[(folder, "rascld", s)], 10000) - lnp(D[(folder, "rsacld", s)], 10000) for s in SEEDS]
    d3[name] = d
    m, lo, hi = boot(d)
    A3[name] = [m, lo, hi]
    note = "  (정밀도 낮아 대조에서 뺌)" if name == "밀도 4" else ""
    print("   %-16s %+9.3f [%+.3f, %+.3f] | %+.3f%s" % (name, m, lo, hi, A1[name]["d10"][0], note))

# 차이의 차이: '밀도 8 → space' 로 D 가 움직인 폭
d1_8 = [lnp(D[("mat_d08", "racld", s)], 10000) - lnp(D[("mat_d08", "rascld", s)], 10000) for s in SEEDS]
d1_sp = [lnp(D[("space", "racld", s)], 10000) - lnp(D[("space", "rascld", s)], 10000) for s in SEEDS]
shift1 = [a - b for a, b in zip(d1_8, d1_sp)]
shift3 = [a - b for a, b in zip(d3["밀도 8"], d3["space(재료없음)"])]
m1s, l1s, h1s = boot(shift1)
m3s, l3s, h3s = boot(shift3)
warn = l3s <= m1s <= h3s
print("   ▸ 밀도 8 → space 로 D 가 움직인 폭 — A1(길이 다른 짝) %+.3f [%+.3f, %+.3f] · A3(길이 같은 짝) %+.3f [%+.3f, %+.3f]"
      % (m1s, l1s, h1s, m3s, l3s, h3s))
print("   ▸ 길이 같은 짝의 구간이 A1 의 값을 품나: %s → %s"
      % ("예" if warn else "아니다",
         "🚩 밀도 효과를 길이 탓으로 읽으면 안 된다" if warn else "길이 같은 짝은 덜 움직인다 — A1·A2 읽기 유지"))
out["A3"] = {"per_rung": A3, "shift_A1": [m1s, l1s, h1s], "shift_A3": [m3s, l3s, h3s], "warn": warn}

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_dens.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_dens.json"))
