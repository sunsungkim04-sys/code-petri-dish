#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""강화계획 ④ — 교차 μ* 를 A(해부 14)의 U 와 해부 12 의 변종 구성으로 예측해 B(해부 15)의 실측과 대조한다. 🟡 사후 서술 · 판정 아님.
  단순 모형(해부 12 M1 의 뼈대): ΔΔ(μ) ≈ ΔΔ(0) + G · λ · (U_a(μ) − U_c(μ)),  U(μ) = U(0.3%) · μ/0.003
    a = 짧은 고리 · c = 긴 고리 · G = 3000/420 = 7.14 · λ = 잃는 변종 몫 (불임 몫 σ, 또는 1)
  → μ* = 0.003 · (−ΔΔ(0)) / (G · λ · (U_a − U_c))
  대조: B 의 ΔΔ(0.3%) − ΔΔ(0) 실측 대 모형 G·λ·(U_a − U_c) · 교차 μ* 실측(선형 보간) 대 예측
"""
import json, math
p14 = json.load(open("_RESULT_pairs14.json")); p15 = json.load(open("_RESULT_pairs15.json")); m12 = json.load(open("_RESULT_model12.json"))
G = 3000 / 420
sigma = {b: m12["flow|%s" % b]["cls"]["불임"] for b in ["racld", "rascld", "rsacld", "racldx"]}
sigma_default = sum(sigma.values()) / len(sigma)
rows = []
for k, v in p15.items():
    if not k.startswith("B1|"): continue
    _, b, a, c = k.split("|")
    Ua = p14["code|%s|%s" % (b, a)]["U_copy"]; Uc = p14["code|%s|%s" % (b, c)]["U_copy"]
    dd0 = v["dd_by_mu"]["0.0"]; dd3 = v["dd_by_mu"]["0.003"]
    lam = sigma.get(b, sigma_default)
    rise_model = G * lam * (Ua - Uc); rise_meas = dd3 - dd0
    mu_pred = 0.003 * (-dd0) / rise_model if rise_model > 0 else float("nan")
    rows.append((b, a, c, dd0, rise_meas, rise_model, v["cross"], mu_pred, lam))
print("== 교차 μ* — 단순 모형 대 실측 (λ = 가족의 불임 몫 · G %.2f) · 🟡 사후 서술" % G)
print("   %-8s %-9s %-9s %8s %10s %10s %9s %9s" % ("밑", "짧은", "긴", "ΔΔ(0)", "상승 실측", "상승 모형", "μ* 실측", "μ* 모형"))
byb = {}
for b, a, c, dd0, rm, rmod, cross, mp, lam in rows:
    print("   %-8s %-9s %-9s %+8.3f %+10.3f %+10.3f %9s %9s" % (b, a, c, dd0, rm, rmod, ("%.4f" % cross) if cross else "밖", ("%.4f" % mp) if mp == mp else "—"))
    byb.setdefault(b, []).append((rm, rmod, cross, mp))
print("\n== 밑 코드별 요약: 상승 실측/모형 비(중앙) · μ* 실측(중앙 · 구간 안 짝 수) · μ* 모형(중앙)")
for b, L in byb.items():
    ratio = sorted(rm / rmod for rm, rmod, _, _ in L if rmod > 0); cr = sorted(c for _, _, c, _ in L if c); mp = sorted(m for _, _, _, m in L if m == m)
    print("   %-8s 비 %.2f (범위 %.2f~%.2f) · μ* 실측 %s (%d/%d 짝) · μ* 모형 %s" % (b, ratio[len(ratio) // 2], ratio[0], ratio[-1], ("%.4f" % cr[len(cr) // 2]) if cr else "—", len(cr), len(L), ("%.4f" % mp[len(mp) // 2]) if mp else "—"))
allr = sorted(rm / rmod for L in byb.values() for rm, rmod, _, _ in L if rmod > 0)
print("\n   전체 60 짝: 상승 실측/모형 중앙 %.2f · 사분위 [%.2f, %.2f] — 단순 모형은 충실도 몫을 그만큼 덜 잡는다(해부 12 M1 과 같은 방향)" % (allr[len(allr) // 2], allr[len(allr) // 4], allr[3 * len(allr) // 4]))
json.dump({"G": G, "sigma": sigma, "rows": rows, "ratio_median": allr[len(allr) // 2]}, open("_RESULT_mustar.json", "w"), indent=1, ensure_ascii=False)
