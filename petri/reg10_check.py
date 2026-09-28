# -*- coding: utf-8 -*-
"""해부 10 회귀 검사 — sim v0.3.5 기본값(cosmic · diffuse · age0/ageVar) · lineage v1.5.0 · dms(--cosmic 0 · --ancestor) 가 저장된 기록과 맞는가(서버 stage10/ 에서 실행).

  replay.js (v0.3.5 기본값)   : 되감기 2 mat 시드 101~103 기록과 체크섬 · 기록 전부 같다
  lineage.js v1.5.0 기본값     : 해부 7 계보 lin7/orig 시드 701 기록과 같다 (opts 필드는 새로 더함)
  dms.js 기본값 (--exact 1)    : 해부 6C inv6 배경 3 μ 0.5% 팔 14 계통 기록이 같다
  새 옵션 (서술)               : lineage --cosmic 0 / --diffuse 0.05 · 0.5 / --age0 600 · 1200 의 R · F · 부풀림 · 순도 · 멸종 —
                                기본값 파일과 체크섬이 달라야 한다(규칙이 다르므로) · cosmic 0 이면 cosmic_wt_changed = 0
실행: cd ~/petri/stage10 && python3 reg10_check.py
"""
import glob
import json

bad = 0


def chk(name, ok):
    global bad
    print(("  ✅ " if ok else "  🚨 ") + name)
    if not ok:
        bad += 1


print("== replay.js v0.3.5 기본값 = 되감기 2 기록")
for s in (101, 102, 103):
    a = json.load(open("reg_replay/mat_d8_mu0p01_s%05d.json" % s))
    b = json.load(open("rewind2/mat_d8_mu0p01_s%05d.json" % s))
    same = all(a[k] == b[k] for k in ("checksum", "samples", "tops", "specials", "snap", "final_counts", "births", "deaths"))
    chk("시드 %d · 체크섬 %s = %s · 기록 전부 같음 %s · sim %s" % (s, a["checksum"], b["checksum"], same, a["sim_version"]), same and a["sim_version"] == "0.3.5")

print("== lineage.js v1.5.0 기본값 = 해부 7 기록")
a = json.load(open("reg_lin/lin_rascld_mu0p01_s00701.json"))
b = json.load(open("lin7/orig/lin_rascld_mu0p01_s00701.json"))
keys = ("checksum", "classes", "series", "watch_births_by_parent", "watch_class", "same_kids_wt", "same_kids_watch", "wt_parent_births", "wt_parent_mut_births", "top_end", "missed_births")
diff = [k for k in keys if a[k] != b[k]]
chk("시드 701 씨앗 μ 1%% · 다른 필드 %s · opts %s" % (diff or "없음", {k: a["opts"].get(k) for k in ("cosmic", "diffuse", "age0", "ageVar")}), not diff)

print("== dms.js 기본값 = 해부 6C 기록 (배경 3 · μ 0.5%)")
a = json.load(open("reg_dms/mu0p005/dms_mat_racld_s00003_list.json"))
b = json.load(open("inv6/mu0p005/dms_mat_racld_s00003_list.json"))
same = sum(1 for x, y in zip(a["arms"], b["arms"]) if x["series"] == y["series"] and x["key"] == y["key"])
chk("배경 체크섬 %s · 팔 계통 기록 같음 %d/%d · 새 표시 없음 %s" % (a["checksum_match"], same, len(b["arms"]), not any(k in a for k in ("cosmic_off", "ancestor"))),
    a["checksum_match"] and same == len(b["arms"]) == 14 and not any(k in a for k in ("cosmic_off", "ancestor")))

print("== 새 옵션 (서술 · 시드 9501 · μ 1% · 밀도 8)")
base = None
for f in sorted(glob.glob("pilot10/lin/*.json")):
    z = json.load(open(f))
    sp = z["spec"]
    adv = sp["ok"] + sum(sp["sub"]) + sp["del"]
    R = (sp["del"] / adv) / (z["mu"] / 3) if adv else float("nan")
    F = (1 + sum(sp["ins"]) / sp["del"] + sum(sp["sub"]) / sp["del"]) / (2 + 30 / 11) if sp["del"] else float("nan")
    u = z["wt_parent_mut_births"] / max(1, z["wt_parent_births"])
    L = len(z["code"])
    unom = 1 - (1 - z["mu"] * (2 / 3 + 10 / 11)) ** L
    pur = (z["top_end"][0][1] if z["top_end"] and z["top_end"][0][0] == z["code"] else 0) / max(1, sum(c for _, c in z["top_end"]))
    o = z["opts"]
    tag = "cosmic %s · diffuse %s · age %s+%s" % (o.get("cosmic"), o.get("diffuse"), o.get("age0"), o.get("ageVar"))
    print("   %-38s %s · R %.2f · F %.2f · 부풀림 %.2f · 우주선 변화 %d · 끝 WT 몫(상위8 중) %.2f · 멸종 %d · chk %s"
          % (f.split("/")[-1], tag, R, F, u / unom, z["cosmic_wt_changed"], pur, z["extinct_at"], z["checksum"]))
    if f.endswith("s09501.json"):
        base = z["checksum"]
    if "_c0" in f:
        chk("cosmic 0 → 우주선 변화 0", z["cosmic_wt_changed"] == 0)
    if base and not f.endswith("s09501.json"):
        chk("%s 체크섬이 기본값과 다름" % f.split("/")[-1], z["checksum"] != base)

print("\n회귀 검사: %s" % ("통과" if bad == 0 else "🚨 실패 %d" % bad))
