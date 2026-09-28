# -*- coding: utf-8 -*-
"""해부 9 회귀 검사 — sim v0.3.4 기본값 · lineage v1.4.0 · dms --remember-die 가 저장된 기록과 맞는가(서버 stage9/ 에서 실행).

  replay.js (v0.3.4 기본값)        : 되감기 2 mat 시드 101~103 기록과 체크섬 · 기록 전부 같다
  lineage.js v1.4.0 (--remember-die 끔): 해부 7 계보 lin7/orig 시드 701(씨앗 μ 1%) 기록과 체크섬 · 계보 · 기록이 같다
  dms.js (기본값 · --exact 1)      : 해부 6C inv6 배경 3 μ 0.5% 의 팔 14 계통 기록이 같다
  dms.js --remember-die 1          : 배경 체크섬은 그대로(배경은 원래 규칙) · 팔 궤적은 바뀐다(μ > 0 에서 달라야 하고, μ 0 에서도 난수 소비가 달라 갈라진다 — 서술)
  lineage.js --remember-die --spec : 돌아가고, R(주사위 몫)이 1 근처로 내려간다(서술 — 판정선은 노트)
실행: cd ~/petri/stage9 && python3 reg9_check.py
"""
import json
import glob

bad = 0


def chk(name, ok):
    global bad
    print(("  ✅ " if ok else "  🚨 ") + name)
    if not ok:
        bad += 1


print("== replay.js v0.3.4 기본값 = 되감기 2 기록")
for s in (101, 102, 103):
    a = json.load(open("reg_replay/mat_d8_mu0p01_s%05d.json" % s))
    b = json.load(open("rewind2/mat_d8_mu0p01_s%05d.json" % s))
    same = all(a[k] == b[k] for k in ("checksum", "samples", "tops", "specials", "snap", "final_counts", "births", "deaths"))
    chk("시드 %d · 체크섬 %s = %s · 기록 전부 같음 %s · sim %s" % (s, a["checksum"], b["checksum"], same, a["sim_version"]), same and a["sim_version"] == "0.3.4")

print("== lineage.js v1.4.0 기본값 = 해부 7 기록")
a = json.load(open("reg_lin/lin_rascld_mu0p01_s00701.json"))
b = json.load(open("lin7/orig/lin_rascld_mu0p01_s00701.json"))
keys = ("checksum", "classes", "series", "watch_births_by_parent", "watch_class", "same_kids_wt", "same_kids_watch", "wt_parent_births", "wt_parent_mut_births", "top_end", "missed_births")
diff = [k for k in keys if a[k] != b[k]]
chk("시드 701 씨앗 μ 1%% · 다른 필드 %s · remember_die %s" % (diff or "없음", a["remember_die"]), not diff and a["remember_die"] is False)

print("== dms.js 기본값 = 해부 6C 기록 (배경 3 · μ 0.5%)")
a = json.load(open("reg_dms/mu0p005/dms_mat_racld_s00003_list.json"))
b = json.load(open("inv6/mu0p005/dms_mat_racld_s00003_list.json"))
same = sum(1 for x, y in zip(a["arms"], b["arms"]) if x["series"] == y["series"] and x["key"] == y["key"])
chk("배경 체크섬 %s · 팔 계통 기록 같음 %d/%d · remember_die 없음 %s" % (a["checksum_match"], same, len(b["arms"]), "remember_die" not in a),
    a["checksum_match"] and same == len(b["arms"]) == 14 and "remember_die" not in a)

print("== dms.js --remember-die 1 (배경 3 · μ 0.5%)")
c = json.load(open("reg_dms_rd/mu0p005/dms_mat_racld_s00003_list.json"))
diffn = sum(1 for x, y in zip(c["arms"], b["arms"]) if x["series"] != y["series"])
chk("배경 체크섬 그대로 %s · remember_die 표시 %s · 팔 궤적 다름 %d/14 (μ > 0 이라 달라야 함)" % (c["checksum_match"], c.get("remember_die"), diffn),
    c["checksum_match"] and c.get("remember_die") is True and diffn >= 13)
c0 = json.load(open("reg_dms_rd/mu0/dms_mat_racld_s00003_list.json"))
b0 = json.load(open("inv6/mu0/dms_mat_racld_s00003_list.json"))
same0 = sum(1 for x, y in zip(c0["arms"], b0["arms"]) if x["series"] == y["series"])
# 09-17 첫 실행에서 'μ 0 이면 궤적이 같아야 한다' 로 걸었다가 0/14 로 실패했다 — 내 기대가 틀렸다: 기억하는 주사위 세계는 헛손질 동안
# 난수를 안 굴리므로 결과(글자)는 같아도 난수 소비가 달라 궤적이 갈라진다. 판정선이 아니라 서술로 둔다(집합 수준 점검은 inv9_analyze M3).
print("   (서술) μ 0 에서 팔 궤적이 원래와 같은 팔 %d/14 — 다르게 나오는 것이 맞다(헛손질 동안 난수를 안 굴려 난수 소비가 다르다)" % same0)

print("== lineage.js --remember-die --spec (서술)")
for f in sorted(glob.glob("reg_lin_rd/*.json")):
    z = json.load(open(f))
    sp = z["spec"]
    adv = sp["ok"] + sum(sp["sub"]) + sp["del"]
    R = (sp["del"] / adv) / (z["mu"] / 3) if adv else float("nan")
    F = (1 + sum(sp["ins"]) / sp["del"] + sum(sp["sub"]) / sp["del"]) / (2 + 30 / 11) if sp["del"] else float("nan")
    u = z["wt_parent_mut_births"] / max(1, z["wt_parent_births"])
    L = len(z["code"])
    unom = 1 - (1 - z["mu"] * (2 / 3 + 10 / 11)) ** L
    print("   %s · remember_die %s · R %.2f · F %.2f · I_계보 %.2f · 헛손질/맞게 %.1f · 멸종 %d" % (f.split("/")[-1], z["remember_die"], R, F, u / unom, sp["stall"] / max(1, sp["ok"]), z["extinct_at"]))

print("\n회귀 검사: %s" % ("통과" if bad == 0 else "🚨 실패 %d" % bad))
