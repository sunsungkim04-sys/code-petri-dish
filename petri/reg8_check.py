# -*- coding: utf-8 -*-
"""해부 8 회귀 검사 — 새 계측기가 저장된 기록과 맞는가(서버 stage8/ 에서 실행).

  replay.js (--find-first 추가) 기본값  : 되감기 2 mat 시드 101~103 기록과 체크섬 · samples · tops · final_counts 가 같다
  lifehist.js --ticks 3000              : 해부 1(dms_main) 배경 3 · 4 의 ref · neutral 10 · rascld 팔 계통 기록(series)이 같다
                                          + 칸 기록의 셈이 맞다: 칸 끝 표지1 = 칸 시작 + 출생 − 사망 (칸마다)
  lineage.js v1.3.0 --spec 끔           : 해부 7 계보(lin7) 기록과 체크섬 · classes · series · watch 가 같다(원래 · 재료 먼저)
  lineage.js v1.3.0 --spec 켬           : 끈 판과 체크섬이 같다(읽기만) · 결과 합이 0 이 아니다
실행: cd ~/petri/stage8 && python3 reg8_check.py  →  _RESULT_reg8_check.txt 로 옮긴다
"""
import json

bad = 0


def chk(name, ok):
    global bad
    print(("  ✅ " if ok else "  🚨 ") + name)
    if not ok:
        bad += 1


print("== replay.js 기본값 = 되감기 2 기록")
for s in (101, 102, 103):
    a = json.load(open("reg_replay/mat_d8_mu0p01_s%05d.json" % s))
    b = json.load(open("rewind2/mat_d8_mu0p01_s%05d.json" % s))
    same = all(a[k] == b[k] for k in ("checksum", "samples", "tops", "specials", "snap", "final_counts", "births", "deaths"))
    chk("시드 %d · 체크섬 %s = %s · 기록 전부 같음 %s · find_first 필드 없음 %s" % (s, a["checksum"], b["checksum"], same, "find_first" not in a),
        same and "find_first" not in a)

print("== lifehist.js 3,000틱 = 해부 1 팔 기록")
for s in (3, 4):
    a = json.load(open("reg_lh/lh_mat_s%05d_t3000.json" % s))
    b = json.load(open("dms_main/dms_mat_racld_s%05d_all.json" % s))
    chk("배경 %d 체크섬 대조 %s" % (s, a["checksum_match"]), a["checksum_match"])
    bb = {}
    for arm in b["arms"]:
        lab = arm["labels"][0] if arm["kind"] != "mut" else arm["key"]
        bb[(arm["kind"], lab)] = arm["series"]
    n_same = 0
    for arm in a["arms"]:
        lab = arm["labels"][0]
        if bb.get((arm["kind"], lab)) == arm["series"]:
            n_same += 1
    chk("배경 %d 팔 계통 기록 같음 %d/%d" % (s, n_same, len(a["arms"])), n_same == len(a["arms"]) == 12)
    ok_id = True
    for arm in a["arms"]:
        ser = {r[0]: r for r in arm["series"]}
        for bn in arm["bins"]:
            t1 = bn["t"] + a["bin"]
            if t1 not in ser:
                continue
            for L in (0, 1):
                start = bn["n0"][L]
                end = ser[t1][2] if L == 1 else ser[t1][1] - ser[t1][2]
                if start + sum(bn["b"][L]) - sum(bn["da"][L]) - sum(bn["do"][L]) != end:
                    ok_id = False
        if arm["lost_parent"]:
            ok_id = False
    chk("배경 %d 칸 셈(끝 = 시작 + 출생 − 사망) 전부 맞음 · 부모 잃음 0" % s, ok_id)

print("== lineage.js v1.3.0")
pairs = [("reg_lin/off/lin_rascld_mu0p01_s00701.json", "lin7/orig/lin_rascld_mu0p01_s00701.json"),
         ("reg_lin/offff/lin_racld_mu0p005_s00702_ff.json", "lin7/ff/lin_racld_mu0p005_s00702_ff.json")]
for new, old in pairs:
    a, b = json.load(open(new)), json.load(open(old))
    keys = ("checksum", "classes", "series", "watch_births_by_parent", "watch_class", "same_kids_wt", "same_kids_watch",
            "wt_parent_births", "wt_parent_mut_births", "top_end", "missed_births")
    diff = [k for k in keys if a[k] != b[k]]
    chk("%s = 해부 7 기록 (다른 필드 %s) · spec 없음 %s" % (new, diff or "없음", a["spec"] is None), not diff and a["spec"] is None)
for on, off in (("reg_lin/on/lin_rascld_mu0p01_s00701.json", "reg_lin/off/lin_rascld_mu0p01_s00701.json"),
                ("reg_lin/onff/lin_racld_mu0p005_s00702_ff.json", "reg_lin/offff/lin_racld_mu0p005_s00702_ff.json")):
    a, b = json.load(open(on)), json.load(open(off))
    sp = a["spec"]
    tot = sp["ok"] + sp["del"] + sp["stall"] + sum(sp["sub"]) + sum(sp["ins"])
    chk("%s 체크섬 = 끈 판 (%s) · 결과 합 %d · 맞게 %d · 빠뜨림 %d · 바뀜 %d · 끼어듦 %d · 헛손질 %d · 그 밖 %d · 죽음 %d"
        % (on, a["checksum"] == b["checksum"], tot, sp["ok"], sp["del"], sum(sp["sub"]), sum(sp["ins"]), sp["stall"], sp["other"], sp["dead"]),
        a["checksum"] == b["checksum"] and a["classes"] == b["classes"] and tot > 0)

print("\n회귀 검사: %s" % ("통과" if bad == 0 else "🚨 실패 %d" % bad))
