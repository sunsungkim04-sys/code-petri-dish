# -*- coding: utf-8 -*-
"""해부 11 회귀 검사 — budget.js v1.2.0 · dms.js(--age0) 기본값이 저장된 기록과 맞는가 (서버 stage11/ 에서 실행)

  budget.js v1.2.0 기본값   : 해부 7 ① bud6c 밀도 8 racld 시드 661 기록과 체크섬 · 칸마다 옛 필드 전부 같다
  dms.js 기본값 (--exact 1) : 해부 6C inv6 배경 3 μ 0.5% 팔 계통 기록이 같다
  dms.js --ancestor racld   : 해부 10 G inv10anc 시드 1101 μ 0 기록과 같다
  새 옵션 (서술)            : budget --age0 1200 은 체크섬이 달라야 한다 · 죽음 수 · 죽은 나이 · 자유 재료가 적힌다
실행: cd ~/petri/stage11 && python3 reg11_check.py
"""
import json

bad = 0
OLD = ['stepped', 'not_stepped', 'unknown', 'copy_ok', 'copy_stall', 'copy_del', 'copy_idle', 'alloc_ok', 'alloc_blocked',
       'alloc_idle', 'alloc_dead', 'div_with_child', 'div_idle', 'births', 't0', 'pop0', 'by_ins']


def chk(name, ok):
    global bad
    print(("  ✅ " if ok else "  🚨 ") + name)
    if not ok:
        bad += 1


print("== budget.js v1.2.0 기본값 = bud6c 기록")
a = json.load(open("reg_bud/bud_mat_d08_racld_s00661.json")); b = json.load(open("bud6c/bud_mat_d08_racld_s00661.json"))
same = len(a["bins"]) == len(b["bins"]) and all(x[k] == y[k] for x, y in zip(a["bins"], b["bins"]) for k in OLD)
chk("체크섬 %s = %s · 칸 %d 옛 필드 전부 같음 %s · b%s" % (a["checksum"], b["checksum"], len(a["bins"]), same, a["budget_version"]),
    same and a["checksum"] == b["checksum"] and a["budget_version"] == "1.2.0")
d = sum(x["deaths"] for x in a["bins"]); bi = sum(x["births"] for x in a["bins"])
chk("새 필드: 죽음 %d · 출생 %d · 끝 개체 %d → 출생 − 죽음 = %d (시조 1 포함 개체 수와 1 이내)" % (d, bi, a["pop_end"], bi - d), abs(bi + 1 - d - a["pop_end"]) <= 1)
f0 = a["bins"][0]["free0"]; chk("자유 재료 첫 칸 합 %d = 처음 총량 %d − 시조 글자 5" % (sum(f0), a["mat_total0"]), sum(f0) == a["mat_total0"] - 5)

print("== dms.js 기본값 = inv6 기록 (배경 3 · μ 0.5%)")
a = json.load(open("reg_dms/dms_mat_racld_s00003_list.json")); b = json.load(open("inv6/mu0p005/dms_mat_racld_s00003_list.json"))
same = a["grow_checksum"] == b["grow_checksum"] and len(a["arms"]) == len(b["arms"]) and all(x["series"] == y["series"] and x["key"] == y["key"] for x, y in zip(a["arms"], b["arms"]))
chk("배경 체크섬 %s · 팔 %d 기록 전부 같음 %s · age0 필드 없음 %s" % (a["grow_checksum"], len(a["arms"]), same, "age0" not in a), same and "age0" not in a)

print("== dms.js --ancestor racld = inv10anc 기록 (시드 1101 · μ 0)")
a = json.load(open("reg_anc/dms_mat_racld_s01101_list.json")); b = json.load(open("inv10anc/mu0/dms_mat_racld_s01101_list.json"))
same = a["grow_checksum"] == b["grow_checksum"] and all(x["series"] == y["series"] for x, y in zip(a["arms"], b["arms"]))
chk("배경 체크섬 %s · 팔 %d 기록 전부 같음 %s" % (a["grow_checksum"], len(a["arms"]), same), same)

print("== 새 옵션 (서술)")
c = json.load(open("reg_bud/bud_mat_d08_racld_s00661_a1200.json")); a = json.load(open("reg_bud/bud_mat_d08_racld_s00661.json"))
chk("budget --age0 1200: 체크섬 %s ≠ 기본 %s · opts.age0 %s" % (c["checksum"], a["checksum"], c["opts"].get("age0")), c["checksum"] != a["checksum"] and c["opts"].get("age0") == 1200)
for z, nm in ((a, "수명 300"), (c, "수명 1200")):
    B = [x for x in z["bins"] if x["t0"] >= 10000]
    ok = sum(x["copy_ok"] for x in B); st = sum(x["copy_stall"] for x in B); bl = sum(x["alloc_blocked"] for x in B); stp = sum(x["stepped"] for x in B)
    de = sum(x["deaths"] for x in B); da = sum(x["death_age_sum"] for x in B)
    print("    %s · 개체 %d · 헛손질/성공 %.2f · 헛손질 몫 %.3f · 막힘 몫 %.3f · 죽은 나이 평균 %.0f" % (nm, B[-1]["pop0"], st / ok, st / stp, bl / stp, da / max(1, de)))
print("문제 %d" % bad)
