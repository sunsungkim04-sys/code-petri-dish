# -*- coding: utf-8 -*-
"""해부 7 회귀 검사 — sim v0.3.3 기본값 · lineage v1.2.0 · dms --find-first · invbud 가 저장된 기록과 맞는가.
   stage7/ 에서 실행(자료는 ../ 에서 읽는다)."""
import json
bad = 0
def chk(name, ok, extra=""):
    global bad
    bad += not ok
    print("  %-52s %s %s" % (name, "같음" if ok else "🚨 다름", extra))

a = json.load(open("out/mono/mono_mat_racld_mu0_s00501.json")); b = json.load(open("../dens/mat_d08/mono_mat_racld_mu0_s00501.json"))
chk("sim 0.3.3 기본값 · mono.js (해부 5A 밀도 8)", a["checksum"] == b["checksum"] and a["samples"] == b["samples"], a["checksum"])
a = json.load(open("out/bud/bud_mat_d08_racld_s00601.json")); b = json.load(open("../bud6/bud_mat_d08_racld_s00601.json"))
chk("sim 0.3.3 기본값 · budget.js (해부 6B)", a["checksum"] == b["checksum"] and a["bins"] == b["bins"], a["checksum"])
a = json.load(open("out/lin/lin_racld_mu0p01_s00501.json")); b = json.load(open("../reg6/lin_racld_mu0p01_s00501.json"))
chk("lineage v1.2.0 기본값 (v1.0.0 기록)", a["checksum"] == b["checksum"] and a["classes"] == b["classes"] and a["series"] == b["series"])
w = json.load(open("out/linw/lin_rascld_mu0p01_s00501.json")); b = json.load(open("../mufine/mono_mat_rascld_mu0p01_s00501.json"))
chk("lineage --watch 는 궤적을 안 바꾼다 (해부 5B 기록)", w["checksum"] == b["checksum"], w["checksum"])
top = dict(w["top_end"]); last = w["series"][-1]
chk("watch 개체 수 = 끝 상위 코드의 racld 수", last[6] == top.get("racld", 0), "%s vs %s" % (last[6], top.get("racld")))
print("     watch 출생 부모:", w["watch_births_by_parent"], "· 처음 새로 생긴 틱", w["watch_first_denovo_tick"], "· 그 밖 부모 상위", w["watch_other_parents_top"][:3])
old = json.load(open("../dms_mu1/dms_mat_racld_s00003_all.json"))
oa = {(x["kind"], x["key"], x["post_seed"]): x["series"] for x in old["arms"]}
d = json.load(open("out/dms/dms_mat_racld_s00003_list.json"))
chk("dms.js 기본값 · sim 0.3.3 (해부 2A 팔 12)", all(oa[(x["kind"], x["key"], x["post_seed"])] == x["series"] for x in d["arms"]))
ff = json.load(open("out/dmsff/dms_mat_racld_s00003_list.json"))
chk("dms --find-first: 배경 체크섬은 그대로", ff["checksum_match"] and ff.get("find_first") is True)
diff = sum(1 for x, y in zip(d["arms"], ff["arms"]) if x["series"] != y["series"])
chk("dms --find-first: 팔 궤적은 달라져야 한다(조작이 먹었나)", diff == len(d["arms"]), "%d/%d 팔 다름" % (diff, len(d["arms"])))
m0 = json.load(open("../dms_main/dms_mat_racld_s00003_all.json"))
m0a = {(x["kind"], x["key"], x["post_seed"]): x["series"] for x in m0["arms"]}
ib = json.load(open("out/ib/ib_inv_s00003_mu0.json"))
same = sum(1 for x in ib["arms"] if m0a.get((x["kind"], x["key"], x["post_seed"])) == x["series"])
chk("invbud 침입 팔의 계통 기록 = 해부 1 (μ0)", ib["checksum_match"] and same == len(ib["arms"]), "%d/%d 팔" % (same, len(ib["arms"])))
for x in ib["arms"]:
    tot = {k: [sum(bb[k][i] for bb in x["bins"]) for i in (0, 1)] for k in ("stepped", "births", "deaths", "copy_ok", "copy_phase", "copy_del", "unknown")}
    per = [tot["copy_phase"][i] / max(tot["copy_ok"][i], 1) for i in (0, 1)]
    print("     %-7s 표지1 출생/틱 %.5f 사망/틱 %.5f 글자당 %.1f틱 · 표지0 출생/틱 %.5f 사망/틱 %.5f 글자당 %.1f틱 · 빠뜨림 %s 모름 %s"
          % (x["key"], tot["births"][1] / tot["stepped"][1], tot["deaths"][1] / tot["stepped"][1], per[1],
             tot["births"][0] / tot["stepped"][0], tot["deaths"][0] / tot["stepped"][0], per[0], tot["copy_del"], tot["unknown"]))
chk("invbud μ0 빠뜨림 0", all(sum(bb["copy_del"][i] for bb in x["bins"]) == 0 for x in ib["arms"] for i in (0, 1)))
lf = json.load(open("out/linff/lin_racld_mu0p01_s00501_ff.json"))
u = lf["wt_parent_mut_births"] / lf["wt_parent_births"]; un = 1 - (1 - 0.01 * (2 / 3 + 10 / 11)) ** 5
print("     재료 먼저 racld μ1%%: u_obs %.4f · u_nom %.4f · 부풀림 %.2f (원래 규칙 4.4)" % (u, un, u / un))
print("회귀 검사:", "전부 통과" if bad == 0 else "🚨 실패 %d" % bad)
