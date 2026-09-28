# -*- coding: utf-8 -*-
"""해부 6 계측기 회귀 검사 — 새 계측기가 궤적을 바꾸지 않는가, dms.js 기본값이 전과 같은가."""
import json
pairs = [
    ("reg6/lin_racld_mu0p01_s00501.json", "mufine/mono_mat_racld_mu0p01_s00501.json"),
    ("reg6/lin_rascld_mu0p005_s00502.json", "mufine/mono_mat_rascld_mu0p005_s00502.json"),
    ("reg6/bud_mat_d04_racld_s00501.json", "dens/mat_d04/mono_mat_racld_mu0_s00501.json"),
    ("reg6/bud_mat_d08_rascld_s00503.json", "dens/mat_d08/mono_mat_rascld_mu0_s00503.json"),
    ("reg6/bud_space_rascld_s00502.json", "dens/space/mono_space_rascld_mu0_s00502.json"),
]
bad = 0
for a, b in pairs:
    x, y = json.load(open(a)), json.load(open(b))
    ok = x["checksum"] == y["checksum"]
    pe = x.get("pop_end", x["series"][-1][1] if "series" in x else None)
    ok2 = pe == y["samples"][-1][1]
    bad += not (ok and ok2)
    print("  %-40s 체크섬 %s(%s) · 끝 개체 %s(%s)" % (a[5:], ok, x["checksum"], ok2, pe))

old = json.load(open("dms_mu1/dms_mat_racld_s00003_all.json"))
pl = json.load(open("reg6/plain/dms_mat_racld_s00003_list.json"))
ex = json.load(open("reg6/exact/dms_mat_racld_s00003_list.json"))
oa = {(a["kind"], a["key"], a["post_seed"]): a for a in old["arms"]}
same = 0
for a in pl["arms"]:
    o = oa[(a["kind"], a["key"], a["post_seed"])]
    same += o["series"] == a["series"]
print("  dms 기본값: 팔 %d 중 옛 기록(dms_mu1)과 같은 궤적 %d · 체크섬 %s" % (len(pl["arms"]), same, pl["grow_checksum"] == old["grow_checksum"]))
bad += same != len(pl["arms"])
strip = lambda z: {k: v for k, v in z.items() if k not in ("elapsed_ms", "arms", "grow_ms")}
e2 = sum(1 for a, b in zip(pl["arms"], ex["arms"]) if [r[:3] for r in b["series"]] == a["series"])
print("  dms --exact: 앞 3칸이 기본값과 같은 팔 %d/%d · 출력 키 차이 %s" % (e2, len(pl["arms"]), sorted(set(ex) - set(pl))))
bad += e2 != len(pl["arms"])
r = [a for a in ex["arms"] if a["kind"] == "mut"][0]["series"]
print("  --exact 예시(씨앗 팔 μ1%): 끝 줄", r[-1], "→ 계통 안 씨앗 그대로 %.1f%%" % (100.0 * r[-1][3] / max(r[-1][2], 1)))
print("회귀 검사:", "전부 통과" if bad == 0 else "🚨 실패 %d" % bad)
