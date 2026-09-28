# -*- coding: utf-8 -*-
"""mono.js v1.1.0 회귀 검사 — 해부 4A 의 한 칸을 다시 돌려 궤적이 같은지 본다."""
import glob
import json
import os

bad = 0
files = sorted(glob.glob("mono_regress5/*.json"))
for f in files:
    b = os.path.basename(f)
    old = json.load(open("mono_curve/" + b))
    new = json.load(open(f))
    same_chk = old["checksum"] == new["checksum"]
    same_rows = len(old["samples"]) == len(new["samples"]) and all(o == n[:6] for o, n in zip(old["samples"], new["samples"]))
    same_top = old["top"] == new["top"]
    ok = same_chk and same_rows and same_top
    bad += (not ok)
    free, tot = new["free_end"], new["mat_total0"]
    print("  %-34s chk %s · 앞6열 %s · top %s · 자유글자 %d/%d (%.1f%%)"
          % (b[:34], same_chk, same_rows, same_top, free, tot, 100.0 * free / tot))
print("파일 %d · 회귀 검사: %s" % (len(files), "통과 — v1.1.0 은 궤적을 바꾸지 않는다" if bad == 0 and files else "🚨 실패 %d" % bad))
