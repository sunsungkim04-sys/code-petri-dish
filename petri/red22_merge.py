#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 22 — 시드별 reduced_d.py 출력(s<seed>.json)을 원래 판(_RESULT_reduced_d.json)과 같은 틀의 한 파일로 합친다.
합친 파일은 고치지 않은 reduced_d_analyze.py 에 그대로 들어간다(원래 판 재현 확인용).
  args 는 첫 파일의 것에 seeds = 파일 수 · seed0 = 가장 작은 시드 · out = 합친 파일 경로 를 넣는다.
  arm · cells 는 시드 오름차순으로 잇는다(원래 판도 시드 순 → 규칙 → μ → 코드 순).
실행: python3 red22_merge.py <폴더> <합친.json>
"""
import glob
import json
import os
import sys

src, dst = sys.argv[1], sys.argv[2]
files = glob.glob(os.path.join(src, "s*.json"))
runs = sorted((json.load(open(f)) for f in files), key=lambda z: z["args"]["seed0"])
if not runs: sys.exit("합칠 파일이 없다: %s" % src)
keys = [k for k in runs[0]["args"] if k not in ("seeds", "seed0", "out")]
for z in runs:
    if z["args"]["seeds"] != 1: sys.exit("시드 %d 파일이 시드 1 개짜리가 아니다" % z["args"]["seed0"])
    for k in keys:
        if z["args"][k] != runs[0]["args"][k]: sys.exit("시드 %d 의 설정 %s 가 다르다" % (z["args"]["seed0"], k))
    if z["codes"] != runs[0]["codes"]: sys.exit("시드 %d 의 코드 목록이 다르다" % z["args"]["seed0"])
args = dict(runs[0]["args"]); args["seeds"] = len(runs); args["seed0"] = runs[0]["args"]["seed0"]; args["out"] = dst
OUT = {"args": args, "cells": [c for z in runs for c in z["cells"]], "codes": runs[0]["codes"], "arm": [a for z in runs for a in z["arm"]]}
json.dump(OUT, open(dst, "w"), indent=1)
print("합침 %d 파일 · 시드 %s · arm %d → %s" % (len(runs), [z["args"]["seed0"] for z in runs], len(OUT["arm"]), dst))
