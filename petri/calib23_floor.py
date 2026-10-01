#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 23 설계 근거(합성 · 실자료 아님): 잡음 바닥 q95 추정기 둘의 보정 — 팔 11 · 배경 20 · 반복 200.
사전등록 해부23 §2 의 합성 점검 숫자의 재생성기(검토 S4 — 원래 scratchpad calib23.py 를 옮김).
판정기 noise23_analyze.py 의 floor_q · floor_q_emp · floor_q_boot · pct 를 그대로 불러 쓴다(같은 폴더의 파일 · 'if SELFTEST:' 앞까지만 실행 — 파일을 읽지 않는 부분).
실행: python3 calib23_floor.py   (파일 읽기 없음 · 수 분 · 단일 스레드)
"""
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = ["calib23_floor", "--noop"]
src = open(os.path.join(HERE, "noise23_analyze.py"), encoding="utf-8").read().split("if SELFTEST:")[0]
ns = {"__file__": os.path.join(HERE, "noise23_analyze.py")}
exec(src, ns)
pct = ns["pct"]; seeds = list(range(20))


def world(rng, het):
    W = {}
    for s in seeds:
        bg = rng.gauss(0, 1); sd = 0.3 * (math.exp(rng.gauss(0, 0.3)) if het else 1)
        W[s] = [bg + rng.gauss(0, sd) for _ in range(11)]
    return W


for het, truth in ((False, 1.959964 * math.sqrt(2 * 0.09 / 20)), (True, 1.959964 * math.sqrt(2 * 0.09 * math.exp(0.18) / 20))):
    rng = random.Random(3 if het else 1); E, P, H = [], [], []
    for rep in range(200):
        W = world(rng, het); E.append(ns["floor_q_emp"](W, seeds)); P.append(ns["floor_q"](W, seeds))
        if rep < 40: H.append(ns["floor_q_boot"](W, seeds, "c%d" % rep))
    print("het=%s 참 q95 %.4f | 경험(55쌍) 중앙 %.4f [%.4f, %.4f] | 합친분산 중앙 %.4f [%.4f, %.4f] | 합친분산 재추출 상한 ≥ 참 %d/40 · 하한 ≤ 참 %d/40 · 상한/참 중앙 %.3f"
          % (het, truth, pct(E, .5), pct(E, .025), pct(E, .975), pct(P, .5), pct(P, .025), pct(P, .975), sum(h >= truth for l, h in H), sum(l <= truth for l, h in H), pct([h / truth for l, h in H], .5)))
