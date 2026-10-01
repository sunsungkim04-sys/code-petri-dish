#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 20 후보 선택 — 기존 여섯 밑 코드와 가족이 다른 밑 코드를 **진화 끝점 기록에서** 규칙으로 고른다.
사전등록: 해부20-사전등록-2026-10-01.md §2 (규칙은 새 코드의 Δ 를 보기 전에 고정 — 이 스크립트가 읽는 파일에는 Δ 가 없다)

  원천(원래 규칙 · 재료 세계 · 밀도 8 의 끝점):  main/mat_d8_*.json(final_counts · μ 1% · 100)  ·  evo17/mat_*.json(final_counts · 80)
                                                lin6/d08/*.json(top_end · 해부 6A)
  세계마다 상위 10 코드만 센다.
  K1  지지: 상위 10 에 든 세계 수 ≥ MIN_WORLDS(20)
  K2  `j`(빌려쓰기) 없음 — 남의 코드를 실행하면 고리 길이가 제 코드의 것이 아니다
  K3  고리 계산 가능(loop_ticks) · `n` 삽입 사다리에 틱 차 1 인 이웃 등급 짝 ≥ 1 · 같은-고리 쌍 ≥ 1(음성 대조)
  K4  기존 여섯과 편집거리(Levenshtein) ≥ 2 — 전부
  K5  고리 모양 (k = 고리 안 `c` 수 · g = 마지막 `c` 와 `l` 사이의 틱 쓰는 글자 수 · 표지 `s` 로 고리 시작 여부) 가
      기존 여섯의 모양 집합에 없다 (여섯은 전부 k 1 · g 0)
  K6  지지 순(세계 수 → 총 개체 수)으로 훑으며, 이미 고른 후보와 **모양이 다르고** 편집거리 ≥ 2 인 것만 받는다 · 최대 MAX_PICK(4)
실행: python3 cand20_select.py [petri 폴더] → _RESULT_cand20.json · 표준 출력
"""
import collections
import glob
import json
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
OLD6 = "racld rascld rsacld racldx acld rascled".split()
MIN_WORLDS = int(os.environ.get("PETRI_MIN_WORLDS", "20")); MAX_PICK, TOPK = 4, 10


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def loop_shape(code):
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    seg = [ch for ch in code[start:li] if ch != "x"]          # l 앞까지 · x 제외
    k = seg.count("c"); last = max(i for i, ch in enumerate(seg) if ch == "c")
    return (k, len(seg) - 1 - last, si >= 0)


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def lev(a, b):
    d = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        p = d[:]; d[0] = i
        for j, cb in enumerate(b, 1):
            d[j] = min(p[j] + 1, d[j - 1] + 1, p[j - 1] + (ca != cb))
    return d[-1]


def ladder_ok(code):
    ins = inserts(code, "n"); Ls = sorted(set(ins.values()))
    adj = any(b - a == 1 for a, b in zip(Ls, Ls[1:]))
    same = any(list(ins.values()).count(L) >= 2 for L in Ls)
    return adj and same, ins


if __name__ == "__main__":
    srcs = [(f, "final_counts") for f in sorted(glob.glob(os.path.join(BASE, "main", "mat_d8_*.json")))]
    srcs += [(f, "final_counts") for f in sorted(glob.glob(os.path.join(BASE, "evo17", "mat_*.json")))]
    srcs += [(f, "top_end") for f in sorted(glob.glob(os.path.join(BASE, "lin6", "d08", "*.json")))]
    W = collections.defaultdict(lambda: [0, 0])
    for f, key in srcs:
        z = json.load(open(f))
        for k, n in (z.get(key) or [])[:TOPK]:
            W[k][0] += 1; W[k][1] += n
    old_shapes = {loop_shape(b) for b in OLD6}
    print("== 해부 20 후보 선택 · 원천 세계 %d · 기존 모양 %s" % (len(srcs), sorted(old_shapes)))
    rows = sorted(W.items(), key=lambda kv: (-kv[1][0], -kv[1][1], kv[0]))
    picked, log = [], []
    for k, (nw, tot) in rows:
        if nw < MIN_WORLDS:
            break
        why = []
        if "j" in k: why.append("K2 j")
        if loop_ticks(k) is None: why.append("K3 고리 없음")
        else:
            ok, _ = ladder_ok(k)
            if not ok: why.append("K3 사다리")
        md = min(lev(k, b) for b in OLD6)
        if md < 2: why.append("K4 거리 %d" % md)
        sh = loop_shape(k) if loop_ticks(k) is not None else None
        if sh is not None and sh in old_shapes: why.append("K5 모양 %s" % (sh,))
        if not why and len(picked) < MAX_PICK:
            if any(loop_shape(p) == sh for p in picked): why.append("K6 모양 중복")
            elif any(lev(k, p) < 2 for p in picked): why.append("K6 거리")
        verdict = "✅ 채택" if not why and len(picked) < MAX_PICK else ("(한도)" if not why else "✗ " + " · ".join(why))
        if not why and len(picked) < MAX_PICK:
            picked.append(k)
        log.append({"code": k, "worlds": nw, "total": tot, "loop": loop_ticks(k), "shape": sh, "verdict": verdict})
        if not why or md >= 2:
            print("  %-10s 세계 %4d · 개체 %8d · 고리 %s · 모양 %s · %s" % (k, nw, tot, loop_ticks(k), sh, verdict))
    print("== 채택 %d: %s" % (len(picked), " ".join(picked)))
    for p in picked:
        ins = inserts(p, "n")
        print("   %s: n 코드 %d · 고리 등급 %s" % (p, len(ins), dict(sorted(collections.Counter(ins.values()).items()))))
    json.dump({"picked": picked, "sources": len(srcs), "log": log, "rule": {"MIN_WORLDS": MIN_WORLDS, "MAX_PICK": MAX_PICK, "TOPK": TOPK}},
              open(os.path.join(BASE, "_RESULT_cand20.json"), "w"), ensure_ascii=False, indent=1)
