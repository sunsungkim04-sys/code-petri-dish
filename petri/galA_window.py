#!/usr/bin/env python3
"""galA_window.py — 갈래 A: mini27 엔진에 '시간 창별 계수기' 만 덧댄 래퍼 (mini27.py 무수정 · import 해서 씀)

왜: 저쪽(Heinrich-Mora & Feldman, arXiv 2511.03073 v2)의 u_{2,t} 는 세대마다 한 유전형 전체가 공유하는 율이다
("Stochasticity acts only through time, not across individuals"). mini27 계수기는 (lineage, L) 별 assay 전체 합계라
창 사이 분포를 못 준다(mini27.py:217 · 714). 이 래퍼는 같은 증분을 (lineage, L, WBIN 틱 칸, 틱 홀짝) 별로도 쌓는다.
홀짝은 '같은 창 안 독립인 두 반쪽' 을 만들려는 것이다 — 두 반쪽 율의 공분산은 표본 잡음 모형 없이 창 간 분산을 준다.

바꿔 끼우는 것(전부 실행 시점에 mini27 모듈 속성을 바꾼다 · 파일은 안 건드린다)
  1. World.counter → Tee: 동결 판과 같은 합계 목록(self.cnt)을 그대로 올리면서 창 목록도 같이 올린다.
     난수를 소비하지 않으므로 동역학과 본 출력 JSON 은 mini27 그대로다 — 가드 G-A1 이 바이트로 확인한다.
  2. build · build_grid · grow · grow_grid → 만든 세계를 붙잡는다.
  3. run_arm · run_arm_grid → 원래 함수를 부른 뒤 붙잡은 세계의 창 계수를 곁 파일 목록에 넣는다.
  4. --design rev : L_SHORT · L_LONG 을 맞바꾼다 → 배경 = 짧은 고리(L2) 상주 · 끼우는 'mut' = 긴 고리(L4).
     mini27 코드 경로는 같고 고리 상수만 바뀐다. 🚩 동일성 가드 밖이다.
  5. --wave A : 양성 대조 전용. assay 세계의 μ 를 --wave-block 틱 덩어리마다 μ̄(1 + A·s), s = ±1 로 흔든다.
     부호는 (seed, ρ, 'galA-wave') 로 키한 따로 흐름이라 세계의 난수 흐름과 섞이지 않는다.

사용 (petri/ 안에서)
  python3 galA_window.py --design std job --rho 8 --seed 41001 --nref 0 --mus 0,0.001625,0.00325,0.0065,0.01625,0.025,0.04 --out galA_std
  python3 galA_window.py --design std job --rho 8 --seed 27001 --out ga1_r8          (G-A1 · mini27 기본 설정 그대로)
출력: --out 안에 mini27 의 본 JSON(mini27_r{ρ}_s{seed}.json · 그대로) + 곁 파일 galA_r{ρ}_s{seed}.win.json
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np

import mini27 as m27

VERSION = "galA_window v1.0.0"
NC = len(m27.CNAMES)

WBIN = 50
WAVE_AMP = 0.0
WAVE_BLOCK = 500
WAVE_SIGNS = None
_LAST = {}
SIDE = {"grow": {}, "arms": []}


class _Tee:
    """c[i] += v 를 동결 합계 목록 a 와 창 목록 b 에 같이 올린다(act() 는 이 꼴로만 계수기를 쓴다)."""
    __slots__ = ("a", "b")

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __getitem__(self, i):
        return self.a[i]

    def __setitem__(self, i, v):
        self.b[i] += v - self.a[i]
        self.a[i] = v


def _counter_tee(self, o):
    key = (o.lin, o.L)
    c = self.cnt.get(key)
    if c is None:
        c = self.cnt[key] = [0] * NC
    t = self.t
    wc = self.__dict__.get("_wc")
    if wc is None:
        wc = self._wc = {}
    wk = (o.lin, o.L, (t - 1) // WBIN, t & 1)
    w = wc.get(wk)
    if w is None:
        w = wc[wk] = [0] * NC
    return _Tee(c, w)


_orig_tick = m27.World.tick


def _tick_wave(self):
    """양성 대조 — 첫 틱에 본 μ 를 μ̄ 로 잡고 덩어리마다 μ̄(1 + A·s) 로 바꾼다. μ̄ = 0 인 세계(성장 단계 등)는 그대로."""
    mb = self.__dict__.get("_mubar")
    if mb is None:
        mb = self._mubar = self.mu
    if mb > 0:
        k = self.t // WAVE_BLOCK                    # 다음 틱 t+1 이 속한 덩어리 = (t+1-1)//블록
        self.mu = mb * (1.0 + WAVE_AMP * WAVE_SIGNS[k])
    _orig_tick(self)


def _ser(w):
    """창 계수 → {"lin_L": [[칸, 홀짝, 계수 10 개 …], …]} (칸 · 홀짝 순)"""
    out = {}
    for (lin, L, b, p), v in sorted(w.__dict__.get("_wc", {}).items()):
        out.setdefault("%d_%d" % (lin, L), []).append([b, p] + v)
    return out


def _wrap_build(orig):
    def f(*a, **k):
        w = orig(*a, **k)
        _LAST["w"] = w
        return w
    return f


def _wrap_grow(orig, geom_of):
    def f(*a, **k):
        res = orig(*a, **k)
        SIDE["grow"][geom_of(a)] = _ser(res[0])
        return res
    return f


def _wrap_arm(orig, geom_of):
    def f(*a, **k):
        _LAST.pop("w", None)
        r = orig(*a, **k)
        w = _LAST.pop("w")
        SIDE["arms"].append({"geom": geom_of(a), "rule": r["rule"], "mu": r["mu"], "kind": r["kind"],
                             "L_inj": r["L_inj"], "win": _ser(w)})
        return r
    return f


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    global WBIN, WAVE_AMP, WAVE_BLOCK, WAVE_SIGNS
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--design", choices=["std", "rev"], required=True)
    ap.add_argument("--wbin", type=int, default=50)
    ap.add_argument("--wave", type=float, default=0.0, help="양성 대조 진폭 A (0 = 끔)")
    ap.add_argument("--wave-block", type=int, default=500)
    w_args, rest = ap.parse_known_args()
    pk = argparse.ArgumentParser(add_help=False)
    pk.add_argument("--rho", type=float, required=True)
    pk.add_argument("--seed", type=int, required=True)
    pk.add_argument("--assay", type=int, default=3000)
    pk.add_argument("--out", required=True)
    k_args, _ = pk.parse_known_args(rest)
    assert w_args.wbin > 0 and w_args.wave_block > 0 and 0.0 <= w_args.wave <= 1.0
    WBIN = w_args.wbin

    m27.World.counter = _counter_tee
    m27.build = _wrap_build(m27.build)
    m27.build_grid = _wrap_build(m27.build_grid)
    m27.grow = _wrap_grow(m27.grow, lambda a: "glob")
    m27.grow_grid = _wrap_grow(m27.grow_grid, lambda a: a[7])
    m27.run_arm = _wrap_arm(m27.run_arm, lambda a: "glob")
    m27.run_arm_grid = _wrap_arm(m27.run_arm_grid, lambda a: a[4])
    if w_args.design == "rev":
        m27.L_SHORT, m27.L_LONG = 4, 2
    if w_args.wave > 0:
        WAVE_AMP, WAVE_BLOCK = w_args.wave, w_args.wave_block
        n = k_args.assay // WAVE_BLOCK + 2
        WAVE_SIGNS = m27.stream(k_args.seed, k_args.rho, "galA-wave").choice([-1.0, 1.0], size=n)
        m27.World.tick = _tick_wave

    sys.argv = ["mini27.py"] + rest
    m27.main()

    side = {"version": VERSION, "mini27_version": m27.VERSION, "mini27_sha256": sha_file(m27.__file__),
            "python": sys.version.split()[0], "numpy": np.__version__,
            "design": w_args.design, "L_res": m27.L_LONG, "L_inj_mut": m27.L_SHORT,
            "wbin": WBIN, "wave": w_args.wave, "wave_block": w_args.wave_block,
            "wave_signs": [] if WAVE_SIGNS is None else [int(s) for s in WAVE_SIGNS],
            "argv": rest, "cnames": m27.CNAMES, "grow": SIDE["grow"], "arms": SIDE["arms"]}
    fn = os.path.join(k_args.out, "galA_r%g_s%d.win.json" % (k_args.rho, k_args.seed))
    with open(fn + ".tmp", "w") as f:
        json.dump(side, f, separators=(",", ":"))
    os.replace(fn + ".tmp", fn)


if __name__ == "__main__":
    sys.exit(main())
