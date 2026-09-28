#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""축소 모형 D — 잘 섞인 세계(공간 없음). 설계 축소모형-D-설계-2026-09-27.md §2.
접시(sim.js)와 같은 것: 글자 보존 · 밀도 · 수명 · 코드 · 고리 길이 L · 주사위 규칙(굴린 뒤 재료) · 오류 셋(빠뜨림 μ/3 · 끼어듦 μ/3 · 바뀜 μ · 글자 11종).
없는 것: 공간(도달 범위 — 재료는 전역 풀) · 먹기/빌려쓰기/빛 · 우주선 돌연변이.
개체 주기: 자리잡기(빈 칸 있어야 · 1틱) → 시도마다 L 틱 → 글자 ℓ 개 다 쓰면 나누기(1틱) → 다시.
  (접시의 고리 밖 명령 몫은 '출생당 준비 틱' pre = ℓ − L 로 뭉뚱그린다)
재는 것(해부 14 와 같은 정의 · 끼운 계통만): R = (ok+stall+del)/ok · T = 복사 중 틱/ok · S = (R−1)·L · U = 부모 그대로인 개체의 자식 중 다른 코드 몫 · ΔΔ(μ).
실행: python3 reduced_d.py --cells 2000 --grow 10000 --seeds 3 --codes-of racld,rascld,rsacld,racldx --out _RESULT_reduced_d.json
"""
import argparse
import json
import math
import random
import sys
import time

LET = "nsracldjehx"; K = len(LET)
MINLEN, MAXLEN, AGE0, AGEVAR, RDEATH = 2, 48, 300, 300, 0.0003


def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None: out.setdefault(k, L)
    return out


enc = lambda s: [LET.index(ch) for ch in s]
dec = lambda g: "".join(LET[k] for k in g)


class Org:
    __slots__ = ("g", "L", "pre", "age", "maxAge", "child", "rh", "wait", "mem", "tag", "exact", "phase")

    def __init__(self, g, rnd, tag=0, exact=True):
        self.g = g; s = dec(g); L = loop_ticks(s)
        self.L = L if L else len(g)           # 고리 없는 코드(c 나 l 없음)는 복제 못 한다 → 아래 act 에서 idle
        self.pre = max(len(g) - self.L, 0)
        self.age = 0; self.maxAge = AGE0 + int(rnd() * AGEVAR)
        self.child = None; self.rh = 0; self.wait = 0; self.mem = None; self.tag = tag; self.exact = exact; self.phase = 0


class World:
    def __init__(self, cells, density, mu, seed, rule="orig"):
        self.cells = cells; self.mu = mu; self.rule = rule; self.rnd = random.Random(seed).random
        tot = density * cells; self.pool = [tot // K] * K
        for k in range(tot - sum(self.pool)): self.pool[k % K] += 1
        self.orgs = []; self.tick_n = 0
        self.cnt = {}

    def free_cells(self):
        return self.cells - len(self.orgs) - sum(1 for o in self.orgs if o.child is not None)

    def seed_with(self, code, n):
        g = enc(code)
        for _ in range(n):
            if all(self.pool[k] > 0 for k in g):
                for k in g: self.pool[k] -= 1
                self.orgs.append(Org(g, self.rnd))

    def count(self, o, key):
        c = self.cnt.setdefault(o.tag, {}); c[key] = c.get(key, 0) + 1

    def die(self, o):
        for k in o.g: self.pool[k] += 1
        if o.child is not None:
            for k in o.child: self.pool[k] += 1
        self.count(o, "deaths")

    def step(self, o, births):
        can_copy = "c" in dec(o.g) and "l" in dec(o.g)
        if o.child is None:
            if not can_copy: return
            if self.free_cells() > 0:
                o.child = []; o.rh = 0; o.wait = o.pre; o.mem = None; self.count(o, "alloc_ok")
            else:
                self.count(o, "alloc_blocked")
            return
        if o.rh < len(o.g):
            self.count(o, "copy_phase")
            if o.wait > 0: o.wait -= 1; return
            o.wait = o.L - 1
            # 시도
            mu = self.mu
            if self.rule == "rem" and o.mem is not None:
                k, adv, dtype = o.mem
            else:
                r = self.rnd()
                if r < mu / 3:
                    o.rh += 1; o.mem = None; self.count(o, "copy_del"); return
                k = o.g[o.rh]; adv = True; dtype = 0
                if r < mu * 2 / 3: k = int(self.rnd() * K); adv = False; dtype = 1
                elif r < mu * 5 / 3: k = int(self.rnd() * K); dtype = 2
                if self.rule == "rem": o.mem = (k, adv, dtype)
            if self.pool[k] > 0:
                self.pool[k] -= 1; o.child.append(k); o.mem = None
                if adv: o.rh += 1
                self.count(o, "copy_ok")
                if len(o.child) > MAXLEN:
                    for kk in o.child: self.pool[kk] += 1
                    o.child = None
            else:
                self.count(o, "copy_stall")
            return
        # 나누기
        g = o.child; o.child = None
        if len(g) < MINLEN:
            for k in g: self.pool[k] += 1
            return
        kid = Org(g, self.rnd, tag=o.tag, exact=(g == o.g and o.exact))
        births.append((o, kid))

    def tick(self):
        self.tick_n += 1
        order = self.orgs[:]; self.rnd_shuffle(order)
        births = []; dead = set()
        for o in order:
            o.age += 1
            if o.age > o.maxAge or self.rnd() < RDEATH:
                self.die(o); dead.add(id(o)); continue
            self.step(o, births)
        if dead: self.orgs = [o for o in self.orgs if id(o) not in dead]
        for parent, kid in births:
            self.orgs.append(kid); self.count(parent, "births")
            if parent.exact:
                c = self.cnt.setdefault(parent.tag, {}); c["b_exact"] = c.get("b_exact", 0) + 1
                if kid.g != parent.g: c["kids_diff"] = c.get("kids_diff", 0) + 1

    def rnd_shuffle(self, a):
        r = self.rnd
        for i in range(len(a) - 1, 0, -1):
            j = int(r() * (i + 1)); a[i], a[j] = a[j], a[i]


def clone(w, seed):
    import copy
    w2 = World.__new__(World)
    w2.cells = w.cells; w2.mu = w.mu; w2.rule = w.rule; w2.rnd = random.Random(seed).random
    w2.pool = w.pool[:]; w2.tick_n = w.tick_n; w2.cnt = {}
    w2.orgs = []
    for o in w.orgs:
        n = Org.__new__(Org)
        for f in Org.__slots__: setattr(n, f, getattr(o, f))
        n.g = o.g[:]; n.child = o.child[:] if o.child is not None else None
        w2.orgs.append(n)
    return w2


def inject(w, code, frac, rnd):
    g = enc(code); idx = list(range(len(w.orgs))); rnd_shuffle(idx, rnd)
    target = round(len(w.orgs) * frac); done = 0
    for i in idx[:target]:
        o = w.orgs[i]
        if o.child is not None:
            for k in o.child: w.pool[k] += 1
            o.child = None
        for k in o.g: w.pool[k] += 1
        if all(w.pool[k] > 0 for k in g):
            for k in g: w.pool[k] -= 1
            o.g = g[:]; o.L = loop_ticks(code) or len(g); o.pre = max(len(g) - o.L, 0); o.rh = 0; o.wait = 0; o.mem = None; o.tag = 1; o.exact = True
            done += 1
        else:
            for k in o.g: w.pool[k] -= 1
    return done


def rnd_shuffle(a, r):
    for i in range(len(a) - 1, 0, -1):
        j = int(r() * (i + 1)); a[i], a[j] = a[j], a[i]


def s_of(n0, m0, nT, mT):
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=2000); ap.add_argument("--density", type=int, default=8)
    ap.add_argument("--grow", type=int, default=10000); ap.add_argument("--ticks", type=int, default=3000)
    ap.add_argument("--seeds", type=int, default=3); ap.add_argument("--seed0", type=int, default=101)
    ap.add_argument("--codes-of", default="racld,rascld,rsacld,racldx"); ap.add_argument("--mus", default="0,0.003")
    ap.add_argument("--rules", default="orig,rem"); ap.add_argument("--resident", default="racld"); ap.add_argument("--frac", type=float, default=0.1)
    ap.add_argument("--out", default="_RESULT_reduced_d.json")
    a = ap.parse_args()
    bases = a.codes_of.split(","); mus = [float(x) for x in a.mus.split(",")]; rules = a.rules.split(",")
    codes = {}
    for b in bases: codes.update(inserts(b, "n"))
    OUT = {"args": vars(a), "cells": [], "codes": {c: L for c, L in codes.items()}}
    t0 = time.time()
    for si in range(a.seeds):
        seed = a.seed0 + si
        w = World(a.cells, a.density, 0.01, seed)
        w.seed_with(a.resident, max(1, a.cells // 20))
        for _ in range(a.grow): w.tick()
        n0 = len(w.orgs); free = sum(w.pool); used_kinds = set(enc(a.resident))
        free_used = sum(w.pool[k] for k in used_kinds); tot_used = a.density * a.cells * len(used_kinds) // K
        print("seed %d · grow %d · pop %d · free letters %d · consumed kinds free %.1f%% · %.0fs" % (seed, a.grow, n0, free, 100 * free_used / tot_used, time.time() - t0), flush=True)
        OUT["cells"].append({"seed": seed, "pop": n0, "free": free, "free_used_frac": free_used / tot_used})
        for rule in rules:
            for mu in mus:
                for c, L in codes.items():
                    w2 = clone(w, seed * 7919 + 1); w2.mu = mu; w2.rule = rule
                    inj = inject(w2, c, a.frac, random.Random(seed * 1000003 + 17).random)
                    m0 = sum(1 for o in w2.orgs if o.tag == 1); n0b = len(w2.orgs)
                    for _ in range(a.ticks): w2.tick()
                    mT = sum(1 for o in w2.orgs if o.tag == 1); nT = len(w2.orgs)
                    cnt = w2.cnt.get(1, {}); ok = cnt.get("copy_ok", 0)
                    R = (ok + cnt.get("copy_stall", 0) + cnt.get("copy_del", 0)) / ok if ok else float("nan")
                    T = cnt.get("copy_phase", 0) / ok if ok else float("nan")
                    U = cnt.get("kids_diff", 0) / cnt["b_exact"] if cnt.get("b_exact") else float("nan")
                    s = s_of(n0b, m0, nT, mT) if m0 else float("nan")
                    OUT.setdefault("arm", []).append({"seed": seed, "rule": rule, "mu": mu, "code": c, "L": L, "inj": inj, "R": R, "T": T, "S": (R - 1) * L, "U": U, "s": s, "m0": m0, "mT": mT})
                print("  rule %s μ %g · %d 코드 · %.0fs" % (rule, mu, len(codes), time.time() - t0), flush=True)
    json.dump(OUT, open(a.out, "w"), indent=1)
    print("기록 →", a.out)


if __name__ == "__main__":
    main()
