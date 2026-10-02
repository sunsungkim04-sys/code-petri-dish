#!/usr/bin/env python3
"""mini26.py — 해부 26: 두 가정만 담은 최소 모형 (독립 재구현)

이 파일은 원 시뮬레이터(sim.js)를 열지도, import 하지도, 번역하지도 않고 썼다.
사양의 출처는 원고 말뿐이다 — report/S2_Model.md Table 1 · Table 4 · report/S5_Methods.md §5.1.
축소 세계(reduced_d.py · 해부 22)도 열지 않았다. 침입 지표의 꼴(s · Δ · 교차 보간)만
판정기 inv21_analyze.py 의 정의를 따른다(그건 mini26_analyze.py 에 있다).

담은 가정 둘
  (가정1) 재료 보존 · 멈춤: 글자는 몸 · 만드는 중인 자식 · 자유 풀 사이를 오갈 뿐 총량이 상수다.
          복사 시도 때 필요한 글자가 자유 풀에 없으면 그 시도는 멈춤(stall) — 틱을 쓰고 아무것도 안 쓴다.
  (가정2) roll-first: 오류 추첨이 재료 탐색보다 먼저 온다 → 멈춘 시도도 추첨 하나를 쓴다.
          대조 규칙 find-first: 필요한 글자가 있을 때만 추첨한다.

일부러 단순하게 한 것 (원고 세계와 다른 점 — 정량 일치를 목표로 하지 않는다)
  · 공간 없음(well-mixed): 칸 K 개는 '자리 수' 로만 있다. 재료는 하나의 자유 풀 — 이웃 · 확산 없음.
    멈춤은 '그 글자 종류가 풀 전체에서 0' 일 때만 일어난다(국소 결핍이 아니라 전역 결핍).
  · 두 종은 같은 글자열(6 글자)을 갖고 고리 길이 L(틱/시도)만 다르다. L 은 글자가 아니라 계통이 물려주는 꼬리표다.
  · 오류가 하나라도 실현된 자식은 불임(복제 안 함 · 자리와 재료만 차지)이다.
  · 죽음은 나이와 무관한 상수율 하나.
  · 우주선(배경 돌연변이) · 먹기 · 빛 · 이웃 실행 없음.

틱 의미
  모든 산 개체는 매 틱 한 명령을 쓴다. 행동 시점만 일정에 올린다:
    claim(1 틱): 빈 자리가 있으면(점유 < K) 하나를 자식 자리로 잡는다 · 없으면 다음 틱 다시.
    copy 시도(L 틱): 시도 하나가 L 틱을 쓴다(멈춰도 L 틱).
    release(1 틱): 자식을 독립 개체로 놓는다.
  같은 틱에 행동하는 개체는 매 틱 새로 섞은 순서로 처리한다(먼저 온 쪽이 풀의 글자 · 빈 자리를 가져간다).
  그 뒤 모든 산 개체가 상수율 d 로 죽는다. 죽으면 몸 · 만드는 중인 자식의 글자가 풀로, 자리가 비워진다.

오류 추첨(Table 1 의 말): 글자당 μ, 셋으로 같게 — 빠뜨림 · 무작위 글자 끼우기 · 무작위 글자로 바꾸기.
  무작위 글자는 11 종에서 고르게(1/11 은 원래 글자와 같다).
  빠뜨림: 재료가 필요 없다 · 위치를 넘긴다.
  끼우기: 글자 y 가 필요 · 위치를 안 넘긴다.        y 가 없으면 멈춤(걸러진 오류)
  바꾸기: 글자 y 가 필요 · 위치를 넘긴다.           y 가 없으면 멈춤(걸러진 오류)
  옳음  : 글자 x(원래 글자) 가 필요 · 위치를 넘긴다.  x 가 없으면 멈춤

난수: (seed, 이름) 으로 키한 독립 흐름(sha256 → SeedSequence). 흐름 이름에 규칙(roll/find)을 넣지 않는다 —
  그래서 μ = 0 에서 두 규칙의 궤적이 같아야 한다(추첨 흐름은 따로라 소비 수가 달라도 다른 흐름에 안 닿는다).
  이것이 G0 의 '규칙 동일성' 가드다.

사용
  python3 mini26.py job --rho 8 --seed 26001 --out mini26_out
  python3 mini26.py job --rho 8 --seed 26901 --mus 0,0.01 --nref 2 --out pilot26   (파일럿)
"""
import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

VERSION = "mini26 v1.0.0"
NKIND = 11                      # 글자 종류 수 (Table 1: 11-letter alphabet)
CODE = (0, 1, 2, 3, 4, 5)       # 두 종 공통 글자열 — 6 글자 · 종류 0~5 (나머지 5 종은 아무도 안 쓴다)
L_SHORT, L_LONG = 2, 4          # 고리 길이: 창시자 꼴 · 상주 꼴 (Table 3: 2 · 4 ticks per letter)
MAXLEN = 48                     # 자식 글자열 상한 (Table 1: 2–48) — 넘으면 자식을 버린다
MUS_DEFAULT = (0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01)   # = petri/inv21_launch.sh:20
RULES = ("roll", "find")

# 행동 단계
IDLE, CLAIM, COPY, RELEASE = 0, 1, 2, 3
# 계수기 칸
C_ATT, C_DRAW, C_POS, C_TICK, C_STALL, C_FILT, C_ERR, C_KOK, C_KBAD, C_ABORT = range(10)
CNAMES = ["attempts", "draws", "positions", "copy_ticks", "stalls", "filtered", "err_realized",
          "kids_ok", "kids_sterile", "aborted"]


def stream(*parts):
    """(seed, 이름 …) 으로 키한 독립 난수 흐름. python hash() 는 쓰지 않는다(실행마다 바뀜)."""
    key = "|".join(str(p) for p in parts).encode()
    h = int.from_bytes(hashlib.sha256(key).digest()[:16], "little")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(h)))


class Buf:
    """흐름에서 균등 난수를 덩어리로 미리 뽑아 하나씩 준다(소비 순서가 같으면 값도 같다)."""
    __slots__ = ("g", "a", "i")

    def __init__(self, g):
        self.g = g
        self.a = g.random(4096)
        self.i = 0

    def next(self):
        if self.i == 4096:
            self.a = self.g.random(4096)
            self.i = 0
        v = self.a[self.i]
        self.i += 1
        return v


class Org:
    __slots__ = ("code", "L", "lin", "phase", "ptr", "buf", "slot", "k", "nxt", "err")

    def __init__(self, code, L, lin):
        self.code = code
        self.L = L          # 0 = 불임
        self.lin = lin      # 0 = 배경 · 1 = 끼운 계통(자손 포함 · 혈통으로 따라간다)
        self.phase = CLAIM if L > 0 else IDLE
        self.ptr = 0
        self.buf = []
        self.slot = False
        self.k = -1         # alive 목록 안 자리
        self.nxt = -1       # 다음 행동 틱
        self.err = False    # 만드는 중인 자식에 실현된 오류가 있나


class World:
    def __init__(self, K, M, rho, death, mu, rule, rng_sched, rng_death, rng_draw):
        self.K = K
        self.M = M
        self.rho = rho
        self.d = death
        self.mu = mu
        self.rule = rule
        self.sched = rng_sched
        self.dg = rng_death
        self.draw = Buf(rng_draw)
        self.free = [0] * NKIND
        self.occ = 0
        self.alive = []
        self.due = {}
        self.t = 0
        self.cnt = {}
        self.cons_checks = 0
        self.cons_bad = []

    # ── 재료 ──
    def take(self, x):
        if self.free[x] > 0:
            self.free[x] -= 1
            return True
        return False

    def give(self, letters):
        f = self.free
        for x in letters:
            f[x] += 1

    # ── 개체 목록 ──
    def add(self, o, when):
        o.k = len(self.alive)
        self.alive.append(o)
        if o.phase != IDLE:
            self.schedule(o, when)

    def schedule(self, o, when):
        o.nxt = when
        self.due.setdefault(when, []).append(o)

    def kill_at(self, k):
        o = self.alive[k]
        self.give(o.code)
        self.give(o.buf)
        self.occ -= 1 + (1 if o.slot else 0)
        o.buf = []
        o.slot = False
        o.k = -1
        last = self.alive.pop()
        if last is not o:
            self.alive[k] = last
            last.k = k

    def counter(self, o):
        key = (o.lin, o.L)
        c = self.cnt.get(key)
        if c is None:
            c = self.cnt[key] = [0] * len(CNAMES)
        return c

    # ── 보존 재계수(계수기를 믿지 않고 구조에서 다시 센다) ──
    def recount(self):
        tot = sum(self.free)
        occ = 0
        for o in self.alive:
            tot += len(o.code) + len(o.buf)
            occ += 1 + (1 if o.slot else 0)
            if (not o.slot) and o.buf:
                self.cons_bad.append("t%d 자리 없는 자식 글자" % self.t)
        self.cons_checks += 1
        if tot != self.M:
            self.cons_bad.append("t%d 글자 합 %d != %d" % (self.t, tot, self.M))
        if occ != self.occ or occ > self.K:
            self.cons_bad.append("t%d 점유 %d / 계수 %d / K %d" % (self.t, occ, self.occ, self.K))
        if min(self.free) < 0:
            self.cons_bad.append("t%d 음수 자유 글자" % self.t)

    # ── 한 행동 ──
    def act(self, o):
        t = self.t
        ph = o.phase
        if ph == CLAIM:
            if self.occ < self.K:
                self.occ += 1
                o.slot = True
                o.phase = COPY
                o.ptr = 0
                o.buf = []
                o.err = False
            self.schedule(o, t + 1)
            return
        if ph == RELEASE:
            ok = (not o.err) and tuple(o.buf) == o.code
            c = self.counter(o)
            child = Org(tuple(o.buf), o.L if ok else 0, o.lin)
            c[C_KOK if ok else C_KBAD] += 1
            child.slot = False          # 자식 자리는 이제 자식 자신의 칸(점유 수는 그대로)
            o.slot = False
            o.buf = []
            o.phase = CLAIM
            self.add(child, t + 1)
            self.schedule(o, t + 1)
            return
        # COPY 시도
        c = self.counter(o)
        L = o.L
        c[C_ATT] += 1
        c[C_TICK] += L
        x = o.code[o.ptr]
        mu = self.mu
        if self.rule == "find" and self.free[x] == 0:
            c[C_STALL] += 1                    # 찾기 먼저: 없으면 추첨 없이 멈춤
            self.schedule(o, t + L)
            return
        u = self.draw.next()
        c[C_DRAW] += 1
        if u < mu / 3.0:                       # 빠뜨림
            o.ptr += 1
            o.err = True
            c[C_POS] += 1
            c[C_ERR] += 1
        elif u < mu:                           # 끼우기(mu/3~2mu/3) · 바꾸기(2mu/3~mu)
            ins = u < 2.0 * mu / 3.0
            y = int(self.draw.next() * NKIND)
            if self.take(y):
                o.buf.append(y)
                if ins:
                    o.err = True
                    c[C_ERR] += 1
                else:
                    o.ptr += 1
                    c[C_POS] += 1
                    if y != x:
                        o.err = True
                        c[C_ERR] += 1
            else:
                c[C_STALL] += 1
                c[C_FILT] += 1
        else:                                  # 옳은 복사
            if self.take(x):
                o.buf.append(x)
                o.ptr += 1
                c[C_POS] += 1
            else:
                c[C_STALL] += 1
        if len(o.buf) > MAXLEN:                # 끼우기가 쌓여 상한을 넘으면 자식을 버린다
            self.give(o.buf)
            o.buf = []
            o.slot = False
            self.occ -= 1
            o.phase = CLAIM
            c[C_ABORT] += 1
            self.schedule(o, t + 1)
            return
        if o.ptr >= len(o.code):
            o.phase = RELEASE
        self.schedule(o, t + L)

    def tick(self):
        self.t += 1
        t = self.t
        lst = self.due.pop(t, None)
        if lst:
            lst = [o for o in lst if o.k >= 0 and o.nxt == t]
            if len(lst) > 1:
                order = self.sched.permutation(len(lst))
                for i in order:
                    self.act(lst[i])
            elif lst:
                self.act(lst[0])
        n = len(self.alive)
        if n:
            dead = np.flatnonzero(self.dg.random(n) < self.d)
            for k in dead[::-1]:               # 뒤에서부터 지워야 자리 바꾸기가 안 섞인다
                self.kill_at(int(k))

    # ── 상태 저장 · 복원(배경을 팔마다 복제) ──
    def snapshot(self):
        orgs = []
        for o in self.alive:
            orgs.append([list(o.code), o.L, o.lin, o.phase, o.ptr, list(o.buf), o.slot,
                         (o.nxt - self.t) if o.phase != IDLE else -1, o.err])
        return {"free": list(self.free), "occ": self.occ, "orgs": orgs}


def snap_checksum(snap):
    return hashlib.sha256(json.dumps(snap, sort_keys=True).encode()).hexdigest()[:16]


def build(snap, K, M, rho, death, mu, rule, rs, rd, rw):
    w = World(K, M, rho, death, mu, rule, rs, rd, rw)
    w.free = list(snap["free"])
    w.occ = snap["occ"]
    for code, L, lin, ph, ptr, buf, slot, off, err in snap["orgs"]:
        o = Org(tuple(code), L, lin)
        o.phase, o.ptr, o.buf, o.slot, o.err = ph, ptr, list(buf), slot, err
        w.add(o, off if ph != IDLE else 0)
    return w


def census(w):
    n = len(w.alive)
    m = mf = nf = 0
    for o in w.alive:
        if o.L > 0:
            nf += 1
        if o.lin == 1:
            m += 1
            if o.L > 0:
                mf += 1
    return [w.t, n, m, mf, nf]


def grow(rho, seed, K, death, T, n0, every):
    M = int(round(rho * K))
    w = World(K, M, rho, death, 0.0, "roll",
              stream(seed, rho, "grow", "sched"), stream(seed, rho, "grow", "death"),
              stream(seed, rho, "grow", "draw"))
    base, rem = divmod(M, NKIND)
    w.free = [base + (1 if i < rem else 0) for i in range(NKIND)]   # 11 종에 고르게(결정적)
    for _ in range(n0):
        if not all(w.free[x] > 0 for x in CODE) or w.occ >= K:
            break
        for x in CODE:
            w.take(x)
        o = Org(CODE, L_LONG, 0)
        w.occ += 1
        w.add(o, 1)
    series = [census(w)]
    w.recount()
    ext = -1
    for _ in range(T):
        w.tick()
        if w.t % every == 0:
            w.recount()
            series.append(census(w))
        if not w.alive and ext < 0:
            ext = w.t
            break
    w.recount()
    return w, M, series, ext


def inject(snap, seed, rho, code_L, frac):
    """배경 사본의 산 개체 frac 를 끼운 코드로 바꾼다 — 재료 보존(돌려준 글자 → 새 글자를 풀에서). 못 하면 건너뛰고 센다.
    고를 개체는 (seed, rho, 'inject') 흐름으로 — 모든 팔 · μ · 규칙에서 같은 개체가 바뀐다."""
    s = json.loads(json.dumps(snap))
    orgs = s["orgs"]
    free = s["free"]
    k = int(round(frac * len(orgs)))
    pick = stream(seed, rho, "inject").permutation(len(orgs))[:k]
    skipped = 0
    for i in sorted(int(j) for j in pick):
        code, L, lin, ph, ptr, buf, slot, off, err = orgs[i]
        back = list(code) + list(buf)
        for x in back:
            free[x] += 1
        need = {}
        for x in CODE:
            need[x] = need.get(x, 0) + 1
        if all(free[x] >= n for x, n in need.items()):
            for x in CODE:
                free[x] -= 1
            if slot:
                s["occ"] -= 1
            orgs[i] = [list(CODE), code_L, 1, CLAIM, 0, [], False, 1, False]
        else:
            for x in back:                     # 되돌린다(방금 돌려준 글자라 반드시 있다)
                free[x] -= 1
            skipped += 1
    return s, k, skipped


def run_arm(snap, M, seed, rho, mu, rule, kind, code_L, K, death, T, every, frac):
    s, k, skipped = inject(snap, seed, rho, code_L, frac)
    tag = "mu%g" % mu
    w = build(s, K, M, rho, death, mu, rule,
              stream(seed, rho, tag, kind, "sched"), stream(seed, rho, tag, kind, "death"),
              stream(seed, rho, tag, kind, "draw"))
    w.t = 0
    w.recount()
    series = [census(w)]
    ext = -1
    for _ in range(T):
        w.tick()
        if w.t % every == 0:
            w.recount()
            series.append(census(w))
        if not w.alive:
            ext = w.t
            break
    w.recount()
    cnt = {"%d_%d" % key: dict(zip(CNAMES, v)) for key, v in sorted(w.cnt.items())}
    return {"rule": rule, "mu": mu, "kind": kind, "L_inj": code_L, "n_inject": k, "skipped": skipped,
            "series": series, "extinct_tick": ext, "counters": cnt,
            "cons_ok": not w.cons_bad, "cons_checks": w.cons_checks, "cons_bad": w.cons_bad[:5],
            "end_free": w.free}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["job"])
    ap.add_argument("--rho", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--K", type=int, default=1000)
    ap.add_argument("--death", type=float, default=1.0 / 450)
    ap.add_argument("--grow", type=int, default=6000)
    ap.add_argument("--assay", type=int, default=3000)
    ap.add_argument("--every", type=int, default=250)
    ap.add_argument("--n0", type=int, default=100)
    ap.add_argument("--frac", type=float, default=0.1)
    ap.add_argument("--nref", type=int, default=5)
    ap.add_argument("--min-pop", type=int, default=50)
    ap.add_argument("--mus", default=",".join("%g" % m for m in MUS_DEFAULT))
    ap.add_argument("--rules", default="roll,find")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mus = [float(x) for x in a.mus.split(",")]
    rules = a.rules.split(",")
    assert all(r in RULES for r in rules)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    w, M, gseries, gext = grow(a.rho, a.seed, a.K, a.death, a.grow, a.n0, a.every)
    snap = w.snapshot()
    gcnt = {"%d_%d" % key: dict(zip(CNAMES, v)) for key, v in sorted(w.cnt.items())}
    n_end = len(w.alive)
    qualified = gext < 0 and n_end >= a.min_pop
    out = {"version": VERSION, "config": {k: getattr(a, k) for k in
                                          ("rho", "seed", "K", "death", "grow", "assay", "every", "n0",
                                           "frac", "nref", "min_pop")},
           "mus": mus, "rules": rules, "M": M,
           "grow": {"series": gseries, "extinct_tick": gext, "n_end": n_end, "qualified": qualified,
                    "checksum": snap_checksum(snap), "counters": gcnt, "cons_ok": not w.cons_bad,
                    "cons_checks": w.cons_checks, "cons_bad": w.cons_bad[:5], "end_free": w.free,
                    "occ_end": w.occ},
           "arms": []}
    if qualified:
        kinds = [("mut", L_SHORT)] + [("ref%d" % i, L_LONG) for i in range(a.nref)]
        for rule in rules:
            for mu in mus:
                for kind, L in kinds:
                    out["arms"].append(run_arm(snap, M, a.seed, a.rho, mu, rule, kind, L, a.K, a.death,
                                               a.assay, a.every, a.frac))
    out["seconds"] = round(time.time() - t0, 1)
    fn = os.path.join(a.out, "mini26_r%g_s%d.json" % (a.rho, a.seed))
    with open(fn + ".tmp", "w") as f:
        json.dump(out, f)
    os.replace(fn + ".tmp", fn)
    print("%s r%g s%d qualified=%s n_end=%d arms=%d %.1fs" % (VERSION, a.rho, a.seed, qualified, n_end,
                                                            len(out["arms"]), out["seconds"]))


if __name__ == "__main__":
    sys.exit(main())
