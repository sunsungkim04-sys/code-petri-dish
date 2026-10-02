#!/usr/bin/env python3
"""mini27.py — 해부 27: 해부 26 최소 모형(mini26.py)에 '국소 결핍(공간)' 만 더한 판

출발점: mini26.py v1.0.0(7d9b9901…)을 복사해 고쳤다. sim.js 의 번역이 아니다 — sim.js 는 읽었지만
(해부 26 독립성은 끝남) 가져온 것은 '이웃 다섯 칸에서 재료' 라는 말뿐이고, 명령 실행 · 방향 · 확산 · 수명 분포 ·
오류 분할은 가져오지 않았다. mini26 에서 바뀐 것은 재료 풀(과 그 풀이 놓일 자리)뿐이다:

  기하(geom) 셋 — 같은 동역학 코드(World.act)를 쓰고 '자리 잡기 · 재료 찾기 · 돌려주기' 고리(hook)만 다르다.
    glob  : 자리 전역 · 재료 전역  = mini26 그대로(코드 경로 · 난수 흐름 이름까지 같다 → 같은 시드면 mini26 출력과 바이트 같음)
    loc   : 자리 국소 · 재료 국소  — NX × NY 주기 격자(40 × 25 = 1,000 칸 = mini26 의 K) · 생물 · 만드는 중인 자식은 한 칸씩
            · 자식 자리 = 비어 있는 4 이웃 칸(von Neumann) 중 하나(고르게 · 'place' 흐름)
            · 복사 시도 = 자기 칸 + 반경 r(맨해튼 · 기본 1 = 자기 칸 + 4 이웃) 안 자유 글자에서만 찾는다(자기 칸 먼저 · 그다음 고정 순서)
            · 죽으면 몸 글자는 자기 칸에, 만드는 중인 자식 글자는 자식 칸에 돌려준다(MAXLEN 초과로 버린 자식도 자식 칸)
            · 처음 재료는 칸마다 ρ 글자 · 칸 c 의 글자 종류 = (ρc + j) mod 11, j < ρ — 종류별 총량이 mini26 의 divmod 분할과 정확히 같다
    place : 자리 국소 · 재료 전역  (서술 대조 — '공간 배치만' 의 효과를 '재료 국소성' 과 가르려고)
  반경 r ≥ 격자 지름(NX//2 + NY//2 = 32)이면 재료 이웃이 격자 전체 = 전역 풀과 같은 뜻이다(--mat-radius all).

mini26 에서 바꾸지 않은 것 (한 요인만 바꾸기)
  오류 분할(μ/3 빠뜨림 · μ/3 끼우기 · μ/3 바꾸기 — 접시의 μ/3 · μ/3 · μ 로 고치지 않는다) · 불임 가정(실현 오류 하나 = 불임) ·
  사망률 d · 두 종 같은 글자열 · L 꼬리표(2 · 4) · claim 1 틱 · 시도 L 틱 · release 1 틱 · 흐름 키 방식 · 확산 없음.
  μ 사다리만 넓혔다: 원 사다리 0 … 1%(mini26 · inv21) + 확장 사다리 = 원 사다리 × 52/32 = 0 … 1.625%
  (접시 글자 변경 52μ/33 · mini 글자 변경 32μ/33 — report/S2_Model.md:95 · 해부26 사전등록 §12c).

난수: (seed, ρ, [geom], 단계, μ, 팔, 용도) 를 sha256 → SeedSequence 로 키한 독립 흐름. glob 은 mini26 의 이름 그대로(geom 없음).
  loc · place 는 이름에 geom 을 넣는다. 규칙(roll/find)은 어느 흐름 이름에도 없다 → μ = 0 두 규칙 궤적 동일(G0 가드).

사용
  python3 mini27.py job --rho 8 --seed 27001 --out mini27_out
  python3 mini27.py job --rho 8 --seed 26001 --geoms glob --mus 0,0.001,0.002,0.003,0.004,0.005,0.0075,0.01 --out ident27   (G2a)
"""
import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

VERSION = "mini27 v1.0.0"
NKIND = 11                      # 글자 종류 수 (Table 1: 11-letter alphabet)
CODE = (0, 1, 2, 3, 4, 5)       # 두 종 공통 글자열 — 6 글자 · 종류 0~5 (나머지 5 종은 아무도 안 쓴다)
L_SHORT, L_LONG = 2, 4          # 고리 길이: 창시자 꼴 · 상주 꼴
MAXLEN = 48                     # 자식 글자열 상한 — 넘으면 자식을 버린다
MUS_ORIG = (0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01)                       # = mini26 · inv21_launch.sh:20
MUS_EXT = (0.0, 0.001625, 0.00325, 0.004875, 0.0065, 0.008125, 0.0121875, 0.01625)      # = MUS_ORIG × 52/32 (접시 글자변경 척도 0–1%)
MUS_ALL = tuple(sorted(set(MUS_ORIG) | set(MUS_EXT)))                                   # 15 점
RULES = ("roll", "find")
NX, NY = 40, 25                 # 주기 격자 — 1,000 칸 = mini26 K
NCELL = NX * NY
GEOMS = ("glob", "loc", "place")

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


# ── 격자 ──
def _xy(c):
    return c % NX, c // NX


def nb4(c):
    x, y = _xy(c)
    return [((y - 1) % NY) * NX + x, y * NX + (x + 1) % NX, ((y + 1) % NY) * NX + x, y * NX + (x - 1) % NX]


def tdist(a, b):
    ax, ay = _xy(a)
    bx, by = _xy(b)
    dx = abs(ax - bx)
    dy = abs(ay - by)
    return min(dx, NX - dx) + min(dy, NY - dy)


def ball(c, r):
    """자기 칸 + 맨해튼 반경 r 안 칸 — 거리 순 · 같은 거리는 칸 번호 순(r = 1 이면 자기 칸 + 4 이웃)"""
    if r == 1:
        return [c] + nb4(c)
    cells = [(tdist(c, d), d) for d in range(NCELL)]
    return [d for dist, d in sorted(cells) if dist <= r]


NB4 = [nb4(c) for c in range(NCELL)]
DIAM = NX // 2 + NY // 2


class Org:
    __slots__ = ("code", "L", "lin", "phase", "ptr", "buf", "slot", "k", "nxt", "err", "pos", "cpos")

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
        self.pos = -1       # (loc · place) 자기 칸
        self.cpos = -1      # (loc · place) 자식 칸


class World:
    """mini26 의 World 그대로 — 재료 · 자리 연산만 고리(hook) 메서드로 뽑았다(glob 기하).
    고리로 뽑은 것이 동역학 · 난수 소비를 바꾸지 않았다는 것은 G2a(같은 시드에서 mini26 출력과 바이트 같음)로 본다."""

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

    # ── 재료 · 자리 고리 (glob: 전역 풀 · 자리 수만) ──
    def take(self, x):
        if self.free[x] > 0:
            self.free[x] -= 1
            return True
        return False

    def give(self, letters):
        f = self.free
        for x in letters:
            f[x] += 1

    def has_for(self, o, x):
        return self.free[x] > 0

    def take_for(self, o, x):
        return self.take(x)

    def claim_slot(self, o):
        if self.occ < self.K:
            self.occ += 1
            return True
        return False

    def drop_child(self, o):
        """MAXLEN 초과로 자식을 버린다 — 글자 · 자리 반납"""
        self.give(o.buf)
        self.occ -= 1

    def place_child(self, o, child):
        pass

    def free_org(self, o):
        """죽음 — 몸 · 만드는 중인 자식 글자 반납 · 자리 반납"""
        self.give(o.code)
        self.give(o.buf)
        self.occ -= 1 + (1 if o.slot else 0)

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
        self.free_org(o)
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

    # ── 한 행동 (mini26 과 같은 논리 · 재료 · 자리만 고리 호출) ──
    def act(self, o):
        t = self.t
        ph = o.phase
        if ph == CLAIM:
            if self.claim_slot(o):
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
            self.place_child(o, child)
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
        if self.rule == "find" and not self.has_for(o, x):
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
            if self.take_for(o, y):
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
            if self.take_for(o, x):
                o.buf.append(x)
                o.ptr += 1
                c[C_POS] += 1
            else:
                c[C_STALL] += 1
        if len(o.buf) > MAXLEN:                # 끼우기가 쌓여 상한을 넘으면 자식을 버린다
            self.drop_child(o)
            o.buf = []
            o.slot = False
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

    def free_total(self):
        return list(self.free)


class GridWorld(World):
    """loc · place 기하 — 자리는 격자 칸 · 자식은 빈 4 이웃 칸 · 재료는 (loc) 칸마다 / (place) 전역 풀.
    동역학(act · tick)은 World 그대로이고 고리만 바꾼다."""

    def __init__(self, K, M, rho, death, mu, rule, rng_sched, rng_death, rng_draw, rng_place, mat_r):
        World.__init__(self, K, M, rho, death, mu, rule, rng_sched, rng_death, rng_draw)
        assert K == NCELL
        self.place = Buf(rng_place)
        self.cell = [0] * NCELL                # 0 빈 칸 · 1 생물 또는 만드는 중인 자식
        self.local = mat_r is not None
        if self.local:
            self.fc = [[0] * NKIND for _ in range(NCELL)]
            self.ballof = [ball(c, mat_r) for c in range(NCELL)] if mat_r != 1 else [[c] + NB4[c] for c in range(NCELL)]
            self.free = None
        self.init_skipped = 0

    # ── 재료 고리 ──
    def has_for(self, o, x):
        if not self.local:
            return self.free[x] > 0
        fc = self.fc
        for c in self.ballof[o.pos]:
            if fc[c][x] > 0:
                return True
        return False

    def take_for(self, o, x):
        if not self.local:
            return self.take(x)
        fc = self.fc
        for c in self.ballof[o.pos]:
            f = fc[c]
            if f[x] > 0:
                f[x] -= 1
                return True
        return False

    def give_at(self, c, letters):
        if not self.local:
            self.give(letters)
            return
        f = self.fc[c]
        for x in letters:
            f[x] += 1

    # ── 자리 고리 ──
    def claim_slot(self, o):
        cell = self.cell
        empt = [n for n in NB4[o.pos] if not cell[n]]
        if not empt:
            return False
        n = empt[0] if len(empt) == 1 else empt[int(self.place.next() * len(empt))]
        cell[n] = 1
        o.cpos = n
        self.occ += 1
        return True

    def drop_child(self, o):
        self.give_at(o.cpos, o.buf)
        self.cell[o.cpos] = 0
        o.cpos = -1
        self.occ -= 1

    def place_child(self, o, child):
        child.pos = o.cpos
        o.cpos = -1

    def free_org(self, o):
        self.give_at(o.pos, o.code)
        self.cell[o.pos] = 0
        if o.slot:
            self.give_at(o.cpos, o.buf)
            self.cell[o.cpos] = 0
        elif o.buf:
            self.cons_bad.append("t%d 자리 없는 자식 글자(죽음)" % self.t)
        self.occ -= 1 + (1 if o.slot else 0)

    def free_total(self):
        if not self.local:
            return list(self.free)
        return [sum(self.fc[c][x] for c in range(NCELL)) for x in range(NKIND)]

    def recount(self):
        tot = sum(self.free_total())
        occ = 0
        seen = [0] * NCELL
        for o in self.alive:
            tot += len(o.code) + len(o.buf)
            occ += 1 + (1 if o.slot else 0)
            if (not o.slot) and o.buf:
                self.cons_bad.append("t%d 자리 없는 자식 글자" % self.t)
            if not (0 <= o.pos < NCELL):
                self.cons_bad.append("t%d 칸 없는 생물" % self.t)
                continue
            seen[o.pos] += 1
            if o.slot:
                if o.cpos not in NB4[o.pos]:
                    self.cons_bad.append("t%d 자식 칸이 이웃 아님" % self.t)
                else:
                    seen[o.cpos] += 1
        self.cons_checks += 1
        if tot != self.M:
            self.cons_bad.append("t%d 글자 합 %d != %d" % (self.t, tot, self.M))
        if occ != self.occ or occ > self.K:
            self.cons_bad.append("t%d 점유 %d / 계수 %d / K %d" % (self.t, occ, self.occ, self.K))
        if max(seen) > 1:
            self.cons_bad.append("t%d 한 칸에 둘" % self.t)
        if seen != self.cell:
            self.cons_bad.append("t%d 칸 표와 생물 위치 불일치" % self.t)
        if self.local:
            if min(min(f) for f in self.fc) < 0:
                self.cons_bad.append("t%d 음수 자유 글자" % self.t)
        elif min(self.free) < 0:
            self.cons_bad.append("t%d 음수 자유 글자" % self.t)

    def snapshot(self):
        orgs = []
        for o in self.alive:
            orgs.append([list(o.code), o.L, o.lin, o.phase, o.ptr, list(o.buf), o.slot,
                         (o.nxt - self.t) if o.phase != IDLE else -1, o.err, o.pos, o.cpos])
        s = {"occ": self.occ, "orgs": orgs}
        if self.local:
            s["fc"] = [list(f) for f in self.fc]
        else:
            s["free"] = list(self.free)
        return s


def snap_checksum(snap):
    return hashlib.sha256(json.dumps(snap, sort_keys=True).encode()).hexdigest()[:16]


def initial_cells(M):
    """칸 c 에 ρ 글자 · 종류 (ρc + j) mod 11 — 종류별 총량 = divmod(M, 11) 분할(mini26 과 같다)"""
    rho = M // NCELL
    assert rho * NCELL == M
    fc = [[0] * NKIND for _ in range(NCELL)]
    for c in range(NCELL):
        for j in range(rho):
            fc[c][(rho * c + j) % NKIND] += 1
    return fc


def build(snap, K, M, rho, death, mu, rule, rs, rd, rw):
    """glob — mini26.build 그대로"""
    w = World(K, M, rho, death, mu, rule, rs, rd, rw)
    w.free = list(snap["free"])
    w.occ = snap["occ"]
    for code, L, lin, ph, ptr, buf, slot, off, err in snap["orgs"]:
        o = Org(tuple(code), L, lin)
        o.phase, o.ptr, o.buf, o.slot, o.err = ph, ptr, list(buf), slot, err
        w.add(o, off if ph != IDLE else 0)
    return w


def build_grid(snap, K, M, rho, death, mu, rule, rs, rd, rw, rp, mat_r):
    w = GridWorld(K, M, rho, death, mu, rule, rs, rd, rw, rp, mat_r)
    if w.local:
        w.fc = [list(f) for f in snap["fc"]]
    else:
        w.free = list(snap["free"])
    w.occ = snap["occ"]
    for code, L, lin, ph, ptr, buf, slot, off, err, pos, cpos in snap["orgs"]:
        o = Org(tuple(code), L, lin)
        o.phase, o.ptr, o.buf, o.slot, o.err, o.pos, o.cpos = ph, ptr, list(buf), slot, err, pos, cpos
        w.cell[pos] = 1
        if slot:
            w.cell[cpos] = 1
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
    """glob — mini26.grow 그대로(흐름 이름 · 순서 같음)"""
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
    return _grow_loop(w, M, T, every)


def _grow_loop(w, M, T, every):
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


def grow_grid(rho, seed, K, death, T, n0, every, geom, mat_r):
    M = int(round(rho * K))
    w = GridWorld(K, M, rho, death, 0.0, "roll",
                  stream(seed, rho, "grow", geom, "sched"), stream(seed, rho, "grow", geom, "death"),
                  stream(seed, rho, "grow", geom, "draw"), stream(seed, rho, "grow", geom, "place"), mat_r)
    fc = initial_cells(M)
    if w.local:
        w.fc = fc
    else:
        w.free = [sum(fc[c][x] for c in range(NCELL)) for x in range(NKIND)]
    cells = stream(seed, rho, "grow", geom, "init").permutation(NCELL)[:n0]
    for c in cells:
        o = Org(CODE, L_LONG, 0)
        o.pos = int(c)
        got = []
        for x in CODE:
            if w.take_for(o, x):
                got.append(x)
            else:
                break
        if len(got) < len(CODE):             # 이웃에 글자가 모자라면 되돌리고 건너뜀(센다)
            w.give_at(o.pos, got)
            w.init_skipped += 1
            continue
        w.cell[o.pos] = 1
        w.occ += 1
        w.add(o, 1)
    return _grow_loop(w, M, T, every)


def inject(snap, seed, rho, code_L, frac):
    """glob — mini26.inject 그대로."""
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


def inject_grid(snap, seed, rho, geom, code_L, frac, mat_r):
    """loc · place — 같은 규칙: 바뀐 개체의 몸 글자는 자기 칸에 · 자식 글자는 자식 칸에 돌려주고
    새 코드 글자를 (loc) 자기 칸 + 반경 r 에서 / (place) 전역 풀에서 가져온다. 못 하면 되돌리고 건너뜀."""
    s = json.loads(json.dumps(snap))
    orgs = s["orgs"]
    local = "fc" in s
    bl = (lambda c: [c] + NB4[c]) if mat_r == 1 else (lambda c: ball(c, mat_r))
    k = int(round(frac * len(orgs)))
    pick = stream(seed, rho, geom, "inject").permutation(len(orgs))[:k]
    skipped = 0
    for i in sorted(int(j) for j in pick):
        code, L, lin, ph, ptr, buf, slot, off, err, pos, cpos = orgs[i]
        if local:
            fc = s["fc"]
            for x in code:
                fc[pos][x] += 1
            if slot:
                for x in buf:
                    fc[cpos][x] += 1
            took = []
            for x in CODE:
                for c in bl(pos):
                    if fc[c][x] > 0:
                        fc[c][x] -= 1
                        took.append((c, x))
                        break
            if len(took) == len(CODE):
                ok = True
            else:
                ok = False
                for c, x in took:
                    fc[c][x] += 1
                for x in code:
                    fc[pos][x] -= 1
                if slot:
                    for x in buf:
                        fc[cpos][x] -= 1
        else:
            free = s["free"]
            back = list(code) + list(buf)
            for x in back:
                free[x] += 1
            if all(free[x] >= 1 for x in CODE):
                for x in CODE:
                    free[x] -= 1
                ok = True
            else:
                for x in back:
                    free[x] -= 1
                ok = False
        if ok:
            if slot:
                s["occ"] -= 1
            orgs[i] = [list(CODE), code_L, 1, CLAIM, 0, [], False, 1, False, pos, -1]
        else:
            skipped += 1
    return s, k, skipped


def _assay(w, T, every):
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
    return series, ext


def run_arm(snap, M, seed, rho, mu, rule, kind, code_L, K, death, T, every, frac):
    """glob — mini26.run_arm 그대로(출력 키도 같다)"""
    s, k, skipped = inject(snap, seed, rho, code_L, frac)
    tag = "mu%g" % mu
    w = build(s, K, M, rho, death, mu, rule,
              stream(seed, rho, tag, kind, "sched"), stream(seed, rho, tag, kind, "death"),
              stream(seed, rho, tag, kind, "draw"))
    series, ext = _assay(w, T, every)
    cnt = {"%d_%d" % key: dict(zip(CNAMES, v)) for key, v in sorted(w.cnt.items())}
    return {"rule": rule, "mu": mu, "kind": kind, "L_inj": code_L, "n_inject": k, "skipped": skipped,
            "series": series, "extinct_tick": ext, "counters": cnt,
            "cons_ok": not w.cons_bad, "cons_checks": w.cons_checks, "cons_bad": w.cons_bad[:5],
            "end_free": w.free}


def run_arm_grid(snap, M, seed, rho, geom, mat_r, mu, rule, kind, code_L, K, death, T, every, frac):
    s, k, skipped = inject_grid(snap, seed, rho, geom, code_L, frac, mat_r)
    tag = "mu%g" % mu
    w = build_grid(s, K, M, rho, death, mu, rule,
                   stream(seed, rho, geom, tag, kind, "sched"), stream(seed, rho, geom, tag, kind, "death"),
                   stream(seed, rho, geom, tag, kind, "draw"), stream(seed, rho, geom, tag, kind, "place"), mat_r)
    series, ext = _assay(w, T, every)
    cnt = {"%d_%d" % key: dict(zip(CNAMES, v)) for key, v in sorted(w.cnt.items())}
    return {"rule": rule, "mu": mu, "kind": kind, "L_inj": code_L, "n_inject": k, "skipped": skipped,
            "series": series, "extinct_tick": ext, "counters": cnt,
            "cons_ok": not w.cons_bad, "cons_checks": w.cons_checks, "cons_bad": w.cons_bad[:5],
            "end_free": w.free_total()}


def place_plan(mus):
    """place 기하(서술 대조)는 roll × 확장 사다리 + find × μ 0(규칙 동일성 가드용)만 돈다"""
    return [("roll", mu) for mu in mus if mu in MUS_EXT] + ([("find", 0.0)] if 0.0 in mus else [])


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
    ap.add_argument("--mat-radius", default="1", help="loc 기하 재료 반경(맨해튼) · 'all' = 격자 전체(전역 풀과 같은 뜻)")
    ap.add_argument("--mus", default=",".join("%g" % m for m in MUS_ALL))
    ap.add_argument("--rules", default="roll,find")
    ap.add_argument("--geoms", default=",".join(GEOMS))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mus = [float(x) for x in a.mus.split(",")]
    rules = a.rules.split(",")
    geoms = a.geoms.split(",")
    assert all(r in RULES for r in rules) and all(g in GEOMS for g in geoms)
    assert a.K == NCELL, "K 는 격자 칸 수(%d)와 같아야 한다" % NCELL
    mat_r = None if a.mat_radius == "all" else int(a.mat_radius)
    if mat_r is not None and mat_r >= DIAM:
        mat_r = None                            # 반경이 격자 지름 이상 = 전역 풀
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    out = {"version": VERSION, "config": {k: getattr(a, k) for k in
                                          ("rho", "seed", "K", "death", "grow", "assay", "every", "n0",
                                           "frac", "nref", "min_pop", "mat_radius")},
           "grid": [NX, NY], "mus": mus, "rules": rules, "geoms": geoms, "M": int(round(a.rho * a.K)), "by_geom": {}}
    kinds = [("mut", L_SHORT)] + [("ref%d" % i, L_LONG) for i in range(a.nref)]
    msg = []
    for geom in geoms:
        tg = time.time()
        if geom == "glob":
            w, M, gseries, gext = grow(a.rho, a.seed, a.K, a.death, a.grow, a.n0, a.every)
            init_skipped = 0
        else:
            gr = mat_r if geom == "loc" else None
            w, M, gseries, gext = grow_grid(a.rho, a.seed, a.K, a.death, a.grow, a.n0, a.every, geom, gr)
            init_skipped = w.init_skipped
        snap = w.snapshot()
        gcnt = {"%d_%d" % key: dict(zip(CNAMES, v)) for key, v in sorted(w.cnt.items())}
        n_end = len(w.alive)
        qualified = gext < 0 and n_end >= a.min_pop
        G = {"grow": {"series": gseries, "extinct_tick": gext, "n_end": n_end, "qualified": qualified,
                      "checksum": snap_checksum(snap), "counters": gcnt, "cons_ok": not w.cons_bad,
                      "cons_checks": w.cons_checks, "cons_bad": w.cons_bad[:5], "end_free": w.free_total(),
                      "occ_end": w.occ},
             "arms": []}
        if geom != "glob":
            G["grow"]["init_skipped"] = init_skipped
        if qualified:
            if geom == "place":
                plan = place_plan(mus)
            else:
                plan = [(rule, mu) for rule in rules for mu in mus]
            for rule, mu in plan:
                for kind, L in kinds:
                    if geom == "glob":
                        G["arms"].append(run_arm(snap, M, a.seed, a.rho, mu, rule, kind, L, a.K, a.death,
                                                 a.assay, a.every, a.frac))
                    else:
                        gr = mat_r if geom == "loc" else None
                        G["arms"].append(run_arm_grid(snap, M, a.seed, a.rho, geom, gr, mu, rule, kind, L, a.K,
                                                      a.death, a.assay, a.every, a.frac))
        G["seconds"] = round(time.time() - tg, 1)
        out["by_geom"][geom] = G
        msg.append("%s:q=%s,n=%d,arms=%d,%.0fs" % (geom, qualified, n_end, len(G["arms"]), G["seconds"]))
    out["seconds"] = round(time.time() - t0, 1)
    fn = os.path.join(a.out, "mini27_r%g_s%d.json" % (a.rho, a.seed))
    with open(fn + ".tmp", "w") as f:
        json.dump(out, f)
    os.replace(fn + ".tmp", fn)
    allq = all(out["by_geom"][g]["grow"]["qualified"] for g in geoms)
    print("%s r%g s%d qualified=%s %s %.1fs" % (VERSION, a.rho, a.seed, allq, " ".join(msg), out["seconds"]))


if __name__ == "__main__":
    sys.exit(main())
