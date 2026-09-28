/* 코드 배양접시 — 시뮬레이션 핵심 (브라우저 · Node 공용, DOM 없음)
 *
 * 세계: N×N 격자 안의 둥근 접시. 칸마다 개체(= 짧은 코드) 하나.
 * 판정자 없음 — 살아남아 복제한 것이 곧 작동한 것.
 * 모자라는 것: 자리(항상) · 재료(opts.mat, 글자 총량 보존) · 에너지(opts.energy, 빛 두 얼룩).
 * 같은 opts(시드 포함)면 틱마다 똑같이 흘러간다 — 난수는 전부 w.rnd 하나에서 나온다.
 *
 * v0.3.0 (2026-09-15) — 11번째 글자 x(무늬) 추가: 실행해도 아무 일 없고 틱 · 에너지를 쓰지 않는다.
 *   코드 길이를 늘리는 것 말고는 효과가 없는 중립 표지 — 사건 문턱의 조건별 중립 대조용.
 *   v0.2.0 과는 돌연변이 글자 풀(10 → 11)과 재료 글자 종류가 달라 결과를 섞어 쓰면 안 된다.
 * v0.3.1 (2026-09-15) — 개체에 계통 표지 tag 추가(자식이 부모 것을 물려받음, 씨앗 0) · mulberry32 내보냄.
 *   난수 소비와 규칙은 그대로라 v0.3.0 과 같은 시드는 같은 궤적 — 본실험 접시 체크섬으로 대조한다(해부 실험 dms.js).
 * v0.3.2 (2026-09-16) — opts.ancestor 로 접시의 첫 코드를 고를 수 있게 함(기본은 ANCESTOR 'rascld').
 *   기본값이면 난수 소비 · 규칙이 v0.3.1 과 같다 — 해부 파일럿 체크섬으로 대조한다.
 * v0.3.3 (2026-09-16) — opts.findFirst(기본 false) 추가: **재료 먼저 세계**(해부 7 ②).
 *   원래 규칙은 복사할 때 오류 주사위를 먼저 굴리고 재료를 찾는다 — 헛손질할 때마다 주사위를 다시 굴린다.
 *   findFirst 면 제 글자가 가까이(제 칸 · 이웃 네 칸) 없을 때 **주사위를 굴리지 않고** 헛손질한다.
 *   재료가 없는 조건(mat 끔)에서는 아무것도 바뀌지 않는다. 기본값이면 v0.3.2 와 궤적이 같다 — 저장된 기록의 체크섬으로 대조한다.
 * v0.3.4 (2026-09-17) — opts.rememberDie(기본 false) 추가: **기억하는 주사위 세계**(해부 9 A).
 *   원래 규칙대로 주사위를 먼저 굴리되, 재료가 없어 헛손질하면 그 결과(맞게 · 끼어듦 k · 바뀜 k)를 기억해 두고
 *   다음 시도에서 다시 굴리지 않는다 — 글자 하나에 주사위는 한 번이다. 빠뜨림은 재료가 필요 없어 그 자리에서 끝난다.
 *   기억한 오류(없는 글자로의 끼어듦 · 바뀜)는 그 글자가 올 때까지 기다렸다 실현되므로 '거름'(원래 규칙에서 다시 굴리며 걸러지던 것)도
 *   거의 사라진다 — 파일럿(09-17 · 시드 9401 · μ 1%)에서 R 1.1 · F 0.94 · 부풀림 0.98. 즉 주사위 한 번 + 거름 없음 = 명목에 가장 가까운 세계.
 *   기억은 자식 자리를 새로 잡거나 나눌 때 지운다. 기본값이면 v0.3.3 과 궤적이 같다 — 체크섬으로 대조한다.
 * v0.3.5 (2026-09-18) — 규칙 의존 사다리용 옵션 셋(해부 10): opts.cosmic(기본 true — false 면 우주선 돌연변이를 끄고 그 주사위도 안 굴린다) ·
 *   opts.diffuse(기본 0.15 — 틱마다 재료를 옮기는 시도 수 = 접시 칸 × 이 값) · opts.age0 / opts.ageVar(기본 300 / 300 — 수명 = age0 + [0, ageVar)).
 *   기본값이면 v0.3.4 와 난수 소비 · 궤적이 같다 — 체크섬으로 대조한다.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.PetriSim = factory();
})(typeof self !== 'undefined' ? self : this, function () {
'use strict';

const VERSION = '0.3.5';
const LET = 'nsracldjehx';
const N = 120, C = N * N, K = LET.length;
const X = LET.indexOf('x');
const NAMES = ['쉬기', '표지', '돌기', '자리잡기', '베끼기', '되돌기', '나누기', '빌려쓰기', '먹기', '빛모으기', '무늬'];
const DESCS = [
  '아무것도 안 한다. 틱 하나를 쓰고, 에너지가 가장 적게 든다.',
  '되돌기가 돌아올 자리를 표시한다.',
  '바라보는 방향을 시계방향으로 돌린다.',
  '바라보는 칸이 비었으면 자식 자리를 잡는다.',
  '내 코드 한 글자를 자식에게 베껴 쓴다. 틀릴 수 있다.',
  '자식을 덜 베꼈으면 앞의 표지로 돌아간다.',
  '베낀 자식을 떼어 내 독립시킨다.',
  '바라보는 이웃의 코드로 건너가 그것을 대신 실행한다.',
  '바라보는 이웃을 부순다. 그 재료와 에너지를 챙긴다.',
  '에너지를 켰을 때 이번 틱 빛을 두 배 더 모은다.',
  '아무 일도 안 하고 틱도 쓰지 않는다. 코드 길이만 늘리는 중립 표지.'
];
const COST = [0.3, 0.6, 0.8, 1.0, 1.5, 0.8, 2.0, 1.0, 3.0, 0.5, 0];
const ANCESTOR = 'rascld';
const MINLEN = 2, MAXLEN = 48, AGE0 = 300, AGEVAR = 300, RDEATH = 0.0003, EMAX = 60;
const DEFAULTS = { seed: 1, mat: false, energy: false, density: 4, mu: 0.01, light: 1, ancestor: null, findFirst: false, rememberDie: false,
  cosmic: true, diffuse: 0.15, age0: AGE0, ageVar: AGEVAR };

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function hueOf(key) {
  let h = 2166136261;
  for (let i = 0; i < key.length; i++) { h ^= key.charCodeAt(i); h = Math.imul(h, 16777619); }
  return (h >>> 0) % 360;
}
function keyOf(g) { let s = ''; for (let i = 0; i < g.length; i++) s += LET[g[i]]; return s; }

/* ---------- 접시 기하 · 빛 지도 (고정) ---------- */
const CX = (N - 1) / 2, R = N / 2 - 1.5;
const dishList = [];
const IS_DISH = new Uint8Array(C);
for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
  if ((x - CX) ** 2 + (y - CX) ** 2 <= R * R) { IS_DISH[y * N + x] = 1; dishList.push(y * N + x); }
}
const DISH = Int32Array.from(dishList);
const DIRS = [[0, -1], [1, 0], [0, 1], [-1, 0]];
const NB = new Int32Array(C * 4).fill(-1);
for (const i of DISH) {
  const x = i % N, y = (i / N) | 0;
  for (let d = 0; d < 4; d++) {
    const nx = x + DIRS[d][0], ny = y + DIRS[d][1];
    if (nx >= 0 && ny >= 0 && nx < N && ny < N && IS_DISH[ny * N + nx]) NB[i * 4 + d] = ny * N + nx;
  }
}
const LIGHT = new Float32Array(C);
const L1 = [0.3 * N, 0.35 * N], L2 = [0.7 * N, 0.65 * N], SIG2 = (0.18 * N) ** 2;
let LMAX = 0;
for (let i = 0; i < C; i++) {
  const x = i % N, y = (i / N) | 0;
  const v = 0.25 + 1.6 * (Math.exp(-((x - L1[0]) ** 2 + (y - L1[1]) ** 2) / SIG2) + Math.exp(-((x - L2[0]) ** 2 + (y - L2[1]) ** 2) / SIG2));
  LIGHT[i] = v; if (v > LMAX) LMAX = v;
}

/* ---------- 세계 ---------- */
function makeWorld(opts) {
  const o = Object.assign({}, DEFAULTS, opts);
  const w = {
    o, rnd: mulberry32(o.seed),
    cell: new Int32Array(C), hueAt: new Int16Array(C),
    mat: o.mat ? new Uint16Array(C * K) : null,
    orgs: new Map(), alive: [], born: [], nextId: 0,
    tick: 0, births: 0, deaths: 0, matTotal0: 0, extinctAt: -1,
    hist: { pop: [], kinds: [], len: [] }, summary: null
  };
  w.cell.fill(-2);
  for (const i of DISH) w.cell[i] = -1;
  if (w.mat) {
    const total = Math.round(DISH.length * o.density);
    for (let t = 0; t < total; t++) {
      const i = DISH[(w.rnd() * DISH.length) | 0];
      w.mat[i * K + ((w.rnd() * K) | 0)]++;
    }
  }
  const start = o.energy ? Math.round(L1[1]) * N + Math.round(L1[0]) : (N >> 1) * N + (N >> 1);
  const g = [...(o.ancestor || ANCESTOR)].map(ch => {
    const k = LET.indexOf(ch);
    if (k < 0) throw new Error(`ancestor 에 없는 글자: ${ch}`);
    return k;
  });
  if (w.mat) for (const k of g) takeAnywhere(w, k);
  spawn(w, start, g, null);
  w.alive = w.born; w.born = []; w.births = 0;
  if (w.mat) w.matTotal0 = countMat(w);
  sample(w);
  return w;
}

function spawn(w, i, g, parent) {
  const id = w.nextId++;
  const key = keyOf(g);
  const org = {
    id, i, g, key, hue: hueOf(key), ip: 0, ctx: id, rh: 0, child: null,
    face: (w.rnd() * 4) | 0, age: 0, maxAge: w.o.age0 + ((w.rnd() * w.o.ageVar) | 0),
    E: 40, dead: false, bornAt: w.tick, last: -1, lastHost: id, tag: parent ? parent.tag : 0,
    dk: -1, dtype: 0, dpos: -1   // v0.3.4 기억한 주사위(rememberDie 에서만 씀): 글자 k · 종류(0 맞게 · 1 끼어듦 · 2 바뀜) · 자리
  };
  if (parent && w.o.energy) { org.E = parent.E * 0.5; parent.E *= 0.5; }
  w.cell[i] = id; w.hueAt[i] = org.hue;
  w.orgs.set(id, org); w.born.push(org); w.births++;
  return org;
}

function has(w, i, k) {
  const m = w.mat;
  if (m[i * K + k] > 0) return true;
  for (let d = 0; d < 4; d++) {
    const t = NB[i * 4 + d];
    if (t >= 0 && m[t * K + k] > 0) return true;
  }
  return false;
}
function take(w, i, k) {
  const m = w.mat;
  if (m[i * K + k] > 0) { m[i * K + k]--; return true; }
  for (let d = 0; d < 4; d++) {
    const t = NB[i * 4 + d];
    if (t >= 0 && m[t * K + k] > 0) { m[t * K + k]--; return true; }
  }
  return false;
}
function takeAnywhere(w, k) {
  for (let t = 0; t < 200000; t++) {
    const i = DISH[(w.rnd() * DISH.length) | 0];
    if (w.mat[i * K + k] > 0) { w.mat[i * K + k]--; return true; }
  }
  return false;
}
function release(w, i, g) { for (let j = 0; j < g.length; j++) w.mat[i * K + g[j]]++; }
function countMat(w) {
  let s = 0; const m = w.mat;
  for (let j = 0; j < m.length; j++) s += m[j];
  for (const org of w.alive) { s += org.g.length; if (org.child) s += org.child.g.length; }
  for (const org of w.born) { s += org.g.length; if (org.child) s += org.child.g.length; }
  return s;
}

function abortChild(w, org) {
  const c = org.child;
  if (w.mat) release(w, c.i, c.g);
  w.cell[c.i] = -1; org.child = null;
}
function kill(w, org, dest) {
  if (org.dead) return;
  org.dead = true; w.deaths++;
  w.cell[org.i] = -1; w.orgs.delete(org.id);
  if (w.mat) release(w, dest >= 0 ? dest : org.i, org.g);
  if (org.child) abortChild(w, org);
}

function copyOne(w, org) {
  const c = org.child;
  if (!c || org.rh >= org.g.length) return;
  if (w.o.findFirst && w.mat && !has(w, org.i, org.g[org.rh])) return;   // v0.3.3 재료 먼저: 주사위도 안 굴린다
  let k, adv = true;
  if (w.o.rememberDie && org.dpos >= 0 && org.dpos === org.rh) {
    // v0.3.4 기억한 주사위: 이 자리에서 이미 굴렸으면 그 결과를 다시 쓴다(주사위 · 글자 뽑기 모두 안 굴린다)
    if (org.dtype === 0) k = org.g[org.rh];
    else { k = org.dk; adv = org.dtype !== 1; }
  } else {
    const mu = w.o.mu, r = w.rnd();
    if (r < mu / 3) { org.rh++; org.dpos = -1; return; }       // 빠뜨리기 — 재료가 필요 없어 그 자리에서 끝난다
    k = org.g[org.rh]; let dtype = 0;
    if (r < mu * 2 / 3) { k = (w.rnd() * K) | 0; adv = false; dtype = 1; }  // 끼어들기
    else if (r < mu * 5 / 3) { k = (w.rnd() * K) | 0; dtype = 2; }          // 바뀜
    if (w.o.rememberDie) { org.dk = k; org.dtype = dtype; org.dpos = org.rh; }
  }
  if (w.mat && !take(w, org.i, k)) return;                      // 재료가 없으면 이번 틱은 헛손질 (기억은 남는다)
  org.dpos = -1;                                                // 성공 — 기억을 지운다
  c.g.push(k);
  if (adv) org.rh++;
  if (c.g.length > MAXLEN) abortChild(w, org);
}

function cosmic(w, org) {
  const p = (w.rnd() * org.g.length) | 0, k = (w.rnd() * K) | 0, old = org.g[p];
  if (k === old) return;
  if (w.mat) { if (!take(w, org.i, k)) return; w.mat[org.i * K + old]++; }
  org.g[p] = k; org.key = keyOf(org.g); org.hue = hueOf(org.key); w.hueAt[org.i] = org.hue;
}

function step(w, org) {
  const o = w.o;
  org.age++;
  if (o.energy) org.E += LIGHT[org.i] * o.light;
  let host = org;
  if (org.ctx !== org.id) {
    const h = w.orgs.get(org.ctx);
    if (h && !h.dead) host = h; else { org.ctx = org.id; org.ip = 0; }
  }
  if (org.ip >= host.g.length) { if (host !== org) { host = org; org.ctx = org.id; } org.ip = 0; }
  let g = host.g, ins = g[org.ip];
  // 무늬(x)는 틱을 쓰지 않고 건너뛴다. 코드가 전부 무늬면 그 틱은 아무 일도 없다.
  for (let guard = 0; ins === X && guard < g.length; guard++) {
    org.ip++;
    if (org.ip >= g.length) { if (host !== org) { host = org; org.ctx = org.id; g = org.g; } org.ip = 0; }
    ins = g[org.ip];
  }
  org.last = org.ip; org.lastHost = host.id; org.ip++;
  switch (ins) {
    case 2: org.face = (org.face + 1) & 3; break;
    case 3:
      if (!org.child) {
        const t = NB[org.i * 4 + org.face];
        if (t >= 0 && w.cell[t] === -1) { w.cell[t] = -3; w.hueAt[t] = org.hue; org.child = { i: t, g: [] }; org.rh = 0; org.dpos = -1; }
      }
      break;
    case 4: copyOne(w, org); break;
    case 5:
      if (org.child && org.rh < org.g.length) {
        let p = org.ip - 2;
        while (p >= 0 && g[p] !== 1) p--;
        org.ip = p >= 0 ? p + 1 : 0;
      }
      break;
    case 6:
      if (org.child) {
        if (org.child.g.length >= MINLEN) { const c = org.child; org.child = null; spawn(w, c.i, c.g, org); }
        else abortChild(w, org);
        org.rh = 0; org.dpos = -1;
      }
      break;
    case 7: {
      const t = NB[org.i * 4 + org.face];
      if (t >= 0) {
        const id = w.cell[t];
        if (id >= 0 && id !== org.id) {
          const h = w.orgs.get(id);
          if (h && !h.dead) { org.ctx = id; const p = h.g.indexOf(1); org.ip = p >= 0 ? p : 0; }
        }
      }
      break;
    }
    case 8: {
      const t = NB[org.i * 4 + org.face];
      if (t >= 0) {
        const id = w.cell[t];
        if (id >= 0) {
          const v = w.orgs.get(id);
          if (v && !v.dead) { if (o.energy) org.E += Math.max(0, v.E) * 0.5; kill(w, v, org.i); }
        }
      }
      break;
    }
    case 9: if (o.energy) org.E += LIGHT[org.i] * o.light * 2; break;
    default: break;
  }
  if (o.energy) {
    org.E -= COST[ins];
    if (org.E > EMAX) org.E = EMAX;
    if (org.E <= 0) { kill(w, org, -1); return; }
  }
  if (org.age > org.maxAge || w.rnd() < RDEATH) { kill(w, org, -1); return; }
  if (o.cosmic && o.mu > 0 && w.rnd() < o.mu * 0.02) cosmic(w, org);   // v0.3.5: cosmic 끄면 주사위도 안 굴린다
}

function diffuse(w) {
  const m = w.mat, n = (DISH.length * w.o.diffuse) | 0;
  for (let q = 0; q < n; q++) {
    const i = DISH[(w.rnd() * DISH.length) | 0], k = (w.rnd() * K) | 0;
    if (m[i * K + k] > 0) {
      const t = NB[i * 4 + ((w.rnd() * 4) | 0)];
      if (t >= 0) { m[i * K + k]--; m[t * K + k]++; }
    }
  }
}

function tick(w) {
  const list = w.alive;
  for (let i = list.length - 1; i > 0; i--) {
    const j = (w.rnd() * (i + 1)) | 0; const t = list[i]; list[i] = list[j]; list[j] = t;
  }
  for (let q = 0; q < list.length; q++) { const org = list[q]; if (!org.dead) step(w, org); }
  const next = [];
  for (const org of list) if (!org.dead) next.push(org);
  for (const org of w.born) if (!org.dead) next.push(org);
  w.alive = next; w.born = [];
  if (w.mat) diffuse(w);
  w.tick++;
  if (w.tick % 10 === 0) sample(w);
  if (w.alive.length === 0 && w.extinctAt < 0) { w.extinctAt = w.tick; summarize(w); }
}

function summarize(w) {
  const counts = new Map();
  let len = 0, borrow = 0, e = 0;
  for (const org of w.alive) {
    counts.set(org.key, (counts.get(org.key) || 0) + 1);
    len += org.g.length; if (org.ctx !== org.id) borrow++; e += org.E;
  }
  const n = w.alive.length;
  w.summary = { tick: w.tick, n, kinds: counts.size, len: n ? len / n : 0, borrow, energy: n ? e / n : 0, counts };
  return w.summary;
}
function pushCap(a, v) { a.push(v); if (a.length > 360) a.shift(); }
function sample(w) {
  const s = summarize(w);
  pushCap(w.hist.pop, s.n); pushCap(w.hist.kinds, s.kinds); pushCap(w.hist.len, s.len);
}

/* 결정론 확인용 상태 지문 — 칸 · 코드 · 나이 · 에너지(0.001 단위) · 배지 · 틱 */
function checksum(w) {
  let h = 2166136261 >>> 0;
  const mix = v => {
    for (let s = 0; s < 32; s += 8) { h ^= (v >>> s) & 0xff; h = Math.imul(h, 16777619); }
  };
  for (let i = 0; i < C; i++) {
    const id = w.cell[i];
    mix(id);
    if (id >= 0) {
      const org = w.orgs.get(id);
      for (const k of org.g) mix(k);
      mix(org.age); mix(Math.round(org.E * 1000));
    }
  }
  if (w.mat) for (let j = 0; j < w.mat.length; j++) mix(w.mat[j]);
  mix(w.tick);
  return (h >>> 0).toString(16).padStart(8, '0');
}

return {
  VERSION, N, C, K, X, LET, NAMES, DESCS, COST, ANCESTOR, EMAX, MINLEN, MAXLEN, DEFAULTS,
  DISH, NB, LIGHT, LMAX, R,
  makeWorld, tick, summarize, countMat, checksum, keyOf, hueOf, mulberry32
};
});
