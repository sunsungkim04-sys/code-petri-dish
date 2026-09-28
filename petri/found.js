#!/usr/bin/env node
/* 코드 배양접시 — 씨앗 시절 창업 경쟁 계측기(해부 2B). 판정은 하지 않는다(found_analyze.py).
 *
 * 1. 본실험 접시(조건 · 시드)를 같은 설정으로 T0 틱까지 다시 키운다 → 본실험 기록의 T0 개체 수 · 상위 20 과 대조.
 * 2. 팔마다 그 상태를 복제해, 살아 있는 개체 **전부**를 무작위 반반으로 나눠 A 코드(표지 0)와 B 코드(표지 1)로 바꾸고
 *    팔 μ 로 ticks 틱 돌려 표지 1 · 표지 0 계통 개체 수를 every 틱마다 센다. 반반 나누기는 배경 시드로만 정한다(팔끼리 같다).
 *    재료가 모자라 바꾸지 못한 개체는 표지 2 로 남긴다(세지 않음).
 * 3. 팔 = μ 마다 × 반복 r = 1..R × { neutral: A = B = 씨앗 · test: A = 씨앗, B = 짝 코드 } — 같은 r 의 두 팔은 뒤 난수가 같다.
 *    두 계통 중 하나가 사라지거나 접시가 멸종하면 그 팔은 거기서 멈춘다.
 * 4. 계측기 점검(프로세스마다): 작은 접시로 복제 충실도 — dms.js 와 같다.
 *
 * 사용: node found.js --cond mat --seed 3 --t0 2000 --b racld [--a rascld] [--reps 5] [--ticks 3000] [--every 250]
 *                     [--mus 0,0.01] [--main-dir ../main] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
const COND = {
  space: { mat: false, energy: false },
  mat: { mat: true, energy: false },
  energy: { mat: false, energy: true }
};
const cond = args.cond;
if (!COND[cond]) { console.error(`알 수 없는 cond: ${cond}`); process.exit(2); }
if (!['0.3.1', '0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.1~0.3.5 가 필요하다`); process.exit(2); }
const seed = parseInt(args.seed, 10);
const T0 = parseInt(args.t0, 10);
const A = args.a || Sim.ANCESTOR;
const B = args.b;
const R = parseInt(args.reps || '5', 10);
const T = parseInt(args.ticks || '3000', 10);
const every = parseInt(args.every || '250', 10);
const MUS = (args.mus || '0,0.01').split(',').map(Number);
const mainDir = args['main-dir'] || path.join(__dirname, 'main');
const outDir = args.out || '.';
const { LET, K, DISH, NB } = Sim;
const enc = s => [...s].map(ch => { const k = LET.indexOf(ch); if (k < 0) throw new Error(`글자 아님: ${ch}`); return k; });
const opts = Object.assign({ seed, density: 8, mu: 0.01, light: 1 }, COND[cond]);
const INJ_SEED = (seed * 1000003 + 29) >>> 0;
const POST = r => (seed * 7919 + 1 + r * 104729) >>> 0;
const t0 = Date.now();

/* ---------- 세계 복제 (dms.js 와 같다) ---------- */
function cloneWorld(w, rngSeed) {
  if (w.born.length) throw new Error('born 이 비어 있지 않은 시점에 복제하려 했다');
  const orgs = new Map();
  for (const [id, o] of w.orgs) {
    orgs.set(id, Object.assign({}, o, { g: o.g.slice(), child: o.child ? { i: o.child.i, g: o.child.g.slice() } : null }));
  }
  const alive = w.alive.map(o => orgs.get(o.id));
  if (alive.length !== orgs.size || alive.some(o => !o)) throw new Error('alive 와 orgs 가 어긋난다');
  const over = {
    o: Object.assign({}, w.o), rnd: Sim.mulberry32(rngSeed),
    cell: w.cell.slice(), hueAt: w.hueAt.slice(), mat: w.mat ? w.mat.slice() : null,
    orgs, alive, born: [],
    hist: { pop: w.hist.pop.slice(), kinds: w.hist.kinds.slice(), len: w.hist.len.slice() }, summary: null
  };
  for (const [k, v] of Object.entries(w)) {
    if (v !== null && typeof v === 'object' && !(k in over)) throw new Error(`복제에서 빠진 참조 필드: ${k}`);
  }
  return Object.assign({}, w, over);
}

function fullState(w) {
  const orgs = [...w.orgs.values()].sort((a, b) => a.id - b.id)
    .map(o => [o.id, o.i, o.key, o.ip, o.ctx, o.rh, o.face, o.age, o.maxAge, +o.E.toFixed(6), o.bornAt, o.last, o.lastHost, o.tag,
      o.child ? [o.child.i, o.child.g.join('')] : null]);
  return JSON.stringify([w.tick, w.nextId, w.births, w.deaths, w.alive.map(o => o.id), orgs,
    Array.from(w.cell), w.mat ? Array.from(w.mat) : null]);
}

function cloneFidelityTest() {
  const w = Sim.makeWorld(Object.assign({}, opts, { seed: 424242 }));
  for (let t = 0; t < 1500 && w.extinctAt < 0; t++) Sim.tick(w);
  const c = cloneWorld(w, 777);
  w.rnd = Sim.mulberry32(777);
  for (let t = 0; t < 500; t++) { Sim.tick(w); Sim.tick(c); }
  return { pop: w.alive.length, same: fullState(w) === fullState(c), checksum: [Sim.checksum(w), Sim.checksum(c)] };
}

function takeNear(w, i, k, rnd) {
  const m = w.mat;
  if (m[i * K + k] > 0) { m[i * K + k]--; return i; }
  for (let d = 0; d < 4; d++) {
    const t = NB[i * 4 + d];
    if (t >= 0 && m[t * K + k] > 0) { m[t * K + k]--; return t; }
  }
  for (let q = 0; q < 200000; q++) {
    const t = DISH[(rnd() * DISH.length) | 0];
    if (m[t * K + k] > 0) { m[t * K + k]--; return t; }
  }
  return -1;
}

/* ---------- 전부 반반 바꾸기 ---------- */
function convertAll(w, GA, GB) {
  const rnd = Sim.mulberry32(INJ_SEED);
  const list = w.alive.slice().sort((a, b) => a.id - b.id);
  for (let q = list.length - 1; q > 0; q--) {
    const j = (rnd() * (q + 1)) | 0;
    const tmp = list[q]; list[q] = list[j]; list[j] = tmp;
  }
  const half = list.length >> 1;
  const keyA = Sim.keyOf(GA), keyB = Sim.keyOf(GB);
  const res = { n: list.length, to_b: half, to_a: list.length - half, done: 0, skipped: 0 };
  for (let q = 0; q < list.length; q++) {
    const org = list[q], toB = q < half, G = toB ? GB : GA;
    if (org.child) {
      if (w.mat) for (const k of org.child.g) w.mat[org.child.i * K + k]++;
      w.cell[org.child.i] = -1; org.child = null;
    }
    if (w.mat) {
      for (const k of org.g) w.mat[org.i * K + k]++;
      const taken = [];
      let ok = true;
      for (const k of G) { const at = takeNear(w, org.i, k, rnd); if (at < 0) { ok = false; break; } taken.push([at, k]); }
      if (!ok) {
        for (const [at, k] of taken) w.mat[at * K + k]++;
        for (const k of org.g) w.mat[org.i * K + k]--;
        org.tag = 2; res.skipped++;
        continue;
      }
    }
    org.g = G.slice(); org.key = toB ? keyB : keyA; org.hue = Sim.hueOf(org.key); w.hueAt[org.i] = org.hue;
    org.ip = 0; org.ctx = org.id; org.rh = 0; org.tag = toB ? 1 : 0;
    res.done++;
  }
  return res;
}

function runArm(base, kind, bKey, mu, r) {
  const ta = Date.now();
  const w = cloneWorld(base, POST(r));
  w.o.mu = mu;
  const conv = convertAll(w, enc(A), enc(bKey));
  const consAfter = w.mat ? Sim.countMat(w) === w.matTotal0 : true;
  const count = () => {
    let m1 = 0, m0 = 0;
    for (const o of w.alive) { if (o.tag === 1) m1++; else if (o.tag === 0) m0++; }
    return [w.alive.length, m1, m0];
  };
  const series = [[0, ...count()]];
  let stopped = -1;
  for (let t = 1; t <= T; t++) {
    Sim.tick(w);
    if (t % every === 0 || w.extinctAt >= 0) {
      const c = count();
      series.push([t, ...c]);
      if (c[1] === 0 || c[2] === 0 || w.extinctAt >= 0) { stopped = t; break; }
    }
  }
  const consEnd = w.mat ? Sim.countMat(w) === w.matTotal0 : true;
  return { kind, a: A, b: bKey, mu, rep: r, post_seed: POST(r), conv, cons_ok: consAfter && consEnd, stopped,
    extinct: w.extinctAt >= 0, series, elapsed_ms: Date.now() - ta };
}

/* ---------- 실행 ---------- */
const fidelity = cloneFidelityTest();
const w = Sim.makeWorld(opts);
while (w.tick < T0 && w.extinctAt < 0) Sim.tick(w);
const counts = new Map();
for (const o of w.alive) counts.set(o.key, (counts.get(o.key) || 0) + 1);
const mine = [...counts.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1)).slice(0, 20);
let early = { checked: false, match: false };
const mainPath = path.join(mainDir, `${cond}_d8_mu0p01_s${String(seed).padStart(5, '0')}.json`);
if (fs.existsSync(mainPath)) {
  const z = JSON.parse(fs.readFileSync(mainPath, 'utf8'));
  const rec = z.samples.find(r => r[0] === T0);
  const top = z.tops.find(r => r[0] === T0);
  early = { checked: true, pop_main: rec ? rec[1] : null, pop_here: w.alive.length,
    top_match: !!top && JSON.stringify(top[1]) === JSON.stringify(mine) };
  early.match = early.pop_main === early.pop_here && early.top_match;
}
const arms = [];
if (w.extinctAt < 0 && w.alive.length >= 2) {
  for (const mu of MUS) {
    for (let r = 1; r <= R; r++) {
      arms.push(runArm(w, 'neutral', A, mu, r));
      arms.push(runArm(w, 'test', B, mu, r));
    }
  }
}
const out = {
  sim_version: Sim.VERSION, cond, seed, t0: T0, a: A, b: B, reps: R, ticks: T, every, mus: MUS, opts, inj_seed: INJ_SEED,
  fidelity, early, grow_pop: w.alive.length, grow_extinct: w.extinctAt, grow_top: mine.slice(0, 5), elapsed_ms: Date.now() - t0, arms
};
fs.mkdirSync(outDir, { recursive: true });
const tag = `found_${cond}_${B}_t${T0}_s${String(seed).padStart(5, '0')}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
console.log([
  `DONE ${tag}`, `v${Sim.VERSION}`, `fidelity ${fidelity.same ? 'OK' : 'FAIL'}`,
  `early ${early.checked ? (early.match ? 'MATCH' : 'MISMATCH') : 'no-main'}`, `pop ${w.alive.length}`, `arms ${arms.length}`,
  arms.length ? `conv ${arms[0].conv.done}/${arms[0].conv.n} skip ${arms[0].conv.skipped}` : 'no-arms',
  `cons ${arms.every(a => a.cons_ok) ? 'OK' : 'FAIL'}`, `total ${((Date.now() - t0) / 1000).toFixed(0)}s`
].join(' · '));
