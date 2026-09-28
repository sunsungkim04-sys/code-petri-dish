#!/usr/bin/env node
/* 코드 배양접시 — 침입 계통의 나이별 인구학 계측기(해부 8 B). 판정은 하지 않는다(lifehist_analyze.py).
 *
 * 해부 7 ③b 는 μ 0 침입에서 씨앗 계통의 이득을 출생 · 사망 · 나머지 공동체로 갈랐다. 여기서는 그 이득이
 *   (1) 3,000틱 동안 꾸준히 쌓이나(복리) 끼운 직후에 생기나 — 팔을 6,000틱까지 늘려 Δ(t) 를 본다
 *   (2) 나이별 출생률 · 사망률이 달라서인가, 나이 구성(더 젊은 계통)이 달라서인가 — 나이 칸별로 센다
 *   를 잰다.
 * 배경 · 끼우기 · 난수는 dms.js 와 같다: 본실험 접시를 50,000틱까지 다시 키워(체크섬 대조) 살아 있는 개체의 10% 에
 *   코드를 끼우고 μ 0 으로 돌린다. 팔 = ref(WT, 뒤 난수 A) + neutral 10(WT, 뒤 난수 B1..B10) + --codes(뒤 난수 A).
 *   → --ticks 3000 이면 팔마다 계통 개체 수 기록(series)이 dms.js(해부 1)의 같은 팔과 같아야 한다(회귀 검사).
 * 부모를 알려고 틱마다 개체의 tag 를 자기 id 로 덮어쓴다(lineage.js 와 같은 방법 — tag 는 동작 · 체크섬에 안 쓰인다).
 *   계통(1 = 끼운 계통 · 0 = 나머지)은 따로 들고 다닌다: 끼운 개체 1, 자식은 부모 것.
 * 칸(bin, 기본 250틱)마다 계통 × 나이 칸(50틱, 0~11)별로:
 *   ot   노출 — 틱 시작에 살아 있던 개체-틱(틱 시작 나이 칸)
 *   b    출생 — 부모의 틱 시작 나이 칸(그 틱에 태어나 살아남은 자식)
 *   da   나이 사망(수명 넘음) · do 그 밖 사망(차례 전 · 무작위 · 먹힘) — 틱 시작 나이 칸
 *   칸마다 계통별 시작 개체 수 n0 · 먹기(e) 글자를 가진 개체 수 ne 와 놓친 출생(태어난 틱에 먹힌 자식 — 계통 모름) 수도 적는다.
 *   (ne 는 마른 실행 뒤 더함 — 파일럿에서 '그 밖 사망'이 씨앗 팔 전체에서 낮았다)
 *
 * 사용: node lifehist.js --seed 24 [--codes rascld] [--ticks 6000] [--bin 250] [--main-dir main] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');
const LIFEHIST_VERSION = '1.0.0';

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다`); process.exit(2); }
const seed = parseInt(args.seed, 10);
const WT = 'racld';
const codes = (args.codes || 'rascld').split(',').map(x => x.trim()).filter(Boolean);
const T = parseInt(args.ticks || '6000', 10);
const BIN = parseInt(args.bin || '250', 10);
const every = 250;
const frac = 0.1;
const grow = 50000;
const N_NEUTRAL = 10;
const NA = 12, AW = 50;
const mainDir = args['main-dir'] || path.join(__dirname, 'main');
const outDir = args.out || '.';
const { LET, K, DISH, NB } = Sim;
const enc = s => [...s].map(ch => { const k = LET.indexOf(ch); if (k < 0) throw new Error(`글자 아님: ${ch}`); return k; });
const opts = { seed, density: 8, mu: 0.01, light: 1, mat: true, energy: false };
const INJ_SEED = (seed * 1000003 + 17) >>> 0;
const POST_A = (seed * 7919 + 1) >>> 0;
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

/* ---------- 끼워 넣기 (dms.js 와 같다) ---------- */
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

function inject(w, G) {
  const rnd = Sim.mulberry32(INJ_SEED);
  const list = w.alive.slice().sort((a, b) => a.id - b.id);
  const n = list.length, target = Math.round(n * frac);
  for (let q = 0; q < target; q++) {
    const j = q + ((rnd() * (n - q)) | 0);
    const tmp = list[q]; list[q] = list[j]; list[j] = tmp;
  }
  const key = Sim.keyOf(G), hue = Sim.hueOf(key);
  let done = 0, skipped = 0;
  for (let q = 0; q < target; q++) {
    const org = list[q];
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
        skipped++;
        continue;
      }
    }
    org.g = G.slice(); org.key = key; org.hue = hue; w.hueAt[org.i] = hue;
    org.ip = 0; org.ctx = org.id; org.rh = 0; org.tag = 1;
    done++;
  }
  return { n, target, done, skipped };
}

/* ---------- 팔 하나 ---------- */
const zeros = () => [new Array(NA).fill(0), new Array(NA).fill(0)];
function newBin(t, n0, ne) { return { t, n0, ne, ot: zeros(), b: zeros(), da: zeros(), do: zeros(), missed: 0 }; }
const ab = a => { const k = (a / AW) | 0; return k < NA ? k : NA - 1; };

function runArm(base, kind, key, labels, postSeed) {
  const ta = Date.now();
  const w = cloneWorld(base, postSeed);
  w.o.mu = 0;
  const inj = inject(w, enc(key));
  const consAfterInj = Sim.countMat(w) === w.matTotal0;
  let cap = Math.max(1 << 16, w.nextId * 2);
  let lin = new Int8Array(cap), age0 = new Int32Array(cap);
  const ensure = id => {
    if (id < cap) return;
    while (cap <= id) cap *= 2;
    const a = new Int8Array(cap); a.set(lin); lin = a;
    const b = new Int32Array(cap); b.set(age0); age0 = b;
  };
  for (const o of w.alive) { ensure(o.id); lin[o.id] = o.tag === 1 ? 1 : 0; }
  const count = () => { let m = 0; for (const o of w.alive) if (lin[o.id] === 1) m++; return [w.alive.length, m]; };
  const series = [[0, ...count()]];
  const bins = [];
  let stopped = -1, missedAll = 0, lostParent = 0;
  for (let t = 1; t <= T; t++) {
    if ((t - 1) % BIN === 0) {
      const c = count(), ne = [0, 0];
      for (const o of w.alive) if (o.key.includes('e')) ne[lin[o.id]]++;
      bins.push(newBin(t - 1, [c[0] - c[1], c[1]], ne));
    }
    const B = bins[bins.length - 1];
    const pre = w.alive.slice();
    for (const o of pre) { ensure(o.id); o.tag = o.id; age0[o.id] = o.age; B.ot[lin[o.id]][ab(o.age)]++; }
    const b0 = w.births;
    Sim.tick(w);
    const bt = w.tick - 1;
    let here = 0;
    for (const o of w.alive) {
      if (o.bornAt !== bt) continue;
      here++;
      ensure(o.id);
      const p = o.tag;
      if (p < 0 || p >= cap) { lostParent++; continue; }
      const L = lin[p];
      lin[o.id] = L;
      B.b[L][ab(age0[p])]++;
    }
    const missed = (w.births - b0) - here;
    B.missed += missed; missedAll += missed;
    for (const o of pre) {
      if (!o.dead) continue;
      const L = lin[o.id], a = ab(age0[o.id]);
      const stepped = o.age === age0[o.id] + 1;
      if (stepped && o.age > o.maxAge) B.da[L][a]++; else B.do[L][a]++;
    }
    if (t % every === 0 || w.extinctAt >= 0) {
      const c = count();
      series.push([t, ...c]);
      if (c[1] === 0 || w.extinctAt >= 0) { stopped = t; break; }
    }
  }
  const consEnd = Sim.countMat(w) === w.matTotal0;
  return { kind, key, labels, post_seed: postSeed, inj, cons_ok: consAfterInj && consEnd, stopped, series, bins,
    missed_births: missedAll, lost_parent: lostParent, elapsed_ms: Date.now() - ta };
}

/* ---------- 실행 ---------- */
const w = Sim.makeWorld(opts);
while (w.tick < grow && w.extinctAt < 0) Sim.tick(w);
const growChecksum = Sim.checksum(w);
const mainPath = path.join(mainDir, `mat_d8_mu0p01_s${String(seed).padStart(5, '0')}.json`);
const main = fs.existsSync(mainPath) ? JSON.parse(fs.readFileSync(mainPath, 'utf8')) : null;
const counts = new Map();
for (const o of w.alive) counts.set(o.key, (counts.get(o.key) || 0) + 1);
const top = [...counts.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1));
const growMs = Date.now() - t0;
const arms = [];
if (w.extinctAt < 0) {
  arms.push(runArm(w, 'ref', WT, ['wt'], POST_A));
  for (let b = 1; b <= N_NEUTRAL; b++) arms.push(runArm(w, 'neutral', WT, [`neutral:${b}`], (POST_A + b * 104729) >>> 0));
  for (const c of codes) arms.push(runArm(w, 'mut', c, [c], POST_A));
}
const out = {
  lifehist_version: LIFEHIST_VERSION, sim_version: Sim.VERSION, cond: 'mat', seed, wt: WT, codes, frac, ticks: T, bin: BIN, every,
  age_width: AW, n_age: NA, main_dir: mainDir, inj_seed: INJ_SEED, post_a: POST_A,
  grow_checksum: growChecksum, main_checksum: main ? main.checksum : null, checksum_match: !!main && main.checksum === growChecksum,
  grow_top: top.slice(0, 5), grow_pop: w.alive.length, grow_extinct: w.extinctAt, grow_ms: growMs,
  arms, elapsed_ms: Date.now() - t0
};
fs.mkdirSync(outDir, { recursive: true });
const tag = `lh_mat_s${String(seed).padStart(5, '0')}_t${T}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
console.log([
  `DONE ${tag}`, `v${Sim.VERSION}/lh${LIFEHIST_VERSION}`,
  `checksum ${main ? (out.checksum_match ? 'MATCH' : 'MISMATCH') : 'no-main'}`,
  `top ${top.length ? top[0][0] : '-'}`, `pop ${w.alive.length}`, `arms ${arms.length}`,
  `missed ${arms.reduce((s, a) => s + a.missed_births, 0)}`, `lost ${arms.reduce((s, a) => s + a.lost_parent, 0)}`,
  `cons ${arms.every(a => a.cons_ok) ? 'OK' : 'FAIL'}`,
  `grow ${(growMs / 1000).toFixed(0)}s`, `total ${((Date.now() - t0) / 1000).toFixed(0)}s`
].join(' · '));
