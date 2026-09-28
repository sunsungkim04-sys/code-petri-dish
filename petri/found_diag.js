#!/usr/bin/env node
/* 해부 2B 사후 진단 — 판정 아님. 창업 경쟁의 팔을 그대로 다시 돌리며 계통(표지)마다 무엇이 되어 가는지 기록한다.
 *
 * found.js 와 같은 성장 · 같은 반반 바꾸기 · 같은 뒤 난수 · 같은 멈춤 규칙이다. 기록만 더하고 난수는 쓰지 않으므로
 * 250틱 계통 수가 found 기록의 같은 팔과 **같아야 한다** — 대조 결과를 함께 출력한다(다르면 이 진단을 믿지 말 것).
 * 500틱마다 표지 0(씨앗 계통) · 표지 1(짝 계통) 각각: 개체 수 · `e` 든 코드 비율 · `h` 든 코드 비율 · 출발 코드 그대로인 비율 · 코드 종류 수 · 상위 3
 *
 * 사용: node found_diag.js --cond space --t0 500 --b reascld --seeds 1,3,4 --mu 0.01 [--rep 1] [--kind test] [--found-dir found] [--out found_diag]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
const COND = { space: { mat: false, energy: false }, mat: { mat: true, energy: false }, energy: { mat: false, energy: true } };
const cond = args.cond, T0 = parseInt(args.t0, 10), B = args.b, A = Sim.ANCESTOR;
const SEEDS = args.seeds.split(',').map(Number);
const MU = parseFloat(args.mu), REP = parseInt(args.rep || '1', 10), KIND = args.kind || 'test';
const T = 3000, every = 250, snapEvery = 500;
const foundDir = args['found-dir'] || path.join(__dirname, 'found');
const outDir = args.out || 'found_diag';
const { LET, K, DISH, NB } = Sim;
const enc = s => [...s].map(ch => LET.indexOf(ch));

function cloneWorld(w, rngSeed) {
  const orgs = new Map();
  for (const [id, o] of w.orgs) orgs.set(id, Object.assign({}, o, { g: o.g.slice(), child: o.child ? { i: o.child.i, g: o.child.g.slice() } : null }));
  const alive = w.alive.map(o => orgs.get(o.id));
  return Object.assign({}, w, {
    o: Object.assign({}, w.o), rnd: Sim.mulberry32(rngSeed), cell: w.cell.slice(), hueAt: w.hueAt.slice(),
    mat: w.mat ? w.mat.slice() : null, orgs, alive, born: [],
    hist: { pop: w.hist.pop.slice(), kinds: w.hist.kinds.slice(), len: w.hist.len.slice() }, summary: null
  });
}
function takeNear(w, i, k, rnd) {
  const m = w.mat;
  if (m[i * K + k] > 0) { m[i * K + k]--; return i; }
  for (let d = 0; d < 4; d++) { const t = NB[i * 4 + d]; if (t >= 0 && m[t * K + k] > 0) { m[t * K + k]--; return t; } }
  for (let q = 0; q < 200000; q++) { const t = DISH[(rnd() * DISH.length) | 0]; if (m[t * K + k] > 0) { m[t * K + k]--; return t; } }
  return -1;
}
function convertAll(w, GA, GB, injSeed) {
  const rnd = Sim.mulberry32(injSeed);
  const list = w.alive.slice().sort((a, b) => a.id - b.id);
  for (let q = list.length - 1; q > 0; q--) { const j = (rnd() * (q + 1)) | 0; const tmp = list[q]; list[q] = list[j]; list[j] = tmp; }
  const half = list.length >> 1;
  for (let q = 0; q < list.length; q++) {
    const org = list[q], toB = q < half, G = toB ? GB : GA;
    if (org.child) { if (w.mat) for (const k of org.child.g) w.mat[org.child.i * K + k]++; w.cell[org.child.i] = -1; org.child = null; }
    if (w.mat) {
      for (const k of org.g) w.mat[org.i * K + k]++;
      const taken = []; let ok = true;
      for (const k of G) { const at = takeNear(w, org.i, k, rnd); if (at < 0) { ok = false; break; } taken.push([at, k]); }
      if (!ok) { for (const [at, k] of taken) w.mat[at * K + k]++; for (const k of org.g) w.mat[org.i * K + k]--; org.tag = 2; continue; }
    }
    org.g = G.slice(); org.key = Sim.keyOf(org.g); org.hue = Sim.hueOf(org.key); w.hueAt[org.i] = org.hue;
    org.ip = 0; org.ctx = org.id; org.rh = 0; org.tag = toB ? 1 : 0;
  }
}
function snapshot(w, tick) {
  const start = { 0: A, 1: KIND === 'test' ? B : A };
  const row = { tick };
  for (const tag of [0, 1]) {
    const counts = new Map();
    let n = 0, e = 0, h = 0, same = 0;
    for (const o of w.alive) {
      if (o.tag !== tag) continue;
      n++; if (o.key.includes('e')) e++; if (o.key.includes('h')) h++; if (o.key === start[tag]) same++;
      counts.set(o.key, (counts.get(o.key) || 0) + 1);
    }
    row[tag] = { n, e: n ? e / n : 0, h: n ? h / n : 0, same: n ? same / n : 0, kinds: counts.size,
      top: [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 3) };
  }
  return row;
}

const results = [];
for (const seed of SEEDS) {
  const fpath = path.join(foundDir, `found_${cond}_${B}_t${T0}_s${String(seed).padStart(5, '0')}.json`);
  if (!fs.existsSync(fpath)) continue;
  const fz = JSON.parse(fs.readFileSync(fpath, 'utf8'));
  const rec = fz.arms.find(a => a.kind === KIND && a.mu === MU && a.rep === REP);
  if (!rec) continue;
  const opts = Object.assign({ seed, density: 8, mu: 0.01, light: 1 }, COND[cond]);
  const w0 = Sim.makeWorld(opts);
  while (w0.tick < T0 && w0.extinctAt < 0) Sim.tick(w0);
  const post = (seed * 7919 + 1 + REP * 104729) >>> 0, inj = (seed * 1000003 + 29) >>> 0;
  const w = cloneWorld(w0, post);
  w.o.mu = MU;
  convertAll(w, enc(A), enc(KIND === 'test' ? B : A), inj);
  const count = () => { let m1 = 0, m0 = 0; for (const o of w.alive) { if (o.tag === 1) m1++; else if (o.tag === 0) m0++; } return [w.alive.length, m1, m0]; };
  const series = [[0, ...count()]];
  const snaps = [snapshot(w, 0)];
  for (let t = 1; t <= T; t++) {
    Sim.tick(w);
    let stop = false;
    if (t % every === 0 || w.extinctAt >= 0) {
      const c = count(); series.push([t, ...c]);
      if (c[1] === 0 || c[2] === 0 || w.extinctAt >= 0) stop = true;
    }
    if (t % snapEvery === 0) snaps.push(snapshot(w, t));
    if (stop) break;
  }
  const match = JSON.stringify(series) === JSON.stringify(rec.series);
  results.push({ seed, match, snaps });
  const last = snaps[snaps.length - 1];
  const fmt = r => `n ${r.n} · e ${(r.e * 100).toFixed(0)}% · h ${(r.h * 100).toFixed(0)}% · 출발코드 ${(r.same * 100).toFixed(0)}% · 종류 ${r.kinds} · ${r.top.map(x => x.join(':')).join(' ')}`;
  console.log(`${cond} ${B} μ${MU} 시드 ${seed} · found 기록과 계통 수 ${match ? '일치' : '🚨 불일치'} · 틱 ${last.tick}`);
  console.log(`   씨앗 계통: ${fmt(last[0])}`);
  console.log(`   짝 계통:   ${fmt(last[1])}`);
}
const med = xs => { const s = xs.slice().sort((a, b) => a - b); return s.length ? s[(s.length - 1) >> 1] : NaN; };
for (const t of [500, 1000, 2000, 3000]) {
  const rows = results.map(r => r.snaps.find(s => s.tick === t)).filter(Boolean);
  if (!rows.length) continue;
  const line = tag => `e ${(med(rows.map(r => r[tag].e)) * 100).toFixed(0)}% · 출발코드 ${(med(rows.map(r => r[tag].same)) * 100).toFixed(0)}%`;
  console.log(`[중앙 ${rows.length}접시] 틱 ${t}: 씨앗 계통 ${line(0)} / 짝 계통 ${line(1)}`);
}
fs.mkdirSync(outDir, { recursive: true });
fs.writeFileSync(path.join(outDir, `diag_${cond}_${B}_mu${MU}_${KIND}.json`), JSON.stringify(results));
