#!/usr/bin/env node
/* 코드 배양접시 — 단일 코드 접시(monoculture) 계측기. 판정은 하지 않는다(mono_analyze.py).
 *
 * 접시에 고른 코드 **하나만** 넣고 돌린다(sim v0.3.2 의 ancestor 옵션). 본실험 접시를 쓰지 않는다 — 빈 접시에서 시작한다.
 * 기록: every 틱마다 [tick, pop, kinds, n_wt(출발 코드 그대로인 개체), births, deaths]
 *       끝: 상위 10 코드 · 체크섬 · 재료 보존 검사 · 멸종 틱
 * 사용: node mono.js --cond mat --code racld --seed 1 --mu 0.01 --ticks 20000 [--every 100] [--density 8] [--out dir]
 *
 * v1.1.0 (해부 5) — 기록에 **자유 글자 수**(아직 아무 개체도 안 쓴 재료)를 더했다. 읽기만 하므로 궤적은 바뀌지 않는다
 *                   — 해부 4A 의 한 칸을 다시 돌려 체크섬 · 앞 6열 동일을 확인한다(사전등록 점검 ④).
 */
'use strict';
const MONO_VERSION = '1.1.0';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
const COND = { space: { mat: false, energy: false }, mat: { mat: true, energy: false }, energy: { mat: false, energy: true } };
const cond = args.cond;
if (!COND[cond]) { console.error(`알 수 없는 cond: ${cond}`); process.exit(2); }
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다(ancestor 옵션)`); process.exit(2); }
const code = args.code;
const seed = parseInt(args.seed, 10);
const mu = parseFloat(args.mu);
const ticks = parseInt(args.ticks || '20000', 10);
const every = parseInt(args.every || '100', 10);
const density = parseFloat(args.density || '8');
const outDir = args.out || '.';
for (const ch of code) if (!Sim.LET.includes(ch)) { console.error(`코드에 없는 글자: ${ch}`); process.exit(2); }

const t0 = Date.now();
const opts = Object.assign({ seed, density, mu, light: 1, ancestor: code }, COND[cond]);
const w = Sim.makeWorld(opts);
const samples = [];
let consChecks = 0, consViol = 0;
function freeMat() {
  if (!w.mat) return null;
  let s = 0;
  for (let j = 0; j < w.mat.length; j++) s += w.mat[j];
  return s;
}
function record() {
  let wt = 0;
  const kinds = new Set();
  for (const o of w.alive) { kinds.add(o.key); if (o.key === code) wt++; }
  samples.push([w.tick, w.alive.length, kinds.size, wt, w.births, w.deaths, freeMat()]);
  if (w.mat) { consChecks++; if (Sim.countMat(w) !== w.matTotal0) consViol++; }
}
record();
while (w.tick < ticks && w.extinctAt < 0) {
  Sim.tick(w);
  if (w.tick % every === 0) record();
}
if (w.tick % every !== 0) record();
const counts = new Map();
for (const o of w.alive) counts.set(o.key, (counts.get(o.key) || 0) + 1);
const top = [...counts.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1)).slice(0, 10);
let freeByLet = null, tiedUp = 0;
if (w.mat) {
  freeByLet = new Array(Sim.K).fill(0);
  for (let i = 0; i < Sim.DISH.length; i++) {
    const c = Sim.DISH[i];
    for (let k = 0; k < Sim.K; k++) freeByLet[k] += w.mat[c * Sim.K + k];
  }
  for (const o of w.alive) { tiedUp += o.g.length; if (o.child) tiedUp += o.child.g.length; }
}
const out = {
  sim_version: Sim.VERSION, mono_version: MONO_VERSION, cond, code, seed, mu, ticks_planned: ticks,
  ticks_run: w.tick, every, opts, density: w.mat ? density : null, mat_total0: w.mat ? w.matTotal0 : null,
  extinct_at: w.extinctAt, births: w.births, deaths: w.deaths,
  conservation_checks: consChecks, conservation_violations: consViol, checksum: Sim.checksum(w),
  free_end: w.mat ? freeMat() : null, free_by_letter_end: freeByLet, tied_up_end: w.mat ? tiedUp : null,
  letters: Sim.LET, dish_cells: Sim.DISH.length,
  header: ['tick', 'pop', 'kinds', 'n_wt', 'births', 'deaths', 'free'], samples, top, elapsed_ms: Date.now() - t0
};
fs.mkdirSync(outDir, { recursive: true });
const tag = `mono_${cond}_${code}_mu${String(mu).replace('.', 'p')}_s${String(seed).padStart(5, '0')}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
const last = samples[samples.length - 1];
const freeTxt = w.mat ? ` · free ${out.free_end}/${w.matTotal0}` : ' · free -';
console.log(`DONE ${tag} · v${Sim.VERSION}/m${MONO_VERSION} · tick ${w.tick} · extinct ${w.extinctAt} · pop ${last[1]} · wt ${last[3]} · kinds ${last[2]}${freeTxt} · cons ${consViol}/${consChecks} · chk ${out.checksum} · ${((Date.now() - t0) / 1000).toFixed(1)}s`);
