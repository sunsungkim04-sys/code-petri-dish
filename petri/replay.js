#!/usr/bin/env node
/* 코드 배양접시 — 접시 하나를 화면 없이 돌려 기록만 남기는 계측기.
 * 판정은 하지 않는다(분석은 별도 스크립트). 기록하는 것:
 *   samples   — every 틱마다 [tick, 개체, 코드종류, 평균길이, 빌려쓰는중, j포함, e포함, h포함, x포함, c없음, 글자11종 합계…]
 *   tops      — topEvery 틱마다 상위 20 코드 [코드, 개체수]
 *   specials  — topEvery 틱마다 조건별 '가장 큰 단일 코드' 를 전수에서: j 포함 · e 포함 · h 포함 · n 포함 · x 포함 · c 없음
 *   snap      — 계획 틱 − 2,000 시점의 코드 전수 (사전등록 Q2 기준자)
 *   final     — 마지막 틱의 코드 전수 · 체크섬 · 재료 보존 위반 횟수 · 소요 시간
 * 사용: node replay.js --cond mat --seed 1 --ticks 50000 [--density 8] [--mu 0.01] [--every 50] [--top-every 500] [--out dir]
 *   cond = space(자리만) | mat(자리+재료) | energy(자리+에너지) | both(자리+재료+에너지)
 *   density 는 재료가 켜진 조건에서만 의미가 있다(칸당 평균 글자 수).
 * 해부 8(09-17 추가): --find-first 1 이면 sim v0.3.3 의 재료 먼저 세계로 처음부터 돌린다(파일 이름 끝 _ff · 출력에 find_first).
 *   기본값(끔)이면 출력이 전과 같다.
 * 해부 17(09-28 추가): --loops 1 이면 topEvery 틱마다 살아 있는 개체의 **복사 고리 길이 히스토그램**을 전수에서 센다(읽기만 · 난수 소비 없음 → 체크섬 그대로).
 *   loops — [tick, {고리 길이: 개체 수}, 고리 없는 개체 수] · 고리 길이는 pairs13_analyze.py 의 loop_ticks 와 같은 정의(라벨 s 뒤부터 l 까지 · x 는 0 틱).
 *   출력에 loops · loops_on 이 더해지고 파일 이름은 그대로다. 기본값(끔)이면 출력이 전과 같다.
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
  energy: { mat: false, energy: true },
  both: { mat: true, energy: true }
};
const cond = args.cond || 'mat';
if (!COND[cond]) { console.error(`알 수 없는 cond: ${cond}`); process.exit(2); }
const seed = parseInt(args.seed || '1', 10);
const ticks = parseInt(args.ticks || '50000', 10);
const every = parseInt(args.every || '50', 10);
const topEvery = parseInt(args['top-every'] || '500', 10);
const density = parseFloat(args.density || '8');
const mu = parseFloat(args.mu || '0.01');
const outDir = args.out || '.';
const findFirst = args['find-first'] === '1';
const loopsOn = args.loops === '1';
function loopTicks(code) {
  if (!code.includes('l') || !code.includes('c')) return null;
  const li = code.indexOf('l');
  const si = code.lastIndexOf('s', li - 1);
  const start = si >= 0 ? si + 1 : 0;
  if (start > li) return null;
  const seg = code.slice(start, li + 1);
  if (!seg.includes('c')) return null;
  let n = 0; for (const ch of seg) if (ch !== 'x') n++;
  return n;
}
const loopCache = new Map();
function loopOf(key) {
  if (!loopCache.has(key)) loopCache.set(key, loopTicks(key));
  return loopCache.get(key);
}
if (findFirst && !['0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error('--find-first 는 sim 0.3.3~0.3.5 가 필요하다'); process.exit(2); }
const opts = Object.assign({ seed, density, mu, light: 1 }, COND[cond]);
if (findFirst) opts.findFirst = true;
const snapTick = ticks - 2000;

const HEADER = ['tick', 'pop', 'kinds', 'mean_len', 'borrowing', 'has_j', 'has_e', 'has_h', 'has_x', 'no_c', ...[...Sim.LET].map(ch => `n_${ch}`)];
const SPECIAL = ['j', 'e', 'h', 'n', 'x', 'noc'];
const t0 = Date.now();
const w = Sim.makeWorld(opts);
const samples = [], tops = [], specials = [], loops = [];
let snap = null, conservationChecks = 0, conservationViolations = 0, orgSteps = 0;

const byCount = (a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1);
function countsOf() {
  const counts = new Map();
  for (const o of w.alive) counts.set(o.key, (counts.get(o.key) || 0) + 1);
  return counts;
}
function record() {
  const n = w.alive.length;
  let len = 0, borrow = 0, hj = 0, he = 0, hh = 0, hx = 0, noc = 0;
  const letters = new Array(Sim.K).fill(0);
  const kinds = new Set();
  for (const o of w.alive) {
    len += o.g.length;
    if (o.ctx !== o.id) borrow++;
    for (const k of o.g) letters[k]++;
    const key = o.key;
    kinds.add(key);
    if (key.includes('j')) hj++;
    if (key.includes('e')) he++;
    if (key.includes('h')) hh++;
    if (key.includes('x')) hx++;
    if (!key.includes('c')) noc++;
  }
  samples.push([w.tick, n, kinds.size, n ? +(len / n).toFixed(3) : 0, borrow, hj, he, hh, hx, noc, ...letters]);
  if (w.mat) { conservationChecks++; if (Sim.countMat(w) !== w.matTotal0) conservationViolations++; }
}
function recordTop() {
  const counts = countsOf();
  tops.push([w.tick, [...counts.entries()].sort(byCount).slice(0, 20)]);
  if (loopsOn) {
    const hist = {}; let none = 0;
    for (const [key, cnt] of counts) { const L = loopOf(key); if (L === null) none += cnt; else hist[L] = (hist[L] || 0) + cnt; }
    loops.push([w.tick, hist, none]);
  }
  const best = SPECIAL.map(() => null);
  for (const [key, cnt] of counts) {
    for (let s = 0; s < SPECIAL.length; s++) {
      const tag = SPECIAL[s];
      const hit = tag === 'noc' ? !key.includes('c') : key.includes(tag);
      if (hit && (!best[s] || cnt > best[s][1] || (cnt === best[s][1] && key < best[s][0]))) best[s] = [key, cnt];
    }
  }
  specials.push([w.tick, ...best]);
}

record(); recordTop();
while (w.tick < ticks && w.extinctAt < 0) {
  orgSteps += w.alive.length;
  Sim.tick(w);
  if (w.tick % every === 0) record();
  if (w.tick % topEvery === 0) recordTop();
  if (w.tick === snapTick) snap = [w.tick, [...countsOf().entries()].sort(byCount)];
}
if (w.tick % every !== 0) record();

const out = {
  sim_version: Sim.VERSION, cond, opts, ticks_planned: ticks, ticks_run: w.tick,
  extinct_at: w.extinctAt, births: w.births, deaths: w.deaths,
  mat_total0: w.matTotal0, conservation_checks: conservationChecks, conservation_violations: conservationViolations,
  checksum: Sim.checksum(w), org_steps: orgSteps, elapsed_ms: Date.now() - t0,
  header: HEADER, samples, tops, special_header: SPECIAL, specials,
  snap, final_counts: [...countsOf().entries()].sort(byCount)
};
if (loopsOn) { out.loops = loops; out.loops_on = true; }
if (findFirst) out.find_first = true;
fs.mkdirSync(outDir, { recursive: true });
const tag = `${cond}_d${String(density).replace('.', 'p')}_mu${String(mu).replace('.', 'p')}_s${String(seed).padStart(5, '0')}${findFirst ? '_ff' : ''}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
const last = samples[samples.length - 1];
console.log([
  `DONE ${tag}`, `v${Sim.VERSION}`, `tick ${w.tick}`, `extinct ${w.extinctAt}`, `pop ${last[1]}`, `kinds ${last[2]}`,
  `len ${last[3]}`, `cons ${conservationViolations}/${conservationChecks}`,
  `chk ${out.checksum}`, `${(out.elapsed_ms / 1000).toFixed(1)}s`, `${(orgSteps / Math.max(1, out.elapsed_ms) * 1000 / 1e6).toFixed(2)}M steps/s`
].join(' · '));
