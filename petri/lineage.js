#!/usr/bin/env node
/* 코드 배양접시 — 누가 누구를 낳았나(계보 추적) 계측기. 판정은 하지 않는다(lineage_analyze.py).
 *
 * 단일 코드 접시(sim v0.3.2 ancestor)를 복제 오류 μ 로 돌리며 **모든 출생의 부모**를 기록한다.
 * sim.js 는 바꾸지 않는다: 틱마다 살아 있는 개체의 표지(tag)를 자기 id 로 덮어쓰면, 그 틱에 태어난 자식은
 * spawn 에서 부모의 표지 = 부모 id 를 물려받는다. tag 는 sim 의 동작과 체크섬에 쓰이지 않는다(궤적 불변 —
 * 같은 설정의 mono.js 기록과 체크섬 대조로 확인한다).
 *
 * 개체마다: 태어난 틱 · 태어날 때 출발 코드(WT)였나 · 깊이(WT 로부터 몇 대째 변종인가: WT=0,
 *           WT 부모의 변종 자식=1, 깊이 d 변종의 변종 자식=d+1) · 낳은 자식 수
 * 부모의 WT 여부는 **낳는 순간의 코드**로 본다(살면서 우주선 돌연변이로 바뀔 수 있다).
 * 창 [win0, T−700] 에 태어난 개체만 센다 — 최대 수명 599틱이라 T 에는 모두 죽어 자식 수가 확정돼 있다.
 *
 * v1.5.0 (해부 10): --cosmic 0 · --diffuse X · --age0 N --age-var N (sim 0.3.5 옵션을 그대로 넘긴다 · 파일 이름 끝 _c0 · _dfX · _aN · 출력 opts).
 *         v1.4.0 의 _rd 꼬리는 실제로는 붙지 않았다(꼬리 패턴이 안 맞아 교체가 빠짐 — 해부 9 파일은 lin_<code>_mu…_s….json) · v1.5.0 부터 붙는다. 끄면 v1.4.0 과 출력이 같다.
 * v1.4.0 (해부 9 A): --remember-die 1 이면 sim v0.3.4 의 기억하는 주사위 세계에서 돌린다(출력 remember_die).
 *         끄면 v1.3.0 과 출력이 같다.
 * v1.3.0 (해부 8 C): --spec 1 이면 WT 코드 개체가 베끼기를 실행한 틱마다 결과를 틱 전후 상태로 가른다(읽기만 한다):
 *         맞게 베낌 · 바뀜(들어간 글자별) · 끼어듦(글자별) · 빠뜨림 · 헛손질 · 그 밖(자식이 사라짐 등).
 *         글자 하나 나아갈 때마다(맞게 · 바뀜 · 빠뜨림) 오류가 몇 번 들어가나를 세어 부풀림을 주사위 횟수 몫과 거름 몫으로 가른다.
 *         끄면 v1.2.0 과 출력 필드가 같다(더한 필드 제외).
 * v1.2.0 (해부 7): --find-first 1 이면 sim v0.3.3 의 재료 먼저 세계에서 돌린다(sim 0.3.3 필요).
 *         --watch CODE 면 그 코드의 출생을 부모 종류(같은 코드 · 출발 코드 · 그 밖)로 나눠 세고,
 *         그 코드로 태어난 개체의 자식 수(창 안) · 처음 새로 생긴 틱 · 500틱마다 개체 수를 적는다.
 *         창 안 개체의 '자기와 같은 코드로 태어난 자식 수' 도 적는다(same_kids_wt · same_kids_watch).
 *         둘 다 끄면 v1.1.0 과 출력 필드가 같다(더한 필드 제외) — 읽기만 하므로 궤적도 같다.
 * v1.1.0: 끝 상위 코드 8개를 기록한다(씨앗 접시가 변종 구름으로 녹았나, `racld` 에 먹혔나를 가르려고). 읽기만 한다.
 *
 * 사용: node lineage.js --code racld --seed 601 --mu 0.01 [--ticks 20000] [--win0 5000] [--density 8] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');
const LINEAGE_VERSION = '1.5.0';

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다`); process.exit(2); }
const code = args.code;
const seed = parseInt(args.seed, 10);
const mu = parseFloat(args.mu);
const T = parseInt(args.ticks || '20000', 10);
const win0 = parseInt(args.win0 || '5000', 10);
const win1 = T - 700;
const density = parseFloat(args.density || '8');
const outDir = args.out || '.';
const findFirst = args['find-first'] === '1';
const watch = args.watch || null;
const spec = args.spec === '1';
const rememberDie = args['remember-die'] === '1';
if (rememberDie && !['0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error('--remember-die 는 sim 0.3.4/0.3.5 가 필요하다'); process.exit(2); }
const cosmicOff = args.cosmic === '0';
const diffuseV = args.diffuse ? parseFloat(args.diffuse) : null;
const age0V = args.age0 ? parseInt(args.age0, 10) : null;
const ageVarV = args['age-var'] ? parseInt(args['age-var'], 10) : null;
if ((cosmicOff || diffuseV !== null || age0V !== null || ageVarV !== null) && Sim.VERSION !== '0.3.5') { console.error('--cosmic/--diffuse/--age0 은 sim 0.3.5 가 필요하다'); process.exit(2); }
const SP = { ok: 0, del: 0, stall: 0, other: 0, dead: 0, borrow: 0, sub: new Array(Sim.K).fill(0), ins: new Array(Sim.K).fill(0) };
const codeArr = [...code].map(ch => Sim.LET.indexOf(ch));
let spPre = [], spAge = [], spCh = [], spLen = [], spRh = [], spLet = [];
if (findFirst && !['0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error('--find-first 는 sim 0.3.3~0.3.5 가 필요하다'); process.exit(2); }
if (watch) for (const ch of watch) if (!Sim.LET.includes(ch)) { console.error(`watch 에 없는 글자: ${ch}`); process.exit(2); }
for (const ch of code) if (!Sim.LET.includes(ch)) { console.error(`코드에 없는 글자: ${ch}`); process.exit(2); }

const t0 = Date.now();
const wopts = { seed, density, mu, light: 1, ancestor: code, mat: true, energy: false };
if (findFirst) wopts.findFirst = true;
if (rememberDie) wopts.rememberDie = true;
if (cosmicOff) wopts.cosmic = false;
if (diffuseV !== null) wopts.diffuse = diffuseV;
if (age0V !== null) wopts.age0 = age0V;
if (ageVarV !== null) wopts.ageVar = ageVarV;
const w = Sim.makeWorld(wopts);

let cap = 1 << 20;
let bornAt = new Int32Array(cap), depth = new Int16Array(cap), kids = new Int32Array(cap);
let bornWT = new Uint8Array(cap), curWT = new Uint8Array(cap), sameKids = new Int32Array(cap);
function grow(need) {
  while (cap <= need) cap *= 2;
  const g = (A, C) => { const b = new C(cap); b.set(A); return b; };
  bornAt = g(bornAt, Int32Array); depth = g(depth, Int16Array); kids = g(kids, Int32Array);
  bornWT = g(bornWT, Uint8Array); curWT = g(curWT, Uint8Array); sameKids = g(sameKids, Int32Array);
  if (bornWatch) { bornWatch = g(bornWatch, Uint8Array); curWatch = g(curWatch, Uint8Array); curKey.length = cap; }
}
for (const o of w.alive) { bornAt[o.id] = 0; depth[o.id] = 0; bornWT[o.id] = o.key === code ? 1 : 0; }

let wtParentBirths = 0, wtParentMutBirths = 0;     // WT 부모의 출생 · 그중 변종
let wtBirths = 0, wtFromMut = 0;                    // WT 로 태어난 출생 · 그중 변종 부모(되돌아옴)
let seen = 0, missed = 0, cosmicWT = 0;
// --watch: 그 코드의 출생을 부모 종류로 · 그 코드로 태어났나 · 창 안 자식 수
const W_PAR = { same: 0, wt: 0, other: 0 };
const W_PAR_OTHER = new Map();                     // 그 밖 부모 코드 → 횟수(상위만 적는다)
let watchFirst = -1;
let bornWatch = watch ? new Uint8Array(cap) : null;
let curWatch = watch ? new Uint8Array(cap) : null;
let curKey = watch ? new Array(cap) : null;
const series = [];
function snap() {
  let nW = 0, n1 = 0, n2 = 0, n3 = 0, nX = 0;
  for (const o of w.alive) {
    if (o.key === code) nW++;
    else { const d = depth[o.id]; if (d <= 1) n1++; else if (d === 2) n2++; else n3++; }
    if (watch && o.key === watch) nX++;
  }
  if (watch) series.push([w.tick, w.alive.length, nW, n1, n2, n3, nX]);
  else series.push([w.tick, w.alive.length, nW, n1, n2, n3]);
}
snap();
while (w.tick < T && w.extinctAt < 0) {
  for (const o of w.alive) {
    o.tag = o.id;
    const isW = o.key === code ? 1 : 0;
    if (!isW && bornWT[o.id] && curWT[o.id]) cosmicWT++;   // WT 로 태어나 살면서 바뀐 개체(처음 바뀐 틱에 한 번)
    curWT[o.id] = isW;
    if (watch) { curWatch[o.id] = o.key === watch ? 1 : 0; curKey[o.id] = o.key; }
  }
  if (spec) {
    spPre = []; spAge = []; spCh = []; spLen = []; spRh = []; spLet = [];
    for (const o of w.alive) {
      if (o.key !== code || !o.child || o.rh >= o.g.length) continue;
      spPre.push(o); spAge.push(o.age); spCh.push(o.child); spLen.push(o.child.g.length); spRh.push(o.rh); spLet.push(o.g[o.rh]);
    }
  }
  const b0 = w.births;
  Sim.tick(w);
  if (spec) {
    for (let q = 0; q < spPre.length; q++) {
      const o = spPre[q];
      if (o.age !== spAge[q] + 1) continue;                 // 차례가 오기 전에 죽음 — 실행 안 함
      if (o.lastHost !== o.id) { SP.borrow++; continue; }
      if (codeArr[o.last] !== 4) continue;                    // 베끼기가 아닌 명령 (틱 시작 코드 = WT 로 읽는다)
      if (o.dead) { SP.dead++; continue; }
      const c = o.child;
      if (c !== spCh[q]) { SP.other++; continue; }
      const dl = c.g.length - spLen[q], dr = o.rh - spRh[q];
      if (dl === 1 && dr === 1) { const k = c.g[spLen[q]]; if (k === spLet[q]) SP.ok++; else SP.sub[k]++; }
      else if (dl === 1 && dr === 0) SP.ins[c.g[spLen[q]]]++;
      else if (dl === 0 && dr === 1) SP.del++;
      else if (dl === 0 && dr === 0) SP.stall++;
      else SP.other++;
    }
  }
  const bt = w.tick - 1;
  let here = 0;
  for (const o of w.alive) {
    if (o.bornAt !== bt) continue;
    here++;
    if (o.id >= cap) grow(o.id);
    const p = o.tag;
    kids[p]++;
    const isW = o.key === code ? 1 : 0;
    const pW = curWT[p];
    bornAt[o.id] = bt; bornWT[o.id] = isW; curWT[o.id] = isW;
    depth[o.id] = isW ? 0 : (pW ? 1 : Math.min(depth[p] + 1, 30000));
    if (watch) {
      const isX = o.key === watch ? 1 : 0;
      bornWatch[o.id] = isX; curWatch[o.id] = isX; curKey[o.id] = o.key;
      if (isX && curWatch[p]) sameKids[p]++;
      if (isX) {
        if (curWatch[p]) W_PAR.same++;
        else {
          if (watchFirst < 0) watchFirst = bt;
          if (pW) W_PAR.wt++;
          else { W_PAR.other++; W_PAR_OTHER.set(curKey[p], (W_PAR_OTHER.get(curKey[p]) || 0) + 1); }
        }
      }
    }
    if (pW) { wtParentBirths++; if (!isW) wtParentMutBirths++; else sameKids[p]++; }
    if (isW) { wtBirths++; if (!pW) wtFromMut++; }
  }
  seen += here;
  missed += (w.births - b0) - here;
  if (w.tick % 500 === 0) snap();
}
if (w.tick % 500 !== 0) snap();

// 창 안 개체의 자식 수 — 종류별
const cls = { wt: [0, 0, 0], d1: [0, 0, 0], d2: [0, 0, 0], d3: [0, 0, 0] };   // [개체 수, 자식 합, 자식 ≥1 개체 수]
const hist1 = new Array(8).fill(0), histW = new Array(8).fill(0);
const clsWatch = [0, 0, 0];
let sameWT = 0, sameWatch = 0;
const maxId = Math.min(w.nextId, cap);             // 놓친 출생(태어난 틱에 먹힌 개체)은 기록이 없다
for (let id = 0; id < maxId; id++) {
  const b = bornAt[id];
  if (id >= 1 && b === 0) continue;             // 기록 안 된 id(창 밖 · 놓친 출생)
  if (b < win0 || b > win1) continue;
  const k = kids[id];
  const c = bornWT[id] ? cls.wt : depth[id] <= 1 ? cls.d1 : depth[id] === 2 ? cls.d2 : cls.d3;
  c[0]++; c[1] += k; if (k > 0) c[2]++;
  if (bornWT[id]) histW[Math.min(k, 7)]++; else if (depth[id] <= 1) hist1[Math.min(k, 7)]++;
  if (bornWT[id]) sameWT += sameKids[id];
  if (watch && bornWatch[id]) { clsWatch[0]++; clsWatch[1] += k; if (k > 0) clsWatch[2]++; sameWatch += sameKids[id]; }
}
const endCounts = new Map();
for (const o of w.alive) endCounts.set(o.key, (endCounts.get(o.key) || 0) + 1);
const topEnd = [...endCounts.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1)).slice(0, 8);
const out = {
  lineage_version: LINEAGE_VERSION, remember_die: rememberDie, opts: wopts, spec: spec ? SP : null, top_end: topEnd, kinds_end: endCounts.size, find_first: findFirst,
  watch, watch_births_by_parent: watch ? W_PAR : null, watch_first_denovo_tick: watch ? watchFirst : null,
  watch_other_parents_top: watch ? [...W_PAR_OTHER.entries()].sort((a, b) => b[1] - a[1]).slice(0, 8) : null,
  watch_class: watch ? clsWatch : null, same_kids_wt: sameWT, same_kids_watch: watch ? sameWatch : null, sim_version: Sim.VERSION, code, seed, mu, density, ticks_planned: T, ticks_run: w.tick,
  win: [win0, win1], extinct_at: w.extinctAt, births: w.births, seen_births: seen, missed_births: missed,
  cosmic_wt_changed: cosmicWT, checksum: Sim.checksum(w),
  wt_parent_births: wtParentBirths, wt_parent_mut_births: wtParentMutBirths, wt_births: wtBirths, wt_from_mut: wtFromMut,
  classes: cls, kids_hist_wt: histW, kids_hist_d1: hist1,
  series_header: watch ? ['tick', 'pop', 'n_wt', 'n_d1', 'n_d2', 'n_d3plus', 'n_watch'] : ['tick', 'pop', 'n_wt', 'n_d1', 'n_d2', 'n_d3plus'], series, elapsed_ms: Date.now() - t0
};
fs.mkdirSync(outDir, { recursive: true });
const tag = `lin_${code}_mu${String(mu).replace('.', 'p')}_s${String(seed).padStart(5, '0')}${findFirst ? '_ff' : ''}${rememberDie ? '_rd' : ''}${cosmicOff ? '_c0' : ''}${diffuseV !== null ? '_df' + String(diffuseV).replace('.', 'p') : ''}${age0V !== null ? '_a' + age0V : ''}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
const kw = cls.wt[0] ? cls.wt[1] / cls.wt[0] : NaN, k1 = cls.d1[0] ? cls.d1[1] / cls.d1[0] : NaN;
console.log(`DONE ${tag} · v${Sim.VERSION}/l${LINEAGE_VERSION}${findFirst ? ' · 재료먼저' : ''}${rememberDie ? ' · 기억주사위' : ''}${watch ? ` · watch ${watch} 새로 ${W_PAR.wt + W_PAR.other} 물림 ${W_PAR.same}` : ''} · tick ${w.tick} · ext ${w.extinctAt} · births ${w.births} seen ${seen} missed ${missed} · kWT ${kw.toFixed(3)} k1 ${k1.toFixed(3)} · chk ${out.checksum} · ${((Date.now() - t0) / 1000).toFixed(1)}s`);
