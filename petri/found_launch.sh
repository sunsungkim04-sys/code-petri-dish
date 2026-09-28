#!/usr/bin/env bash
# 해부 2B — 씨앗 시절 창업 경쟁. 사전등록 해부2-사전등록-2026-09-15.md §3
#   조건마다 본실험 시드 1~100 접시를 T0 까지 다시 키워(본실험 기록과 대조) 살아 있는 개체 전부를 반씩 씨앗 · 짝 코드로 바꾼다.
#   space T0 500 · 짝 reascld(양성 대조) / energy T0 500 · 짝 rahcld(양성 대조) / mat T0 2000 · 짝 racld(시험)
#   팔 = μ {0, 0.01} × 반복 5 × {중립, 시험} = 20 · 3,000틱
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=found.log; : > "$L"
J=found_jobs.txt; : > "$J"
for s in $(seq 1 100); do
  echo "space 500 reascld $s" >> "$J"
  echo "energy 500 rahcld $s" >> "$J"
  echo "mat 2000 racld $s" >> "$J"
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 40 -n 4 sh -c 'node found.js --cond "$0" --t0 "$1" --b "$2" --seed "$3" --out found' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
