#!/usr/bin/env bash
# 해부 8 A — 재료 먼저 세계(sim v0.3.3 findFirst)에서 처음부터 되감기. 사전등록 해부8-사전등록-2026-09-17.md §2
#   mat 밀도 8 · μ 1% · 50,000틱 · 시드 101~200 (되감기 2 mat 과 같은 시드 = 같은 첫 재료 배치) → 100 접시
#   판정은 rewind8_analyze.py (원래 규칙 쪽은 되감기 2 기록 rewind2/).
# 마른 실행: PETRI_SEEDS="101 102" PETRI_ROOT=smoke8a PETRI_TICKS=2000 ./rw8_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[345]'" sim.js || { echo "sim.js 가 v0.3.3~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "find-first" replay.js || { echo "replay.js 에 --find-first 가 없다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 101 200)}
ROOT=${PETRI_ROOT:-rewind8ff}
TICKS=${PETRI_TICKS:-50000}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for s in $SEEDS; do echo "$s $TICKS $ROOT" >> "$J"; done
[ "$(awk 'NF!=3' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 3 이 아니다 — 발사 중단"; exit 1; }
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 50 -n 3 sh -c 'node replay.js --cond mat --density 8 --mu 0.01 --seed "$0" --ticks "$1" --find-first 1 --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
