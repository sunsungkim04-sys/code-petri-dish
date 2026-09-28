#!/usr/bin/env bash
# 되감기 2 — 사전등록 되감기2-사전등록-2026-09-15.md §4
#   조건 space · energy · mat(밀도 8) × 시드 101~200 × 50,000틱 · sim v0.3.1 · μ 0.01  →  300접시
# 판정: rewind2_analyze.py(Q2' · Q3' · Q4') · Q1 재현은 main_analyze.py(같은 정의).
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=rewind2.log; : > "$L"
J=rewind2_jobs.txt; : > "$J"
for s in $(seq 101 200); do for c in space energy mat; do echo "$c 8 $s" >> "$J"; done; done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node replay.js --cond "$0" --density "$1" --seed "$2" --ticks 50000 --mu 0.01 --out rewind2' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
