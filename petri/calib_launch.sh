#!/usr/bin/env bash
# 사건 문턱 고정 — 보정 전용 접시. 사전등록 사건문턱-고정-사전등록-2026-09-15.md
#   조건 space · energy · mat(밀도 8) × 시드 301~500 × 50,000틱 · sim v0.3.1 · μ 0.01  →  600접시
#   이 접시들은 θ(생존 접시 L2_x 95 백분위)를 재는 데만 쓴다 — 사건 판정에는 쓰지 않는다.
#   시드 201~300 은 되감기 2 설계를 바꿀 때 쓰려고 비워 둔다.
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=calib.log; : > "$L"
J=calib_jobs.txt; : > "$J"
for s in $(seq 301 500); do for c in space energy mat; do echo "$c 8 $s" >> "$J"; done; done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 56 -n 3 sh -c 'node replay.js --cond "$0" --density "$1" --seed "$2" --ticks 50000 --mu 0.01 --out calib' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
