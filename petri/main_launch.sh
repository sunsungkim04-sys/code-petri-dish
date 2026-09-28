#!/usr/bin/env bash
# 본실험 — 테이프 되감기. 사전등록 §4 확정값:
#   조건 space · energy · mat(밀도 8) × 시드 1~100 × 50,000틱 · sim v0.3.0 · μ 0.01  →  300접시
# space · energy 에서 밀도는 의미가 없다(파일명 표기만 d8). 판정은 main_analyze.py.
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
# v0.3.1 은 계통 표지만 더했고 궤적이 같다(dms.js 가 본실험 체크섬과 대조) — 둘 다 받는다.
grep -qE "VERSION = '0\.3\.[01]'" sim.js || { echo "sim.js 가 v0.3.0/0.3.1 이 아니다 — 발사 중단"; exit 1; }
L=main.log; : > "$L"
J=main_jobs.txt; : > "$J"
for s in $(seq 1 100); do for c in space energy mat; do echo "$c 8 $s" >> "$J"; done; done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node replay.js --cond "$0" --density "$1" --seed "$2" --ticks 50000 --mu 0.01 --out main' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
