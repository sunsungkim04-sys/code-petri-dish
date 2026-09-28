#!/usr/bin/env bash
# 해부 17 — 진화 궤적. 사전등록 해부17-사전등록-2026-09-28.md §1~§2
#   replay.js --cond mat --density 8 --ticks 50000 --loops 1 · μ {0.001 0.003 0.005 0.01} × 시드 1~20 → 80
#   μ 0.01 · 시드 1~20 은 main/ 과 같은 설정 → 체크섬이 같아야 한다(G0 · 판정기가 대조)
# 마른 실행(파일럿): PETRI_SEEDS="901 902" PETRI_MUS="0.001 0.01" PETRI_ROOT=pilot17 ./evo17_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "loopsOn" replay.js || { echo "replay.js 에 --loops 가 없다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 1 20)}
MUS=${PETRI_MUS:-"0.001 0.003 0.005 0.01"}
PAR=${PETRI_PAR:-60}
ROOT=${PETRI_ROOT:-evo17}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for mu in $MUS; do for s in $SEEDS; do echo "$mu $s $ROOT" >> "$J"; done; done
echo "시드 $SEEDS · μ $MUS · 작업 $(wc -l < "$J") · $(date '+%F %T') · uptime -s $(uptime -s)" | tr '\n' ' ' >> "$L"; echo >> "$L"
t0=$(date +%s)
xargs -P "$PAR" -n 3 sh -c 'node replay.js --cond mat --density 8 --mu "$0" --seed "$1" --ticks 50000 --loops 1 --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
