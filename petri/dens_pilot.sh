#!/usr/bin/env bash
# 해부 5 설계용 탐색 발사 — 판정 아님. 밀도 사다리 + space 가 살아남는지, 재료가 바닥나는지.
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
L=dens_pilot.log; : > "$L"
J=dens_pilot_jobs.txt; : > "$J"
for code in racld rascld; do
  for d in 4 8 14 20 32; do
    for s in 1 2 3; do echo "mat $d $code $s densp/mat_d$(printf %02d $d)" >> "$J"; done
  done
  for s in 1 2 3; do echo "space 8 $code $s densp/space" >> "$J"; done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 36 -n 5 sh -c 'node mono.js --cond "$0" --density "$1" --code "$2" --seed "$3" --mu 0 --ticks 20000 --out "$4"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
