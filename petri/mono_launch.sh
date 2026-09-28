#!/usr/bin/env bash
# 해부 4 A · B — 단일 코드 접시. 사전등록 해부4-사전등록-2026-09-16.md §2 · §3
#   curve     : mat · racld/rascld × μ {0,0.0025,0.005,0.01,0.02,0.04} × 시드 1~30 × 20,000틱
#   neighbors : mat · 두 WT 의 한 글자 이웃 전부(겹치면 한 번) · μ 0 × 시드 1~3 × 4,000틱
#   neighbors10 : 같은 이웃을 시드 4~10 으로 더 (사전등록 규칙이 '보류' 를 냈을 때의 보강 — 판정선은 그대로)
set -u
MODE=${1:?curve · neighbors · neighbors10}
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=mono_$MODE.log; : > "$L"
J=mono_${MODE}_jobs.txt; : > "$J"
if [ "$MODE" = curve ]; then
  for code in racld rascld; do for mu in 0 0.0025 0.005 0.01 0.02 0.04; do for s in $(seq 1 30); do
    echo "$code $mu $s 20000 mono_curve" >> "$J"
  done; done; done
else
  MONO_MODE=$MODE python3 - > "$J" <<'PYEOF'
LET = "nsracldjehx"
def neighbours(wt):
    out = set()
    for p in range(len(wt)):
        for ch in LET:
            if ch != wt[p]:
                out.add(wt[:p] + ch + wt[p+1:])
    for p in range(len(wt)):
        out.add(wt[:p] + wt[p+1:])
    for p in range(len(wt) + 1):
        for ch in LET:
            out.add(wt[:p] + ch + wt[p:])
    return {k for k in out if 2 <= len(k) <= 48}
import os
SEEDS = (1, 2, 3) if os.environ.get("MONO_MODE") == "neighbors" else (4, 5, 6, 7, 8, 9, 10)
codes = neighbours("racld") | neighbours("rascld") | {"racld", "rascld"}
for c in sorted(codes):
    for s in SEEDS:
        print(c, 0, s, 4000, "mono_nb")
PYEOF
fi
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 5 sh -c 'node mono.js --cond mat --code "$0" --mu "$1" --seed "$2" --ticks "$3" --out "$4"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
