#!/usr/bin/env bash
# 이긴 코드 해부(가상 DMS) — 설계 노트 해부-DMS-설계-2026-09-15.md §5 · §7
#   pilot : WT 마다 배경 1개(규칙상 첫 접시) × ref + neutral 10 + d 빠뜨림
#   main  : 배경 space 20 · mat 20 · energy 60 × ref + neutral 10 + 한 글자 돌연변이 전부
#           (energy 60 은 09-15 파일럿 중립 팔 흔들림으로 증원 — 설계 노트 §7)
# 배경 규칙: 본실험 생존 접시 중 끝 우세 코드가 WT 인 접시를 시드 순으로 앞에서부터.
set -u
#   mu1   : 해부 2A — main 과 같은 배경 · 팔 · 돌연변이, 팔을 돌리는 동안 μ 1% (해부2-사전등록-2026-09-15.md)
MODE=${1:?pilot · main · mu1}
EXTRA=""
case "$MODE" in pilot) ARMS=pilot ;; main) ARMS=all ;; mu1) ARMS=all; EXTRA="--mu-assay 0.01" ;; *) echo "모드: pilot | main | mu1"; exit 2 ;; esac
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=dms_$MODE.log; : > "$L"
J=dms_${MODE}_jobs.txt
python3 - "$MODE" > "$J" <<'EOF'
import json, sys
NBG = {"space": 20, "mat": 20, "energy": 60} if sys.argv[1] in ("main", "mu1") else {"space": 1, "mat": 1, "energy": 1}
for cond, wt in (("space", "reasccld"), ("mat", "racld"), ("energy", "reascld")):
    nbg = NBG[cond]
    picked = []
    for s in range(1, 101):
        z = json.load(open(f"main/{cond}_d8_mu0p01_s{s:05d}.json"))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == wt:
            picked.append(s)
        if len(picked) == nbg:
            break
    for s in picked:
        print(cond, wt, s)
EOF
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 64 -n 3 sh -c 'node dms.js --cond "$0" --wt "$1" --seed "$2" --arms '"$ARMS"' '"$EXTRA"' --out dms_'"$MODE" < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
