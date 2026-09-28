#!/usr/bin/env bash
# 해부 4C — 에너지 구역 조성. 사전등록 해부4-사전등록-2026-09-16.md §4
#   본실험 에너지 생존 접시를 시드 순으로 20개 다시 키워(체크섬 대조) 빛 얼룩 · 어두운 곳의 조성을 잰다.
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=region.log; : > "$L"
J=region_jobs.txt; : > "$J"
python3 - > "$J" <<'PYEOF'
import json
n = 0
for s in range(1, 101):
    z = json.load(open(f"main/energy_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0:
        print(s)
        n += 1
    if n == 20:
        break
PYEOF
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 20 -n 1 sh -c 'node region.js --seed "$0" --out region' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
