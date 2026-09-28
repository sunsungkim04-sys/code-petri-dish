#!/usr/bin/env bash
# 해부 3A — 에너지 빈도 해부. 사전등록 해부3-사전등록-2026-09-16.md §2
#   배경: 해부 1 과 같은 energy 60 접시(끝 우세 reascld) · 팔 = ref + neutral 10 + 초점 코드 3 + 양성 대조 reascl
#   끼우는 비율 frac ∈ {0.01, 0.5} × 팔 μ ∈ {0, 0.01}
#   frac 0.1 은 이미 있는 자료를 쓴다 — 해부 1(dms_main · μ 0) · 해부 2A(dms_mu1 · μ 1%)
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
CODES=reashcld,reaschld,reasccld,reascl
L=freq.log; : > "$L"
J=freq_jobs.txt; : > "$J"
python3 - > "$J" <<'PYEOF'
import json
picked = []
for s in range(1, 101):
    z = json.load(open(f"main/energy_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "reascld":
        picked.append(s)
    if len(picked) == 60:
        break
for frac, ftag in (("0.01", "f0p01"), ("0.5", "f0p5")):
    for mu, mtag in (("0", "mu0"), ("0.01", "mu0p01")):
        for s in picked:
            print(frac, mu, f"dms_freq_{ftag}_{mtag}", s)
PYEOF
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 34 -n 4 sh -c 'node dms.js --cond energy --wt reascld --seed "$3" --frac "$0" --mu-assay "$1" --arms list --codes '"$CODES"' --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
