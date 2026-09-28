#!/usr/bin/env bash
# 해부 3B — 평평함 직접 재기. 사전등록 해부3-사전등록-2026-09-16.md §3
#   씨앗 rascld 의 한 글자 이웃 전부를, racld 를 쟀던 것과 **같은 mat 배경 20접시 · 같은 설정**(frac 0.1 · μ 0 · 3,000틱)에서 잰다.
#   배경 규칙은 해부 1 그대로(끝 우세 코드가 racld 인 접시를 시드 순 앞에서 20) — WT 만 rascld 로 바꾼다.
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[12345]'" sim.js || { echo "sim.js 가 v0.3.1~0.3.5 가 아니다 — 발사 중단"; exit 1; }
L=flat.log; : > "$L"
J=flat_jobs.txt; : > "$J"
python3 - > "$J" <<'PYEOF'
import json
picked = []
for s in range(1, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == 20:
        break
for s in picked:
    print(s)
PYEOF
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 20 -n 1 sh -c 'node dms.js --cond mat --wt rascld --seed "$0" --arms all --out dms_flat' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
