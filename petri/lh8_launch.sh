#!/usr/bin/env bash
# 해부 8 B — 침입 Δ 가 쌓이는 경로(나이별 인구학). 사전등록 해부8-사전등록-2026-09-17.md §3
#   본실험 mat 접시 중 끝 우세 racld 인 접시의 21~40 번째 20개 × 팔 ref + neutral 10 + rascld × μ 0 × 6,000틱 → 20 작업
#   판정은 lifehist_analyze.py.
# 마른 실행: PETRI_BGS="24" PETRI_ROOT=smoke8b PETRI_TICKS=500 ./lh8_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "LIFEHIST_VERSION = '1.0.0'" lifehist.js || { echo "lifehist.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
ROOT=${PETRI_ROOT:-lh8}
TICKS=${PETRI_TICKS:-6000}
if [ -n "${PETRI_BGS:-}" ]; then BGS=$PETRI_BGS; else
BGS=$(python3 - <<'PYEOF'
import json
picked = []
for s in range(1, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
print(" ".join(map(str, picked[20:40])))
PYEOF
)
fi
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for s in $BGS; do echo "$s $TICKS $ROOT" >> "$J"; done
[ "$(awk 'NF!=3' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 3 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 20 -n 3 sh -c 'node lifehist.js --seed "$0" --ticks "$1" --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
