#!/usr/bin/env bash
# 해부 28 — 다시 굴림 확률 p 의 용량–반응. 사전등록 사전등록/해부28-사전등록-2026-10-02.md
#   계측기 · 팔 · 배경 · μ 사다리는 해부 9 와 같다(inv9_launch.sh):
#     node dms.js --cond mat --wt racld --seed S --arms list --codes rascld,rsacld,racldx --mu-assay μ --exact 1 --redraw-p P
#     배경 = main/ 1~100 중 '끝 우세 racld · 멸종 없음' 앞 20 · μ 0 · 0.1 · 0.3 · 0.5 · 1%
#   새로 돌리는 것: p = 0.25 · 0.5 · 0.75 (300 작업). 양 끝은 기존 기록을 쓴다 — p = 1 ≡ inv6/(원래 규칙) · p = 0 ≡ inv9mem/(기억 주사위).
#   회귀(${ROOT}reg/): 첫 배경 · μ 0.3% · 1% 에서
#     p1   = --redraw-p 1  → inv6/ 팔 궤적과 같아야 한다
#     p0   = --redraw-p 0  → inv9mem/ 팔 궤적과 같아야 한다
#     none = 깃발 없음      → inv6/ 과 같아야 한다(sim 0.3.6 패치가 기본 궤적을 안 바꿨다)
#     det_a · det_b = --redraw-p 0.5 두 번 → 서로 같아야 한다(동전 흐름 결정론)
# 마른 실행: PETRI_NBG=1 PETRI_PS=0.5 PETRI_MUS="0 0.01" PETRI_ROOT=dry28 ./inv28_launch.sh
# 본 발사:   PETRI_ROOT=inv28 setsid nohup ./inv28_launch.sh > inv28.nohup 2>&1 &
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
node --version | grep -q '^v20\.' || { echo "node 가 v20 이 아니다($(node --version)) — 발사 중단"; exit 1; }
grep -q "VERSION = '0.3.6'" sim.js || { echo "sim.js 가 v0.3.6 이 아니다 — 발사 중단"; exit 1; }
grep -q "redraw-p" dms.js || { echo "dms.js 에 --redraw-p 가 없다 — 발사 중단"; exit 1; }
[ -d inv6/mu0p003 ] && [ -d inv9mem/mu0p003 ] || { echo "회귀 대조 기록(inv6 · inv9mem)이 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.001 0.003 0.005 0.01"}
PS=${PETRI_PS:-"0.25 0.5 0.75"}
PAR=${PETRI_PAR:-20}
ROOT=${PETRI_ROOT:-inv28}
[ "$PAR" -le 24 ] || { echo "병렬 $PAR > 코어 예산 24(공용 서버 · 해부 27 동시 실행) — 발사 중단"; exit 1; }
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
BGS=$(NBG=$NBG python3 - <<'PYEOF'
import json, os
nbg = int(os.environ["NBG"]); picked = []
for s in range(1, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == nbg:
        break
print(" ".join(map(str, picked)))
PYEOF
)
S0=$(echo $BGS | awk '{print $1}')
mt() { echo "mu$(echo "$1" | sed 's/\./p/')"; }
for mu in 0.003 0.01; do
  echo "$mu $S0 1 ${ROOT}reg/p1/$(mt $mu)" >> "$J"
  echo "$mu $S0 0 ${ROOT}reg/p0/$(mt $mu)" >> "$J"
  echo "$mu $S0 none ${ROOT}reg/none/$(mt $mu)" >> "$J"
done
echo "0.003 $S0 0.5 ${ROOT}reg/det_a/mu0p003" >> "$J"
echo "0.003 $S0 0.5 ${ROOT}reg/det_b/mu0p003" >> "$J"
for p in $PS; do
  for mu in $MUS; do
    for s in $BGS; do echo "$mu $s $p ${ROOT}/p$(echo "$p" | sed 's/\./p/')/$(mt $mu)" >> "$J"; done
  done
done
[ "$(awk 'NF!=4' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") · 병렬 $PAR · $(date '+%F %T') · uptime -s $(uptime -s) · node $(node --version)" >> "$L"
t0=$(date +%s)
xargs -P "$PAR" -n 4 sh -c 'if [ "$2" = none ]; then F=""; else F="--redraw-p $2"; fi; node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 $F --out "$3" || echo "ERROR 작업 $0 $1 $2 $3"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
