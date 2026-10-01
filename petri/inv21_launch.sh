#!/usr/bin/env bash
# 해부 21 — §3.2 교차(씨앗 침입 Δ 의 부호가 바뀌는 μ, 해부 6C 0.18%)를 새 배경에서 재현. 사전등록 해부21-사전등록-2026-10-01.md
#   계측기 · 팔 · μ 사다리는 해부 6C(inv_launch.sh)와 같다:
#     node dms.js --cond mat --wt racld --seed S --arms list --codes rascld,rsacld,racldx --mu-assay μ --exact 1
#     μ 0 · 0.1 · 0.2 · 0.3 · 0.4 · 0.5 · 0.75 · 1%
#   바뀐 것은 배경뿐: 해부 6C 는 본실험 시드 3~22(main/ 기록이 있어 체크섬 대조). 여기는 **한 번도 안 쓴 시드 4201~4250**
#   (main/ 기록이 없다 → checksum 'no-main' 이 정상). 배경은 dms.js 가 본실험과 같은 설정(mat · 밀도 8 · μ 0.01 · 50,000틱)으로 키운다.
#   후보 50 전부를 8 μ 로 돌리고, '끝 우세 racld · 멸종 없음' 인 배경을 시드 순으로 앞 40 개까지 판정에 쓴다(inv21_analyze.py).
#   회귀(R): 해부 6C 배경 시드 3 을 μ 0.2% · 0.3% 로 다시 돌려 inv6/ 기록과 팔 궤적 전부가 같아야 한다(같은 계측기 확인).
# 마른 실행: PETRI_CANDS="4291 4292" PETRI_ROOT=dry21 PETRI_PAR=18 ./inv21_launch.sh
# 본 발사:   PETRI_ROOT=inv21 setsid nohup ./inv21_launch.sh > inv21.nohup 2>&1 &
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
node --version | grep -q '^v20\.' || { echo "node 가 v20 이 아니다($(node --version)) — nvm 을 확인. 발사 중단"; exit 1; }
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "const exact = args.exact === '1';" dms.js || { echo "dms.js 에 --exact 가 없다 — 발사 중단"; exit 1; }
[ -d inv6/mu0p003 ] && [ -d inv6/mu0p002 ] || { echo "회귀 대조 inv6/ 기록이 없다 — 발사 중단"; exit 1; }
CANDS=${PETRI_CANDS:-$(seq 4201 4250)}
MUS=${PETRI_MUS:-0 0.001 0.002 0.003 0.004 0.005 0.0075 0.01}
REG_MUS=${PETRI_REG_MUS:-0.002 0.003}
PAR=${PETRI_PAR:-28}
ROOT=${PETRI_ROOT:-inv21}
[ "$PAR" -le 28 ] || { echo "병렬 $PAR > 코어 예산 28 — 발사 중단"; exit 1; }
for s in $CANDS; do [ -e "main/mat_d8_mu0p01_s$(printf %05d "$s").json" ] && { echo "시드 $s 에 본실험 기록이 있다 — 새 배경이 아니다. 발사 중단"; exit 1; }; done
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
mt() { echo "mu$(echo "$1" | sed 's/\./p/')"; }
for mu in $REG_MUS; do echo "$mu 3 ${ROOT}reg/$(mt "$mu")" >> "$J"; done
for mu in $MUS; do for s in $CANDS; do echo "$mu $s $ROOT/$(mt "$mu")" >> "$J"; done; done
echo "후보 배경 $(echo $CANDS | tr '\n' ' ')" >> "$L"
echo "작업 $(wc -l < "$J") · 병렬 $PAR · $(date '+%F %T') · uptime -s $(uptime -s) · node $(node --version)" >> "$L"
t0=$(date +%s)
xargs -P "$PAR" -n 3 sh -c 'node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 --out "$2" || echo "ERROR 작업 $0 $1 $2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
