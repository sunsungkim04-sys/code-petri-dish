#!/usr/bin/env bash
# 해부 25 — 둘째 밀도의 배경 성장. 사전등록 해부25-사전등록-2026-10-01.md §2
#   밀도 8 배경(main/mat_d8_mu0p01_s*.json)을 만든 명령과 **같은 규칙 · 같은 틱 수**, 밀도만 바꾼다:
#     main_launch.sh :  node replay.js --cond mat --density 8 --seed S --ticks 50000 --mu 0.01 --out main
#     여기          :  node replay.js --cond mat --density D --seed S --ticks 50000 --mu 0.01 --out <ROOT>
#   기록 이름은 replay.js 규칙 그대로 mat_dD_mu0p01_sSSSSS.json — dms.js / invbud.js 가 --main-dir <ROOT> 로
#   같은 배경을 다시 키워 체크섬을 대조한다(inv25_launch.sh · G0).
#   자격(끝 우세 racld · 멸종 없음)은 여기서 판단하지 않는다 — pairs25_analyze.py qual 이 Δ 없이 센다.
# 본 발사:  PETRI_DENS=<사전등록 §2 에서 정한 밀도> PETRI_ROOT=main25 setsid nohup ./grow25_launch.sh > grow25.out 2>&1 < /dev/null &
# 파일럿:   PETRI_DENS="12 16" PETRI_SEEDS="$(seq 2591 2602)" PETRI_ROOT=pilot25g ./grow25_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
node --version | grep -q '^v20\.' || { echo "node 가 v20 이 아니다($(node --version)) — nvm 확인 · 발사 중단"; exit 1; }
# 계측기는 새 패치 없이 그대로(사전등록 §1) — 해시로 묶는다
echo "b459bb4a9933a38bbc3e06e3681aeb977dfa0a7cf680dd4d375d744415cf06b2  sim.js" | sha256sum -c --quiet - || { echo "sim.js 해시가 dms20/21 동결과 다르다 — 발사 중단"; exit 1; }
echo "991c9b9341f5b2886fad35c1cf8e717bdac7e7aaff99cc5b7fee494b4a2a4e2c  replay.js" | sha256sum -c --quiet - || { echo "replay.js 해시가 10-01 사본과 다르다 — 발사 중단"; exit 1; }
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
# 무인 재부팅 확인(lab101 메모) — 다르면 보고만 하고 멈춘다
UP=$(uptime -s)
[ "$UP" = "${PETRI_UPTIME:-2026-09-26 23:34:08}" ] || { echo "uptime -s $UP ≠ 2026-09-26 23:34:08 — 무인 재부팅 가능성 · 발사 중단(보고)"; exit 1; }
# 🔒 사전등록 §2a 의 밀도 규칙으로 파일럿 뒤 정한 밀도(검토 개정 R8 — inv25_launch.sh 의 DENS_FIXED 와 같은 손잡이 닫기)
DENS_FIXED=16
DENS=${PETRI_DENS:?PETRI_DENS 를 준다(본 발사는 사전등록 §2 에서 정한 밀도 하나)}
SEEDS=${PETRI_SEEDS:-$(seq 2501 2550)}
ROOT=${PETRI_ROOT:-main25}
PAR=${PETRI_PAR:-24}
[ "$PAR" -le 24 ] || { echo "병렬 $PAR > 코어 예산 24 — 발사 중단"; exit 1; }
for d in $DENS; do case $d in 12|16) ;; *) echo "밀도 $d — 12 또는 16 만 · 발사 중단"; exit 1;; esac; done
# 본 발사(main25)는 밀도 하나 · 시드 2501~2550 그대로여야 한다(판정선을 움직이는 손잡이를 닫는다)
if [ "$ROOT" = main25 ]; then
  [ "$(echo $DENS | wc -w)" = 1 ] || { echo "본 발사는 밀도 하나 — 발사 중단"; exit 1; }
  [ "$DENS" = "$DENS_FIXED" ] || { echo "본 발사 밀도 $DENS ≠ 고정 $DENS_FIXED — 발사 중단"; exit 1; }
  [ -z "${PETRI_SEEDS:-}" ] || { echo "본 발사에서 PETRI_SEEDS 를 바꿀 수 없다 — 발사 중단"; exit 1; }
fi
# 파일럿 시드와 본 시드는 겹치지 않는다
for s in $SEEDS; do
  if [ "$ROOT" = main25 ] && [ "$s" -ge 2591 ] && [ "$s" -le 2602 ]; then echo "본 시드 $s 가 파일럿 범위 — 발사 중단"; exit 1; fi
  for d in $DENS; do
    f="$ROOT/mat_d${d}_mu0p01_s$(printf %05d "$s").json"
    [ -e "$f" ] && { echo "$f 가 이미 있다(덮어쓰지 않는다) — 발사 중단"; exit 1; }
  done
done
# 코어 예산(검토 개정 R7): 다른 node(해부 24 등) + 이 발사의 병렬 ≤ 36. 넘으면 1분마다 다시 보고, PETRI_BUDGET_WAIT 초(기본 3600) 뒤에도 넘으면 멈춘다.
#   판정과 무관한 손잡이(기다리는 시간)만 있다 — 예산 36 은 고정.
BW=${PETRI_BUDGET_WAIT:-3600}; waited=0
while :; do
  OTHER=$(pgrep -c -x node || true); OTHER=${OTHER:-0}
  [ $(( OTHER + PAR )) -le 36 ] && break
  [ "$waited" -ge "$BW" ] && { echo "다른 node $OTHER + 병렬 $PAR > 36 — ${waited}초 기다려도 넘는다 · 발사 중단"; exit 1; }
  echo "다른 node $OTHER + 병렬 $PAR > 36 — 60초 뒤 다시 본다($(date '+%F %T'))"; sleep 60; waited=$(( waited + 60 ))
done
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for d in $DENS; do for s in $SEEDS; do echo "$d $s" >> "$J"; done; done
echo "해부 25 배경 성장 · 밀도 $DENS · 시드 $(echo $SEEDS | tr '\n' ' ')" >> "$L"
echo "작업 $(wc -l < "$J") · 병렬 $PAR · 다른 node $OTHER(합 $(( OTHER + PAR )) ≤ 36 · 기다림 ${waited}초) · $(date '+%F %T') · uptime -s $UP · node $(node --version)" >> "$L"
t0=$(date +%s)
R="$ROOT" xargs -P "$PAR" -n 2 sh -c 'node replay.js --cond mat --density "$0" --seed "$1" --ticks 50000 --mu 0.01 --out "$R" || echo "ERROR 작업 $0 $1"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
