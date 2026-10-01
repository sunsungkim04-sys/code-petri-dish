#!/usr/bin/env bash
# 해부 22 — 축소 모형 D(잘 섞인 세계 · 공간 없음) 씨앗 2 → 20. 사전등록 해부22-사전등록-2026-10-01.md
#   같은 코드(reduced_d.py · sha256 d9473aa0… = 09-27 원래 판)와 같은 설정으로 시드마다 1 작업:
#   reduced_d.py --cells 1000 --density 8 --grow 6000 --ticks 3000 --seeds 1 --seed0 <s>
#                --codes-of racld,rascld --mus 0,0.003 --rules orig,rem --resident racld --frac 0.1
#   판정 시드 2201~2220(20 작업). 판정: red22_analyze.py · 원래 판 재현 확인: red22_merge.py + reduced_d_analyze.py
# 마른 실행(원래 판 재현 · 시드 101 102): PETRI_SEEDS="101 102" PETRI_ROOT=pilot22 PETRI_PAR=2 ./red22_launch.sh
# 본 발사:                               PETRI_ROOT=red22 PETRI_PAR=8 setsid nohup ./red22_launch.sh > /dev/null 2>&1 &
# 순수 python(단일 스레드) — 작업 하나 = 코어 하나. 코어 예산 8(해부 22 몫).
set -u
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
PY=${PETRI_PY:-/usr/bin/python3}
IMPL_SHA=d9473aa0c1faaecaeb36df11d6cb3e6a62194f86e65c34ab0ec5a258da70ea43
got=$(sha256sum reduced_d.py | awk '{print $1}')
[ "$got" = "$IMPL_SHA" ] || { echo "reduced_d.py sha256 $got ≠ 원래 판 $IMPL_SHA — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 2201 2220)}
PAR=${PETRI_PAR:-8}
[ "$PAR" -le 8 ] || { echo "PETRI_PAR=$PAR — 코어 예산 8 초과 · 발사 중단"; exit 1; }
ROOT=${PETRI_ROOT:-red22}
[ -e "$ROOT" ] && { echo "$ROOT 가 이미 있다 — 덮어쓰지 않는다 · 발사 중단"; exit 1; }
mkdir -p "$ROOT"
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for s in $SEEDS; do echo "$s $ROOT/s$s.json" >> "$J"; done
[ "$(awk 'NF!=2' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 2 가 아니다 — 발사 중단"; exit 1; }
echo "$got  reduced_d.py" > "$ROOT/impl.sha256"
echo "해부 22 · $($PY --version 2>&1) · 시드 $(echo $SEEDS | tr '\n' ' ')· 작업 $(wc -l < "$J") · 병렬 $PAR · $(date '+%F %T') · uptime -s $(uptime -s) · impl $got" >> "$L"
t0=$(date +%s)
PY="$PY" xargs -P "$PAR" -L 1 sh -c '
  ts=$(date +%s)
  "$PY" reduced_d.py --cells 1000 --density 8 --grow 6000 --ticks 3000 --seeds 1 --seed0 "$0" \
      --codes-of racld,rascld --mus 0,0.003 --rules orig,rem --resident racld --frac 0.1 --out "$1" > "$1.log" 2>&1
  rc=$?
  echo "시드 $0 · rc $rc · $(( $(date +%s) - ts ))초 · $(date +%T)"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
