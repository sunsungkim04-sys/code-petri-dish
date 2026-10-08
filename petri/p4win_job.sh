#!/bin/sh
# p4win_job.sh — P4 같은 창 재계수 작업 하나 (p4win_launch.sh 가 xargs -n 3 으로 부른다 · 작업 줄 = 연구 ρ 시드)
#   연구 26 : 해부26 배경 — mini27 glob 경로(= mini26 · 해부27 G2a 바이트 동일) · roll · μ 0 · 기준 팔 0
#   연구 27 : 해부27 배경 — glob · loc · roll · μ 0 · 기준 팔 0
# mut 팔의 난수 흐름은 (시드, ρ, μ, 팔, 용도) 로만 키한다 — 기준 팔 · 다른 μ · 다른 규칙을 빼도 mut 팔은 원 실행과 같다
# (p4win_analyze.py 의 G1 이 원 출력과 정규형 JSON 으로 대조한다).
set -eu
S=$1; RHO=$2; SEED=$3
: "${P4W_PY:?}" "${P4W_ROOT:?}"
case "$S" in
  26) G=glob ;;
  27) G=glob,loc ;;
  *) echo "모르는 연구 $S"; exit 2 ;;
esac
exec nice -n 10 "$P4W_PY" galA_window.py --design std job --rho "$RHO" --seed "$SEED" --geoms "$G" --rules roll --mus 0 --nref 0 --out "$P4W_ROOT/h$S"
