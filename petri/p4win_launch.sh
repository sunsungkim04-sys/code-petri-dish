#!/bin/sh
# p4win_launch.sh — 원고 Table S10 · S11 P4 의 같은 창 재계수 발사기 (사후 · 판정 아님 · 2026-10-07)
#   본 실행: P4W_ROOT=p4w P4W_PY=$HOME/miniforge3/bin/python3 ./p4win_launch.sh
#            → p4w/h26 (해부26 ρ8 26001~26020 · glob) · p4w/h27 (해부27 ρ8 27001~27020 · glob · loc)
#   시운전:  P4W_ROOT=p4w_pilot P4W_SEEDS26="26901" P4W_SEEDS27="27901" ./p4win_launch.sh   (파일럿 시드만)
#   분리 실행: … nohup ./p4win_launch.sh < /dev/null > p4w.nohup 2>&1 &   (stdin 을 닫지 않으면 ssh 가 안 끝난다)
# 같은 폴더에 mini27.py(해부27 동결본) · galA_window.py(갈래 A 동결본) · p4win_job.sh 가 있어야 한다. 둘 다 고치지 않고 쓴다.
# 🔑 작업 줄은 3 칸(연구 ρ 시드) → xargs -n 3.
# 🔑 판 고정 — python 3.13.12 · numpy 2.4.3 (해부26 · 27 원 출력 · 해부27 G2a 와 같은 판).
# 🔑 끝 줄의 곁 파일 수로 완료를 센다 — mini27 의 'qualified=' 줄은 래퍼가 곁 파일을 쓰기 전에 찍힌다.
set -eu
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
P4W_PY=${P4W_PY:-python3}
P4W_ROOT=${P4W_ROOT:?P4W_ROOT 필요}
PAR=${P4W_PAR:-8}
export P4W_PY P4W_ROOT
if [ "$PAR" -gt 16 ]; then echo "병렬 $PAR > 16 — 거부(lab101 공유 중)"; exit 1; fi
for f in mini27.py galA_window.py p4win_job.sh; do [ -f "$f" ] || { echo "$f 가 없다 — 올릴 것"; exit 1; }; done

NEED_NP=${P4W_NUMPY:-2.4.3}
GOT=$($P4W_PY -c 'import numpy;print(numpy.__version__)' 2>/dev/null || echo none)
if [ "$GOT" != "$NEED_NP" ]; then
  echo "numpy $GOT ≠ $NEED_NP ($P4W_PY) — 거부. lab101: P4W_PY=\$HOME/miniforge3/bin/python3"
  exit 1
fi

J="$P4W_ROOT.jobs"
L="$P4W_ROOT.log"
JN="$J.new.$$"
: > "$JN"
trap 'rm -f "$JN"' EXIT
case "$P4W_ROOT" in
  *pilot*)
    for s in ${P4W_SEEDS26:-}; do
      [ "$s" -ge 26901 ] && [ "$s" -le 26910 ] || { echo "파일럿 26 시드는 26901~26910 만: $s"; exit 1; }
      echo "26 8 $s" >> "$JN"
    done
    for s in ${P4W_SEEDS27:-}; do
      [ "$s" -ge 27901 ] && [ "$s" -le 27910 ] || { echo "파일럿 27 시드는 27901~27910 만: $s"; exit 1; }
      echo "27 8 $s" >> "$JN"
    done ;;
  *)
    [ -z "${P4W_SEEDS26:-}${P4W_SEEDS27:-}" ] || { echo "본 실행은 시드를 고정한다 — P4W_SEEDS* 를 쓰지 말 것"; exit 1; }
    for s in $(seq 26001 26020); do echo "26 8 $s" >> "$JN"; done
    for s in $(seq 27001 27020); do echo "27 8 $s" >> "$JN"; done ;;
esac
[ -s "$JN" ] || { echo "작업 없음"; exit 1; }

# ── 여기부터만 디스크의 옛 판을 건드린다 ──
mkdir -p "$P4W_ROOT/h26" "$P4W_ROOT/h27"
[ -f "$L" ] && mv "$L" "$L.prev-$(date '+%H%M%S')"
mv "$JN" "$J"
trap - EXIT
echo "===== 시작 $(date '+%F %T') · 작업 $(wc -l < "$J" | tr -d ' ') · 병렬 $PAR · $($P4W_PY -c 'import numpy,sys;print("python",sys.version.split()[0],"numpy",numpy.__version__)') · 부팅 $(uptime -s 2>/dev/null || echo ?)" > "$L"
{ sha256sum mini27.py galA_window.py p4win_job.sh p4win_launch.sh 2>/dev/null || shasum -a 256 mini27.py galA_window.py p4win_job.sh p4win_launch.sh; } >> "$L"
xargs -P "$PAR" -n 3 sh -c 'sh ./p4win_job.sh "$0" "$1" "$2" || echo "ERROR 작업 $0 $1 $2"' < "$J" >> "$L" 2>&1
NW=$(ls "$P4W_ROOT"/h26/*.win.json "$P4W_ROOT"/h27/*.win.json 2>/dev/null | wc -l | tr -d ' ')
echo "===== 끝 $(date '+%F %T') · 작업 $(wc -l < "$J" | tr -d ' ') · 곁 파일 $NW · ERROR $(grep -c '^ERROR' "$L") · 부팅 $(uptime -s 2>/dev/null || echo ?) · P4W_END" >> "$L"
