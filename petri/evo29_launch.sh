#!/usr/bin/env bash
# 해부 29 — 진화 궤적을 한 번 굴림 규칙에서 (M3). 사전등록 해부29-사전등록-2026-10-03.md §1~§2
#   해부 17(evo17/ · draw-first)과 같은 세계 · 같은 시드 1~20 · 5만 틱 · --loops 1 에서 복사 규칙만 바꾼다.
#     do 팔 : replay.js --remember-die 1 (기억하는 주사위 · draw-once) · μ {0.001 0.003 0.005 0.01 0.03 0.05} × 시드 → 120
#     ff 팔 : replay.js --find-first 1  (재료 먼저 · find-first)      · μ {0.001 0.003 0.005 0.01}           × 시드 →  80
#   draw-first 쪽은 새로 돌리지 않는다 — evo17/ 기록을 쓴다(G0a 가 sim 0.3.6 에서도 같은 궤적임을 보인다).
#   회귀(G0 · 판정기 evo29_analyze.py 가 대조) — 본 작업보다 먼저 같은 xargs 에 넣는다:
#     G0a 깃발 없음 · --loops 1 · μ 0.01 시드 1 · 2(초기 멸종) · μ 0.001 시드 1 → ${ROOT}reg/df/ ≡ evo17/ (기록 전부)
#     G0b --find-first 1 · μ 0.01 시드 101 · 102 → ${ROOT}reg/ff/ ≡ rewind8ff/ (sim 0.3.3 기록 · 기록 전부)
#     G0c --remember-die 1 · μ 0.01 시드 901 · 5,000 틱 두 번 → ${ROOT}reg/do_a/ ≡ ${ROOT}reg/do_b/ (결정론)
#     G0d dms.js 식 켜기(makeWorld 뒤 w.o.rememberDie = true) · 같은 칸 · 5,000 틱 → 체크섬을 ${ROOT}reg/dms_style.txt 에 → G0c 와 같아야 한다
# 마른 실행(파일럿 · 판정 밖 시드):
#   PETRI_SEEDS="901 902" PETRI_MUS_DO="0.001 0.01 0.05" PETRI_MUS_FF="0.001 0.01" PETRI_ROOT=pilot29 ./evo29_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.6'" sim.js || { echo "sim.js 가 v0.3.6 이 아니다 — 발사 중단"; exit 1; }
grep -q "loopsOn" replay.js || { echo "replay.js 에 --loops 가 없다 — 발사 중단"; exit 1; }
grep -q "args\['remember-die'\]" replay.js || { echo "replay.js 에 --remember-die 가 없다 — 발사 중단"; exit 1; }
[ -d evo17 ] && [ "$(ls evo17/*.json | wc -l)" = 80 ] || { echo "evo17/ 기록 80 개가 없다 — 발사 중단"; exit 1; }
[ -d rewind8ff ] || { echo "rewind8ff/ 가 없다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 1 20)}
MUS_DO=${PETRI_MUS_DO:-"0.001 0.003 0.005 0.01 0.03 0.05"}
MUS_FF=${PETRI_MUS_FF:-"0.001 0.003 0.005 0.01"}
PAR=${PETRI_PAR:-16}
[ "$PAR" -le 16 ] || { echo "병렬 상한 16 (공용 서버) — 발사 중단"; exit 1; }
ROOT=${PETRI_ROOT:-evo29}
REG=${ROOT}reg
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
# 줄 = 깃발 μ 시드 틱 출력폴더 (깃발: df · ff · do)
for s in 1 2; do echo "df 0.01 $s 50000 $REG/df" >> "$J"; done
echo "df 0.001 1 50000 $REG/df" >> "$J"
for s in 101 102; do echo "ff 0.01 $s 50000 $REG/ff" >> "$J"; done
echo "do 0.01 901 5000 $REG/do_a" >> "$J"
echo "do 0.01 901 5000 $REG/do_b" >> "$J"
for mu in $MUS_DO; do for s in $SEEDS; do echo "do $mu $s 50000 $ROOT/do" >> "$J"; done; done
for mu in $MUS_FF; do for s in $SEEDS; do echo "ff $mu $s 50000 $ROOT/ff" >> "$J"; done; done
[ "$(awk 'NF!=5' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 5 가 아니다 — 발사 중단"; exit 1; }
echo "시드 $SEEDS · do μ $MUS_DO · ff μ $MUS_FF · 작업 $(wc -l < "$J") · 병렬 $PAR · $(date '+%F %T') · uptime -s $(uptime -s) · load $(cut -d' ' -f1 /proc/loadavg)" | tr '\n' ' ' >> "$L"; echo >> "$L"
mkdir -p "$REG"
node -e '
const Sim = require("./sim.js");
const w = Sim.makeWorld({ seed: 901, density: 8, mu: 0.01, light: 1, mat: true, energy: false });
w.o.rememberDie = true;   // dms.js 202 행과 같은 켜기
while (w.tick < 5000 && w.extinctAt < 0) Sim.tick(w);
console.log(Sim.VERSION, w.tick, Sim.checksum(w));
' > "$REG/dms_style.txt" 2>&1
t0=$(date +%s)
xargs -P "$PAR" -n 5 sh -c '
case "$0" in df) F="";; ff) F="--find-first 1";; do) F="--remember-die 1";; *) echo "깃발 $0?"; exit 1;; esac
node replay.js --cond mat --density 8 --mu "$1" --seed "$2" --ticks "$3" --loops 1 $F --out "$4"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
