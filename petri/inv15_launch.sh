#!/usr/bin/env bash
# 해부 15 — 짝 안의 μ 사다리(B) · 세 규칙 × 짝(C) · 새 배경 재현(E). 사전등록 해부15-사전등록-2026-09-27.md §1~§2
#   B   dms.js --wt B --arms list --codes <n 코드> --mu-assay μ   μ {0 0.001 0.003 0.005 0.01} · 밑 코드 6 · 배경 3~22           → 6 × 5 × 20 = 600
#       (μ 0 여섯 · 둘째 가족 0.003 은 inv12/inv13 이웃 기록과 겹친다 — 목록 팔 = 이웃 팔 회귀 검사)
#   C   같은 명령 + --find-first 1 | --remember-die 1 · μ {0 0.003} · 배경 3~22                                      → 6 × 2 × 2 × 20 = 480
#   E   같은 명령(원래 규칙 · μ 0) + invbud.js --mu 0 · 새 배경(26 부터 racld 정상 20 = 26~50)                     → 6 × 20 × 2 = 240
#   계측기는 해부 12 · 14 것 그대로 — 새 패치 없음. 판정: pairs15_analyze.py
# 마른 실행(파일럿): PETRI_NBG=2 PETRI_BG_FROM=23 PETRI_EBG_FROM=51 PETRI_NBASES="racld rascled" PETRI_ROOT=pilot15 ./inv15_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "INVBUD_VERSION = '1.0.0'" invbud.js || { echo "invbud.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
grep -q "rememberDie\|remember-die" dms.js || { echo "dms.js 에 --remember-die 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
NBASES=${PETRI_NBASES:-"racld rascld rsacld racldx acld rascled"}
PARTS=${PETRI_PARTS:-"B C E"}
MUS_B=${PETRI_MUS_B:-"0 0.001 0.003 0.005 0.01"}
MUS_C=${PETRI_MUS_C:-"0 0.003"}
ROOT=${PETRI_ROOT:-inv15}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
pick_bgs() {  # $1 = 첫 시드 · $2 = 개수
  FIRST=$1 NBG=$2 python3 - <<'PYEOF'
import json, os
nbg = int(os.environ["NBG"]); first = int(os.environ["FIRST"]); picked = []
for s in range(first, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == nbg:
        break
print(" ".join(map(str, picked)))
PYEOF
}
BGS=$(pick_bgs "${PETRI_BG_FROM:-1}" "$NBG")
EBGS=$(pick_bgs "${PETRI_EBG_FROM:-26}" "$NBG")
CODEMAP=$(NBASES="$NBASES" python3 - <<'PYEOF'
import os
def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None
def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None: out.setdefault(k, L)
    return out
for b in os.environ["NBASES"].split():
    print(b + ":" + ",".join(inserts(b, "n")))
PYEOF
)
mt() { echo "mu$(echo "$1" | sed 's/\./p/')"; }
for b in $NBASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for p in $PARTS; do case $p in
    B) for mu in $MUS_B; do for s in $BGS; do echo "D $mu $s $b $codes orig ${ROOT}B_$b/$(mt $mu)" >> "$J"; done; done;;
    C) for rule in ff rem; do for mu in $MUS_C; do for s in $BGS; do echo "D $mu $s $b $codes $rule ${ROOT}C_$b/$rule/$(mt $mu)" >> "$J"; done; done; done;;
    E) for s in $EBGS; do echo "D 0 $s $b $codes orig ${ROOT}E_$b/mu0" >> "$J"; echo "IB 0 $s $b $codes orig ${ROOT}Eib_$b/mu0" >> "$J"; done;;
  esac; done
done
[ "$(awk 'NF!=7' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 7 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS · 새 배경 $EBGS" >> "$L"
echo "코드 목록: $CODEMAP" | tr '\n' ' ' >> "$L"; echo >> "$L"
echo "작업 $(wc -l < "$J") (B $(grep -c "${ROOT}B_" "$J") · C $(grep -c "${ROOT}C_" "$J") · E $(grep -c "${ROOT}E_\|${ROOT}Eib_" "$J")) · $(date '+%F %T') · uptime -s $(uptime -s)" >> "$L"
t0=$(date +%s)
xargs -P 60 -L 1 sh -c '
  case "$5" in ff) R="--find-first 1";; rem) R="--remember-die 1";; *) R="";; esac
  if [ "$0" = IB ]; then node invbud.js --seed "$2" --mu "$1" --codes "$4" --out "$6";
  else node dms.js --cond mat --wt "$3" --seed "$2" --mu-assay "$1" --exact 1 --arms list --codes "$4" $R --out "$6"; fi' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
