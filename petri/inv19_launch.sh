#!/usr/bin/env bash
# 해부 19 — 꽉 찬 접시 가설의 전향 검정: 순수 racld 공동체(성장 μ 0) × 밀도 {8 32} × 수명 {150 600} × 짝 · μ 0 만. 사전등록 해부19-사전등록-2026-09-28.md §4
#   D  dms.js --ancestor racld --grow-mu 0 --density D --age0 A --age-var A --grow 50000·f --ticks 3000·f --every 250·f (f = A/300 · 해부 11 관례)
#      --wt B --arms list --codes <n 코드> --mu-assay 0 --exact 1                × 밑 코드 6 × 시드 3~22 × 밀도 2 × 수명 2 → 480
#   수명 300 칸은 해부 18 기록(inv18D_<b>/d<D>/mu0)을 그대로 쓴다 — 다시 돌리지 않는다. 판정: pairs19_analyze.py
# 파일럿: PETRI_NBG=2 PETRI_BG_FROM=23 PETRI_NBASES="racld rascled" PETRI_ROOT=pilot19 PETRI_PAR=16 ./inv19_launch.sh
# 회귀(명시한 수명 300 = 기본값): PETRI_AGES=300 PETRI_DENSES=8 PETRI_NBASES=racld PETRI_NBG=1 PETRI_BG_FROM=24 PETRI_ROOT=pilot19reg ./inv19_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "age0V" dms.js || { echo "dms.js 에 --age0 가 없다 — 발사 중단"; exit 1; }
grep -q "densityV" dms.js && grep -q "growMu" dms.js || { echo "dms.js 에 --density / --grow-mu 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
NBASES=${PETRI_NBASES:-"racld rascld rsacld racldx acld rascled"}
DENSES=${PETRI_DENSES:-"8 32"}
AGES=${PETRI_AGES:-"150 600"}
ANC=${PETRI_ANC:-racld}
GMU=${PETRI_GMU:-0}
PAR=${PETRI_PAR:-60}
ROOT=${PETRI_ROOT:-inv19}
for a in $AGES; do [ $(( (250 * a) % 300 )) = 0 ] || { echo "수명 $a: 250·f 가 정수가 아니다 — 발사 중단"; exit 1; }; done
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
BGS=$(NBG=$NBG python3 - <<'PYEOF'
import json, os
nbg = int(os.environ["NBG"]); first = int(os.environ.get("PETRI_BG_FROM", "1")); picked = []
for s in range(first, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == nbg:
        break
print(" ".join(map(str, picked)))
PYEOF
)
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
for b in $NBASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for d in $DENSES; do for a in $AGES; do for s in $BGS; do
    echo "D $a $s $b $codes $d ${ROOT}D_$b/d$d/a$a" >> "$J"
  done; done; done
done
[ "$(awk 'NF!=7' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 7 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS · 밀도 $DENSES · 수명 $AGES · 시조 $ANC · 성장 μ $GMU · μ 0 · 작업 $(wc -l < "$J") · $(date '+%F %T') · uptime -s $(uptime -s)" >> "$L"
t0=$(date +%s)
ANC=$ANC GMU=$GMU xargs -P "$PAR" -L 1 sh -c '
  G=$(( 50000 * $1 / 300 )); T=$(( 3000 * $1 / 300 )); E=$(( 250 * $1 / 300 ));
  node dms.js --cond mat --ancestor "$ANC" --grow-mu "$GMU" --density "$5" --age0 "$1" --age-var "$1" --grow "$G" --ticks "$T" --every "$E" --wt "$3" --seed "$2" --mu-assay 0 --exact 1 --arms list --codes "$4" --out "$6"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
