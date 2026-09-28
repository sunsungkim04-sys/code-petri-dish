#!/usr/bin/env bash
# 해부 18 — 순수 racld 공동체(성장 μ 0) × 밀도 {8 32} × 짝. 사전등록 해부18-사전등록-2026-09-28.md §4b(설계 변경)
#   R  invbud.js --ancestor racld --grow-mu 0 --density D --codes <n 코드> --mu 0                     × 밑 코드 6 × 시드 3~22 × 밀도 2 → 240
#   D  dms.js --ancestor racld --grow-mu 0 --density D --wt B --arms list --codes <n 코드> --mu-assay {0 0.003} × 6 × 20 × 2 × 2 → 480
#   회귀: 두 계측기의 grow_checksum 이 같고 dms mut 팔 [t, n, m] = invbud 팔. 판정: pairs18_analyze.py
# 마른 실행(파일럿): PETRI_NBG=2 PETRI_BG_FROM=23 PETRI_NBASES="racld rascled" PETRI_ROOT=pilot18c ./inv18_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "densityV" dms.js || { echo "dms.js 에 --density 가 없다 — 발사 중단"; exit 1; }
grep -q "densityV" invbud.js || { echo "invbud.js 에 --density 가 없다 — 발사 중단"; exit 1; }
grep -q "growMu" dms.js && grep -q "ancestorV" invbud.js || { echo "--grow-mu / --ancestor 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
NBASES=${PETRI_NBASES:-"racld rascld rsacld racldx acld rascled"}
DENSES=${PETRI_DENSES:-"8 32"}
ANC=${PETRI_ANC:-racld}
GMU=${PETRI_GMU:-0}
MUS=${PETRI_MUS:-"0 0.003"}
PAR=${PETRI_PAR:-60}
ROOT=${PETRI_ROOT:-inv18}
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
mt() { echo "mu$(echo "$1" | sed 's/\./p/')"; }
for b in $NBASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for d in $DENSES; do for s in $BGS; do
    echo "IB 0 $s $b $codes $d ${ROOT}ib_$b/d$d/mu0" >> "$J"
    for mu in $MUS; do echo "D $mu $s $b $codes $d ${ROOT}D_$b/d$d/$(mt $mu)" >> "$J"; done
  done; done
done
[ "$(awk 'NF!=7' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 7 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS · 밀도 $DENSES · 시조 $ANC · 성장 μ $GMU · μ $MUS · 작업 $(wc -l < "$J") · $(date '+%F %T') · uptime -s $(uptime -s)" >> "$L"
t0=$(date +%s)
ANC=$ANC GMU=$GMU xargs -P "$PAR" -L 1 sh -c '
  if [ "$0" = IB ]; then node invbud.js --ancestor "$ANC" --grow-mu "$GMU" --density "$5" --seed "$2" --mu "$1" --codes "$4" --out "$6";
  else node dms.js --cond mat --ancestor "$ANC" --grow-mu "$GMU" --density "$5" --wt "$3" --seed "$2" --mu-assay "$1" --exact 1 --arms list --codes "$4" --out "$6"; fi' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
