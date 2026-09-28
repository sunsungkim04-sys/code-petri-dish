#!/usr/bin/env bash
# 해부 8 파일럿 요약 재생성기 — 서버 ~/petri/stage8 에서 실행 (판정 아님 · 노트 §5)
#   A: 재료 먼저 시드 9101~9110 끝 상위 코드 + 판정기 마른 실행(9101~9103, 원래 규칙 짝 pilotA_orig)
#   B: 되감기 2 배경 101~103 판정기 마른 실행
#   C: 시드 9201~9203 · racld · rascld 판정기 마른 실행
cd "$HOME/petri/stage8" || exit 1
# 판정기는 자기 폴더에 _RESULT_*.json 을 쓴다 — 서버 루트(본 판정 자리)를 더럽히지 않게 stage8 에 복사해 돌린다
cp ../rewind8_analyze.py ../lifehist_analyze.py ../spec_analyze.py .
echo "== 🟡 해부 8 파일럿 요약 — 판정 아님"
echo; echo "== A 재료 먼저 시드 9101~9110 끝 상위"
python3 - <<'PYEOF'
import json, glob
for f in sorted(glob.glob("pilotA/*.json")):
    z = json.load(open(f))
    print("  ", f.split("/")[-1], "멸종", z["extinct_at"], "끝 상위", z["final_counts"][:4])
PYEOF
echo; PETRI_ROOT=pilotA PETRI_SEEDS=9101-9103 PETRI_ORIG=pilotA_orig PETRI_ORIG_V=0.3.3 python3 rewind8_analyze.py "$HOME/petri/stage8" | sed '$d'
echo; PETRI_ROOT=pilotB PETRI_BGS=101,102,103 PETRI_MAIN=rewind2 python3 lifehist_analyze.py "$HOME/petri/stage8" | sed '$d'
echo; PETRI_ROOT=pilotC PETRI_SEEDS=9201,9202,9203 PETRI_CODES=racld,rascld python3 spec_analyze.py "$HOME/petri/stage8" | sed '$d'
