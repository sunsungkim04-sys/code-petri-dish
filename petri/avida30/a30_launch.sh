#!/bin/sh
# 해부 30 — confirmatory launch of the Avida conserved-material replication (E1 + E2).
# Refuses to run unless (1) every frozen hash verifies, (2) the Avida build is the patched branch
# with exactly the frozen patch applied, (3) the regression guards G0 pass.
# Usage (on lab101, from ~/projects/avida-port/e30):  sh a30_launch.sh  [e1|e2|all]
set -eu
cd "$(dirname "$0")"
mkdir -p e1 e2 g0
WHAT=${1:-all}
PAR=${A30_PAR:-16}
AV=$HOME/projects/avida-port/avida

# (1) frozen hashes: dms30_frozen.sha256 is a copy of the vault's petri/dms30_frozen.sha256 (paths avida30/...);
#     the prereg note line is checked in the vault, every other line here
sed 's#  avida30/#  #' dms30_frozen.sha256 | grep -v '^#' | grep -v '사전등록' | sha256sum -c

# (2) the build is the frozen patch on devosoft/avida 47f13da, branch material-draw, installed
[ "$(git -C $AV rev-parse --abbrev-ref HEAD)" = "material-draw" ]
[ "$(git -C $AV rev-parse --short=7 HEAD)" = "47f13da" ]
( cd $AV && git diff; for f in avida-core/source/main/cMaterial.h avida-core/source/main/cMaterial.cc; do diff -u --label /dev/null --label b/$f /dev/null $f || true; done ) > /tmp/a30_current.patch
cmp /tmp/a30_current.patch material-draw.patch
cmp $AV/cbuild/bin/avida $AV/cbuild/work/avida

# (3) G0 regression guards (pilot seeds only; results go to g0/)
A30_ROOT=$PWD/g0 python3 a30_e0.py c > g0_c.txt
grep -c "spop_identical True census_identical True" g0_c.txt | grep -qx 3

if [ "$WHAT" = e1 ] || [ "$WHAT" = all ]; then
  A30_ROOT=$PWD/e1 A30_PAR=$PAR A30_SEEDS="1 2 3 4 5 6 7 8" nice python3 a30_e1.py run > e1_run.log 2>&1
  A30_ROOT=$PWD/e1 A30_SEEDS="1 2 3 4 5 6 7 8" python3 a30_e1.py summarize
  python3 a30_analyze.py e1 e1/e1_rows.json > _RESULT_a30_e1.txt
fi
if [ "$WHAT" = e2 ] || [ "$WHAT" = all ]; then
  # 20 backgrounds (seeds 1..20), 4 at a time x 4 arms in parallel = 16 processes
  export A30_ROOT=$PWD/e2 A30_G=8000 A30_A=10000 A30_MUS="0 0.0003 0.001 0.003 0.01 0.03" A30_NARMS=3 A30_PAR=4
  seq 1 20 | xargs -P 4 -I{} sh -c 'nice python3 a30_e2.py run {} > e2/run_bg{}.log 2>&1'
  python3 a30_e2.py table $PWD/e2
  python3 a30_analyze.py e2 e2/arms_table.json > _RESULT_a30_e2.txt
fi
echo "해부 30 launch finished: $WHAT"
