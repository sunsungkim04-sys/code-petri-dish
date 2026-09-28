#!/usr/bin/env bash
# 해부 7 ② ③ 본 발사 — μ 0 침입(③b) · 재료 먼저 침입(②) · 계보 두 규칙(② ③a) 차례로
cd "$HOME/petri" || exit 1
./ib_launch.sh
./ff_inv_launch.sh
./ff_lin_launch.sh
echo "===== 해부7 끝 $(date "+%F %T") =====" > run7.done
