#!/usr/bin/env bash
# 해부 6 본 발사 — C(침입) · B(예산) · A(계보) 차례로. 코어를 나눠 쓰지 않는다.
cd "$HOME/petri" || exit 1
./inv_launch.sh
./bud_launch.sh
./lin_launch.sh
echo "===== 해부6 전부 끝 $(date "+%F %T") =====" > run6.done
