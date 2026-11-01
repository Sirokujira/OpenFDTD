#!/bin/sh
# hdf5_post_check.sh — 後処理を HDF5 から読んでも ofd.out から読んだのと同じになるかの検証 (CI 用)
#
# ofd_post は既定で ofd.out (post/readout.c) を読み、-h5 を付けると
# time_series_data.h5 (post/readhdf5.c) を読む。HDF5 には /metadata・/geometry・
# /convergence に加えて、後処理用の /post (余白込みの近傍界 DFT・給電・入力電力。
# include/ofd_hdf5.h) をソルバーが書く。両方の経路で ofd_post の出力
# (*.log と作図の ev.ev2 / ev.ev3) がバイト単位で一致することを確かめる。
#
# サンプルは経路の違う 4 つ:
#   dipole        給電 (入力インピーダンス・遠方界の利得・近傍界の面)
#   dipole_couple 給電 + 観測点 (結合度)
#   loop_load     負荷付きの給電
#   sphere        平面波 (遠方界の散乱パターン・入射界を足した近傍界)
#
# 使い方 : hdf5_post_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (dipole_gain_check.sh と同じ):
#   OFD_LAUNCHER … ソルバーの前に置くコマンド (例: "mpirun -n 2")
#   OFD_ARGS     … ソルバーに渡す追加オプション (例: "-cpu")
#   OFD_POST     … 後処理 (既定: ソルバーと同じ場所の ofd_post)

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OFD" ]; then
	echo "Usage: hdf5_post_check.sh <ofd> [workdir]" >&2
	exit 2
fi

POST="${OFD_POST:-$(dirname "$OFD")/ofd_post}"
if [ ! -x "$POST" ] && [ -x "$POST.exe" ]; then
	POST="$POST.exe"
fi
if [ ! -x "$POST" ]; then
	echo "*** ofd_post not found : $POST" >&2
	exit 1
fi

status=0
for b in dipole dipole_couple loop_load sphere; do
	d="$WORK/$b"
	mkdir -p "$d/out" "$d/h5"
	cp "$SRCDIR/$b.ofd" "$d/out/run.ofd"
	# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
	(cd "$d/out" && $LAUNCHER "$OFD" $EXTRA -n 2 run.ofd > /dev/null)
	grep -q "normal end" "$d/out/ofd.log" || { echo "*** $b : solver failed" >&2; status=1; continue; }
	cp "$d/out/run.ofd" "$d/out/time_series_data.h5" "$d/h5/"
	(cd "$d/out" && "$POST" -n 2 run.ofd > post.txt 2>&1) || { echo "*** $b : ofd_post (ofd.out) failed" >&2; status=1; continue; }
	(cd "$d/h5"  && "$POST" -h5 -n 2 run.ofd > post.txt 2>&1) || { echo "*** $b : ofd_post -h5 failed" >&2; cat "$d/h5/post.txt" >&2; status=1; continue; }
	n=0
	bad=""
	for f in "$d/out"/*.log "$d/out"/ev.ev2 "$d/out"/ev.ev3; do
		[ -f "$f" ] || continue
		g=$(basename "$f")
		[ "$g" = "ofd.log" ] && continue     # ソルバーのログ (後処理の出力ではない)
		n=$((n + 1))
		cmp -s "$f" "$d/h5/$g" || bad="$bad $g"
	done
	if [ "$n" -eq 0 ]; then
		echo "*** $b : no post output to compare" >&2
		status=1
	elif [ -n "$bad" ]; then
		echo "$b : $n files, differ :$bad -> NG"
		status=1
	else
		echo "$b : $n files identical (ofd.out vs -h5) -> OK"
	fi
done

if [ "$status" -ne 0 ]; then
	echo "*** HDF5 post check FAILED"
else
	echo "HDF5 post check passed"
fi
exit $status
