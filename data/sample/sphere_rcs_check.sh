#!/bin/sh
# sphere_rcs_check.sh — 平面波散乱の検証 (CI 用)
#
# data/sample/sphere_rcs.ofd を実行し、ofd.log の "=== cross section ===" の
# 後方 / 前方散乱断面積を完全導体球の Mie 厳密解と比較する。
#
#   ka = 3.0, a = 0.05 m
#   後方 σ = 4.0901e-03 m²   前方 σ = 8.4797e-02 m²
#   (導出は sphere_rcs.ofd の先頭コメント参照)
#
# 球を直交格子で階段近似しているため厳密一致はしない。ここで見たいのは
# 「平面波励振が生きていて桁と傾向が合っているか」なので、許容は 0.625〜1.6 倍と
# 広めに取る。0 / NaN / 桁違い (以前の不具合はすべて 0 だった) は確実に落ちる。
#
# 使い方 : sphere_rcs_check.sh <ofd 実行ファイル(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (tpa_slab_check.sh と同じ):
#   OFD_LAUNCHER … 実行ファイルの前に置くコマンド (MPI 版の検証用)
#   OFD_ARGS     … 実行ファイルに渡す追加オプション

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRC="$(cd "$(dirname "$0")" && pwd)/sphere_rcs.ofd"

# Mie 厳密解 [m^2] と許容比
MIE_BACK=4.0901e-03
MIE_FWD=8.4797e-02
LO=0.625
HI=1.6

if [ -z "$OFD" ]; then
	echo "Usage: sphere_rcs_check.sh <ofd> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
cp "$SRC" "$WORK/sphere_rcs.ofd"

# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
(cd "$WORK" && $LAUNCHER "$OFD" $EXTRA -n 2 sphere_rcs.ofd > /dev/null)

if ! grep -q "normal end" "$WORK/ofd.log"; then
	echo "*** no 'normal end' in $WORK/ofd.log" >&2
	exit 1
fi

line=$(grep -A2 "=== cross section ===" "$WORK/ofd.log" | tail -1)
back=$(echo "$line" | awk '{print $2}')
fwd=$(echo "$line" | awk '{print $3}')
if [ -z "$back" ] || [ -z "$fwd" ]; then
	echo "*** no cross section line in $WORK/ofd.log" >&2
	exit 1
fi

status=0
for pair in "backward $back $MIE_BACK" "forward $fwd $MIE_FWD"; do
	name=$(echo "$pair" | awk '{print $1}')
	val=$(echo "$pair"  | awk '{print $2}')
	ref=$(echo "$pair"  | awk '{print $3}')
	res=$(awk -v v="$val" -v r="$ref" -v lo="$LO" -v hi="$HI" 'BEGIN{
		# NaN は自分自身と等しくないことで検出する (0 除算の症状)
		if (v != v)          { printf "NG (NaN)"; exit }
		if (v <= 0)          { printf "NG (<= 0)"; exit }
		q = v / r
		printf "%s (Mie 比 %.3f)", ((q >= lo) && (q <= hi)) ? "OK" : "NG", q
	}')
	echo "$name : sigma=$val  Mie=$ref -> $res"
	case "$res" in NG*) status=1 ;; esac
done

if [ "$status" -ne 0 ]; then
	echo "*** plane wave RCS check FAILED (許容 ${LO}〜${HI} 倍)" >&2
else
	echo "plane wave RCS check passed (許容 ${LO}〜${HI} 倍)"
fi
exit $status
