#!/bin/sh
# tpa_slab_check.sh — TPA スラブ検証 (CI 用)
#
# data/sample/tpa_slab.ofd の waveamp を 3 点掃引して実行し、
# ofd.log の "TPA: transmission = ..." を解析解 T = 1/(1 + β I0 L) と
# 比較する (許容誤差 ±7%、期待値の導出は tpa_slab.ofd のコメント参照)。
#
# 使い方 : tpa_slab_check.sh <ofd 実行ファイル(絶対パス)> [作業ディレクトリ]

set -e

OFD="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)/tpa_slab.ofd"
TOL=0.07

if [ -z "$OFD" ]; then
	echo "Usage: tpa_slab_check.sh <ofd> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"

status=0
# "waveamp 解析解T" の組 (tpa_slab.ofd : β=2000cm/GW, L=2um, I0=(1/2)ε0c E0^2)
for pair in "5.0e7 0.882830" "1.0e8 0.653217" "1.5e8 0.455687"; do
	amp=${pair% *}
	texp=${pair#* }
	sed "s/^waveamp = .*/waveamp = $amp/" "$SRC" > "$WORK/tpa_run.ofd"
	(cd "$WORK" && "$OFD" -n 2 tpa_run.ofd > /dev/null)
	t=$(grep "TPA: transmission" "$WORK/ofd.log" | tail -1 | awk '{print $4}')
	if [ -z "$t" ]; then
		echo "*** no 'TPA: transmission' line in $WORK/ofd.log" >&2
		exit 1
	fi
	res=$(awk -v t="$t" -v e="$texp" -v tol="$TOL" \
		'BEGIN{d=(t-e)/e; a=(d<0)?-d:d; printf "%s %+.2f%%", (a<=tol)?"OK":"NG", d*100}')
	echo "waveamp=$amp T_fdtd=$t T_analytic=$texp -> $res"
	case "$res" in NG*) status=1 ;; esac
done

if [ "$status" -ne 0 ]; then
	echo "*** TPA validation FAILED (tolerance ${TOL})" >&2
else
	echo "TPA validation passed (tolerance ${TOL})"
fi
exit $status
