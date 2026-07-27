#!/bin/sh
# thermal_material_check.sh — 熱解析レイヤのセル毎材料参照の検証 (CI 用)
#
# data/sample/thermal_slab.ofd を材料定数だけ変えて 4 通り実行し、
# ofd.log の "Thermal: dissipated[0] = ..." が次を満たすか判定する。
#
#   (1) σ = 0.05 S/m の導電スラブがある     -> 正の値
#   (2) 同じ材料を別の material id に付け替え -> (1) と完全一致
#   (3) スラブの σ を 0 にする               -> 0
#   (4) スラブの σ を 2 倍にする             -> (1) の 2 倍 (±1%)
#
# (1)(3) は「セルが真空 (材料 0) 固定ではない」こと、(2) は「id そのもの
# ではなくセルが指す材料が使われている」こと、(4) は損失が弱く場がほぼ
# 変わらない領域で p = (1/2)σ|E|^2 の σ 依存性が出ていることを見る。
# いずれも DFT の正規化定数に依らない性質なので、絶対値の較正なしに
# 判定できる。
#
# 使い方 : thermal_material_check.sh <ofd 実行ファイル(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (tpa_slab_check.sh と同じ):
#   OFD_LAUNCHER … 実行ファイルの前に置くコマンド (MPI 版の検証用)
#   OFD_ARGS     … 実行ファイルに渡す追加オプション

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRC="$(cd "$(dirname "$0")" && pwd)/thermal_slab.ofd"
TOL=0.01

if [ -z "$OFD" ]; then
	echo "Usage: thermal_material_check.sh <ofd> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"

# $1 = 入力ファイル名 -> "Thermal: dissipated[0]" の値を返す
run_case() {
	# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
	(cd "$WORK" && $LAUNCHER "$OFD" $EXTRA -n 2 "$1" > /dev/null)
	if ! grep -q "normal end" "$WORK/ofd.log"; then
		echo "*** $1 : no 'normal end' in ofd.log" >&2
		exit 1
	fi
	v=$(grep "Thermal: dissipated\[0\]" "$WORK/ofd.log" | tail -1 | awk '{print $4}')
	if [ -z "$v" ]; then
		echo "*** $1 : no 'Thermal: dissipated[0]' line in ofd.log" >&2
		exit 1
	fi
	echo "$v"
}

status=0

# (1) 基準 : σ = 0.05 S/m
cp "$SRC" "$WORK/th_base.ofd"
p_base=$(run_case th_base.ofd)
res=$(awk -v v="$p_base" 'BEGIN{print (v > 0) ? "OK" : "NG"}')
echo "(1) sigma=0.05        : P=$p_base -> $res (正であること)"
[ "$res" = OK ] || status=1

# (2) 同じ材料を別の id (3) に付け替える
sed 's/^geometry = 2 /geometry = 3 /' "$SRC" > "$WORK/th_id3.ofd"
p_id3=$(run_case th_id3.ofd)
res=$([ "$p_id3" = "$p_base" ] && echo OK || echo NG)
echo "(2) 同一材料を id=3 へ : P=$p_id3 -> $res ((1) と一致すること)"
[ "$res" = OK ] || status=1

# (3) σ = 0 -> 発熱なし
sed 's/^material = 1 4 0.05 1 0/material = 1 4 0 1 0/' "$SRC" > "$WORK/th_sig0.ofd"
p_sig0=$(run_case th_sig0.ofd)
res=$(awk -v v="$p_sig0" 'BEGIN{print (v == 0) ? "OK" : "NG"}')
echo "(3) sigma=0           : P=$p_sig0 -> $res (0 であること)"
[ "$res" = OK ] || status=1

# (4) σ を 2 倍 -> 発熱も 2 倍 (損失が弱いので場はほぼ変わらない)
sed 's/^material = 1 4 0.05 1 0/material = 1 4 0.10 1 0/' "$SRC" > "$WORK/th_sig2.ofd"
p_sig2=$(run_case th_sig2.ofd)
res=$(awk -v a="$p_sig2" -v b="$p_base" -v tol="$TOL" \
	'BEGIN{if (b <= 0) {print "NG"; exit} r=a/b; d=(r-2)/2; if (d<0) d=-d; printf "%s (P2/P1=%.4f)", (d<=tol)?"OK":"NG", r}')
echo "(4) sigma=0.10        : P=$p_sig2 -> $res ((1) の 2 倍であること)"
case "$res" in NG*) status=1 ;; esac

if [ "$status" -ne 0 ]; then
	echo "*** thermal per-cell material check FAILED" >&2
else
	echo "thermal per-cell material check passed"
fi
exit $status
