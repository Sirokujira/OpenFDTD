#!/bin/sh
# pml_reflection_check.sh — PML (abc = 1) の吸収性能の検証 (CI 用)
#
# data/sample/pml_reflection.ofd (小さい領域 + PML) と、そこから作る
# 大きい参照領域 (Mur) を同じソルバーで走らせ、観測点の時系列を比べる。
# 参照の境界反射が戻る前 (150 ステップ) までの差が、PML 側の境界が作る
# 反射。各点で max|V_pml - V_ref| / max|V_ref| <= TOL を判定する。
# 方法・ゲートの根拠・実測値は pml_reflection.ofd のコメント参照。
#
# 使い方 : pml_reflection_check.sh <ofd 実行ファイル(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (tpa_slab_check.sh と同じ):
#   OFD_LAUNCHER … 実行ファイルの前に置くコマンド (例: "mpirun -n 2")
#   OFD_ARGS     … 実行ファイルに渡す追加オプション (例: "-cpu")
#   OFD_H5CHECK  … HDF5 の時系列を読む道具 (既定: ソルバーと同じ場所の
#                  ofd_h5check。h5dump は Windows / macOS に無いので使わない)

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
SRC="$(cd "$(dirname "$0")" && pwd)/pml_reflection.ofd"
TOL=0.025
NSTEP=150

if [ -z "$OFD" ]; then
	echo "Usage: pml_reflection_check.sh <ofd> [workdir]" >&2
	exit 2
fi

H5CHECK="${OFD_H5CHECK:-$(dirname "$OFD")/ofd_h5check}"
if [ ! -x "$H5CHECK" ] && [ -x "$H5CHECK.exe" ]; then
	H5CHECK="$H5CHECK.exe"
fi
if [ ! -x "$H5CHECK" ]; then
	echo "*** ofd_h5check not found : $H5CHECK (cmake --build build --target ofd_h5check)" >&2
	exit 1
fi

mkdir -p "$WORK/pml" "$WORK/ref"

# 参照 : セル幅 (1 cm) と給電・観測点はそのまま、領域を ±0.52 m に広げ Mur にする
cp "$SRC" "$WORK/pml/pml.ofd"
sed -e 's/^xmesh = .*/xmesh = -0.52 104 0.52/' \
    -e 's/^ymesh = .*/ymesh = -0.52 104 0.52/' \
    -e 's/^zmesh = .*/zmesh = -0.52 104 0.52/' \
    -e 's/^abc = .*/abc = 0/' \
    "$SRC" > "$WORK/ref/ref.ofd"

for run in pml ref; do
	# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
	(cd "$WORK/$run" && $LAUNCHER "$OFD" $EXTRA -n 2 $run.ofd > /dev/null)
	grep -q "normal end" "$WORK/$run/ofd.log" || { echo "*** $run run failed" >&2; exit 1; }
	# データセットのパスは標準入力で渡す (Git Bash のパス変換対策, tools/h5check.c)
	echo /metadata/VPoint | "$H5CHECK" -d "$WORK/$run/time_series_data.h5" > "$WORK/$run/vpoint.txt"
done

grep -q "ABC = PML" "$WORK/pml/ofd.log" || { echo "*** PML is not active in pml run" >&2; exit 1; }

# VPoint は観測点ごとに (NSTEP + 1) 個ずつ並ぶ
paste "$WORK/ref/vpoint.txt" "$WORK/pml/vpoint.txt" | awk -v ns=$((NSTEP + 1)) -v tol=$TOL '
{
	p = int((NR - 1) / ns)
	if (p + 1 > np) np = p + 1
	d = $2 - $1; if (d < 0) d = -d
	a = $1;      if (a < 0) a = -a
	if (d > dm[p]) dm[p] = d
	if (a > am[p]) am[p] = a
	if (($2 ~ /[nN][aA][nN]|[iI][nN][fF]/) || ($2 + 0 > 1e30) || ($2 + 0 < -1e30)) bad = 1
}
END {
	status = 0
	if (bad) { print "*** non-finite value in PML run"; status = 1 }
	if (np != 4) { printf "*** expected 4 points, got %d\n", np; status = 1 }
	for (p = 0; p < np; p++) {
		r = (am[p] > 0) ? dm[p] / am[p] : 1
		ok = (r <= tol)
		printf "point %d : max|ref| = %.3e  max|pml-ref| = %.3e  ratio = %.4f -> %s\n", p + 1, am[p], dm[p], r, ok ? "OK" : "NG"
		if (!ok) status = 1
	}
	exit status
}' && status=0 || status=1

if [ "$status" -ne 0 ]; then
	echo "*** PML reflection check FAILED (tolerance ${TOL})" >&2
else
	echo "PML reflection check passed (tolerance ${TOL})"
fi
exit $status
