#!/bin/sh
# hdf5_layout_check.sh — time_series_data.h5 の構成が契約どおりかの検証 (CI 用)
#
# このファイルには読み手が 2 つある:
#   1. 表示側 (OpenFDTD-X) — include/ofd_hdf5.h に書いた構成が契約
#   2. post/readhdf5.c     — HDF5 移行用の読み込み (既定では呼ばれない)
#
# どちらも「どのグループのどの名前にあるか」を決め打ちで開くので、
# 書き手 (sol/outputHdf5.c) 側でグループを移したり名前を変えたりすると
# 静かに壊れる。post/readhdf5.c は既定で呼ばれないため、実際に
# /geometry や /convergence へ移したときも長い間気付けなかった。
#
# ここではソルバーを 1 回走らせ、**必要なデータセットが期待するパスに
# 存在するか**を h5ls で確認する。値の正しさは別の検証 (mpi_network_check.sh
# など) が見るので、こちらは構成だけを見る。
#
# 使い方 : hdf5_layout_check.sh <ソルバー(絶対パス)> <作業ディレクトリ> [追加オプション...]
#   例 : hdf5_layout_check.sh /path/to/bin/ofd /tmp/w
#
# 環境変数:
#   OFD_LAUNCHER … 実行ファイルの前に置くコマンド (MPI 実行など)
#   OFD_ARGS     … 追加オプション

set -e

if [ $# -lt 2 ]; then
	echo "Usage: hdf5_layout_check.sh <solver> <workdir> [extra solver options...]" >&2
	exit 2
fi

SOLVER="$1"
WORK="$2"
shift 2
EXTRA="$*"

if ! command -v h5ls > /dev/null 2>&1; then
	echo "*** h5ls not found (hdf5-tools). skipped." >&2
	exit 0
fi

SRC=$(dirname "$0")/dipole.ofd
if [ ! -f "$SRC" ]; then
	echo "*** not found : $SRC" >&2
	exit 2
fi

mkdir -p "$WORK"
cp "$SRC" "$WORK"/
cd "$WORK"

# shellcheck disable=SC2086
$OFD_LAUNCHER "$SOLVER" $OFD_ARGS $EXTRA -n 1 dipole.ofd > /dev/null 2>&1

H5=time_series_data.h5
if [ ! -f "$H5" ]; then
	echo "*** $H5 was not created" >&2
	exit 1
fi

# 存在するデータセットの一覧 (フルパス)
h5ls -r "$H5" 2>/dev/null | awk '$2 == "Dataset" { print $1 }' | sort > have.txt

status=0
check() {
	if grep -qx "$1" have.txt; then
		echo "  OK   $1"
	else
		echo "  MISS $1"
		status=1
	fi
}

echo "--- 表示側の契約 (include/ofd_hdf5.h) ---"
for ds in \
	/geometry/Xn /geometry/Yn /geometry/Zn \
	/geometry/Xc /geometry/Yc /geometry/Zc \
	/timeseries/E /timeseries/H \
	/timeseries/itime /timeseries/time /timeseries/time_H \
	/freqdomain/E /freqdomain/H /freqdomain/freq \
	/convergence/iter /convergence/E /convergence/H
do
	check "$ds"
done

echo "--- post/readhdf5.c が開くもの ---"
for ds in \
	/metadata/Title /metadata/Dt \
	/metadata/Nx /metadata/Ny /metadata/Nz \
	/metadata/Ni /metadata/Nj /metadata/Nk /metadata/N0 /metadata/NN \
	/metadata/NFreq1 /metadata/NFreq2 /metadata/NFeed /metadata/NPoint \
	/metadata/Niter /metadata/Ntime /metadata/NGline \
	/metadata/Solver_maxiter /metadata/Solver_nout \
	/metadata/IPlanewave /metadata/Planewave /metadata/Planewave_pol \
	/metadata/Freq1 /metadata/Freq2 \
	/metadata/VFeed /metadata/IFeed \
	/metadata/NSurface
do
	check "$ds"
done

# 給電がある構成では Zin の表も書かれる (sol/outputZin.c が追記する)
check /metadata/input_impedance

if [ $status -ne 0 ]; then
	echo "*** HDF5 layout check failed" >&2
	echo "    sol/outputHdf5.c (書き手) と include/ofd_hdf5.h / post/readhdf5.c" >&2
	echo "    (読み手) のどちらかだけを変えていないか確認すること。" >&2
	exit 1
fi

echo "HDF5 layout check passed"
