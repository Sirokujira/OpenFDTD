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

# 検査の実体は 2 通り。どちらも同じパスの一覧を見る。
#   h5ls        … hdf5-tools が入っていれば使う
#   ofd_h5check … 入っていないときの代替 (ソルバーと同じ HDF5 を
#                 リンクした小さな道具。tools/h5check.c)
# Windows の CI は vcpkg の静的ビルドでツール類が PATH に来ないので、
# 後者が無いと**中身が何も検証されない**まま通ってしまう。
# 既定はソルバーと同じディレクトリを見る (bin/ に並んで置かれる)。
H5CHECK=${OFD_H5CHECK:-$(dirname "$SOLVER")/ofd_h5check}
case "$SOLVER" in *.exe) H5CHECK="$H5CHECK.exe" ;; esac

if command -v h5ls > /dev/null 2>&1; then
	MODE=h5ls
elif [ -x "$H5CHECK" ]; then
	MODE=h5check
else
	echo "*** h5ls も $H5CHECK も見つからない。検査できない。" >&2
	echo "    hdf5-tools を入れるか、ofd_h5check をビルドすること" >&2
	echo "    (cmake --build <build> --target ofd_h5check)。" >&2
	exit 1
fi
echo "検査方法 : $MODE"

# dipole.ofd を使う。形状 (geometry) を 1 つ持つので
# /geometry/Gline と /geometry/MGline も作られる
# (この 2 つは NGline == 0 のとき作られない。3D シーンで構造を
#  重ねて描くのに使うデータなので、欠けても実行は正常終了してしまい
#  「表示だけ出ない」形で静かに壊れる)。
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

# 検査するパスの一覧。両方式で同じものを使う。
#
# 表示側の契約 (include/ofd_hdf5.h)。ここが正で、GUI (OpenFDTD-X) は
# このパスを決め打ちで開く。/geometry/Gline と MGline は 3D シーンに
# 構造を重ねて描くのに使う (欠けても実行は normal end で終わるので、
# 「表示だけ出ない」形で静かに壊れる)。
CONTRACT="
/geometry/Xn /geometry/Yn /geometry/Zn
/geometry/Xc /geometry/Yc /geometry/Zc
/geometry/Gline /geometry/MGline
/timeseries/E /timeseries/H
/timeseries/itime /timeseries/time /timeseries/time_H
/freqdomain/E /freqdomain/H /freqdomain/freq
/convergence/iter /convergence/E /convergence/H
"

# post/readhdf5.c が開くもの。
#
# 注意 : readhdf5() は /metadata/Surface も開こうとするが、書き手
# (sol/outputHdf5.c) は Surface を書いていない (表示に不要なため PR #18 で
# 意図的に見送った)。したがって readhdf5() は現状そのままでは完結しない。
# ここに Surface を足すと「書いていないものを要求する」検査になるので入れない。
# HDF5 移行を完了させるときは、Surface を書く側に足してから追加すること。
#
# 末尾の input_impedance は給電がある構成でのみ書かれる
# (sol/outputZin.c が追記する)。dipole.ofd は給電があるので必ず出る。
READER="
/metadata/Title /metadata/Dt
/metadata/Nx /metadata/Ny /metadata/Nz
/metadata/Ni /metadata/Nj /metadata/Nk /metadata/N0 /metadata/NN
/metadata/NFreq1 /metadata/NFreq2 /metadata/NFeed /metadata/NPoint
/metadata/Niter /metadata/Ntime /metadata/NGline
/metadata/Solver_maxiter /metadata/Solver_nout
/metadata/IPlanewave /metadata/Planewave /metadata/Planewave_pol
/metadata/Freq1 /metadata/Freq2
/metadata/VFeed /metadata/IFeed
/metadata/NSurface
/metadata/input_impedance
"

status=0

if [ "$MODE" = "h5ls" ]; then
	# 存在するデータセットの一覧 (フルパス) を作って突き合わせる
	h5ls -r "$H5" 2>/dev/null | awk '$2 == "Dataset" { print $1 }' | sort > have.txt
	check_all() {
		for ds in $2; do
			if grep -qx "$ds" have.txt; then
				echo "  OK   $ds"
			else
				echo "  MISS $ds"
				status=1
			fi
		done
	}
else
	# ofd_h5check にまとめて渡す (OK/MISS を同じ書式で出す)。
	#
	# **パスは必ず標準入力で渡すこと。** Git Bash (MSYS) は先頭が '/' の
	# 引数をネイティブ Windows 実行ファイルへ渡すとき勝手に Windows パスへ
	# 書き換えるので、コマンドライン引数だと "/metadata/Niter" が
	# "C:/Program Files/Git/metadata/Niter" になって全部 MISS になる
	# (CI の Windows ジョブで実際に踏んだ)。標準入力は変換されない。
	check_all() {
		printf '%s\n' $2 | "$H5CHECK" "$H5" || status=1
	}
fi

echo "--- 表示側の契約 (include/ofd_hdf5.h) ---"
check_all "contract" "$CONTRACT"

echo "--- post/readhdf5.c が開くもの ---"
check_all "reader" "$READER"

if [ $status -ne 0 ]; then
	echo "*** HDF5 layout check failed" >&2
	echo "    sol/outputHdf5.c (書き手) と include/ofd_hdf5.h / post/readhdf5.c" >&2
	echo "    (読み手) のどちらかだけを変えていないか確認すること。" >&2
	exit 1
fi

echo "HDF5 layout check passed"
