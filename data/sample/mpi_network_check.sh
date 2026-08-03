#!/bin/sh
# mpi_network_check.sh — ノード間 (PC 跨ぎ) MPI 相当の経路で結果が変わらないことの検証 (CI 用)
#
# MPI は同一ノード内では共有メモリ (vader/sm) でメッセージを渡すため、
# ノードを跨いだときだけ通る TCP の経路が単一マシンの試験では踏まれない。
# --mca btl tcp,self で共有メモリ転送を無効化すると、同じマシンのまま
# **物理的に別ノードへ送るときと同じ転送経路**を通せる。
#
# ここで捕まえたいのは、実際に PC を跨いだときだけ出る類のバグ:
#   - 送受信のバイト数・型の食い違い (共有メモリだと偶然通ることがある)
#   - 送信バッファの詰め方と受信側の解釈のずれ
#   - ランク条件の中に集団操作が入っていることによるデッドロック
# スナップショットの集約 (mpi/comm.c の comm_snapshot) は 1 ステップごとに
# 全ランクから rank 0 へ実データを送るので、ここが主な検査対象になる。
#
# 判定 : 1 プロセス実行を基準に、TCP 経路の複数プロセス実行で
#        ofd.log の収束履歴と入力インピーダンス表が完全一致すること。
#        HDF5 が読める環境なら /timeseries と /freqdomain も h5diff で比較する。
#
# 使い方 : mpi_network_check.sh <ソルバー(絶対パス)> <作業ディレクトリ> [追加オプション...]
#   例 : mpi_network_check.sh /path/to/bin/ofd_mpi      /tmp/w
#        mpi_network_check.sh /path/to/bin/ofd_cuda_mpi /tmp/w -cpu
#
# 環境変数:
#   MPIRUN      … mpirun コマンド (既定 mpirun)
#   MPIRUN_OPTS … mpirun のオプション (既定 --oversubscribe)
#   BTL_OPTS    … 転送経路の指定 (既定 --mca btl tcp,self)
#   NPROC_LIST  … 試すプロセス数 (既定 "2 4")

set -e

if [ $# -lt 2 ]; then
	echo "Usage: mpi_network_check.sh <solver> <workdir> [extra solver options...]" >&2
	exit 2
fi

SOLVER="$1"
WORK="$2"
shift 2
EXTRA="$*"

MPIRUN=${MPIRUN:-mpirun}
MPIRUN_OPTS=${MPIRUN_OPTS:---oversubscribe}
BTL_OPTS=${BTL_OPTS:---mca btl tcp,self}
NPROC_LIST=${NPROC_LIST:-2 4}

SRC=$(dirname "$0")/dipole.ofd
if [ ! -f "$SRC" ]; then
	echo "*** not found : $SRC" >&2
	exit 2
fi

mkdir -p "$WORK"
cp "$SRC" "$WORK"/
cd "$WORK"

# ログから比較対象を取り出す
#   収束履歴 : "<itime> <E> <H>" の 3 数値行
#   Zin 表   : 数値だけの行 (周波数・R・X ...)
extract() {
	grep -E '^[[:space:]]*[0-9]' "$1" | sed 's/[[:space:]]\+/ /g;s/^ //;s/ $//'
}

# 基準 : 1 プロセス (共有メモリも TCP も使わない)
$MPIRUN $MPIRUN_OPTS -n 1 "$SOLVER" $EXTRA -n 1 dipole.ofd > /dev/null 2>&1
extract ofd.log > ref.txt
if [ -f time_series_data.h5 ]; then
	cp time_series_data.h5 ref.h5
fi
echo "reference : 1 process ($(wc -l < ref.txt) lines)"

status=0
for np in $NPROC_LIST; do
	$MPIRUN $MPIRUN_OPTS $BTL_OPTS -n "$np" "$SOLVER" $EXTRA -n 1 dipole.ofd > /dev/null 2>&1
	extract ofd.log > "out_$np.txt"

	if cmp -s ref.txt "out_$np.txt"; then
		echo "tcp -n $np : ofd.log 一致 -> OK"
	else
		echo "tcp -n $np : ofd.log が基準と異なる -> NG"
		diff ref.txt "out_$np.txt" | head -10 || true
		status=1
	fi

	# HDF5 も比較する (h5diff があるときだけ)
	if [ -f ref.h5 ] && [ -f time_series_data.h5 ] && command -v h5diff > /dev/null 2>&1; then
		for ds in /timeseries/E /timeseries/H /freqdomain/E /freqdomain/H; do
			if h5diff ref.h5 time_series_data.h5 "$ds" "$ds" > /dev/null 2>&1; then
				echo "           $ds 一致 -> OK"
			else
				echo "           $ds が基準と異なる -> NG"
				status=1
			fi
		done
	fi
done

if [ $status -ne 0 ]; then
	echo "*** MPI network-path check failed" >&2
	exit 1
fi

echo "MPI network-path check passed"
