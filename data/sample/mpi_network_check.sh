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
#        ofd.log の収束履歴と入力インピーダンス表、および (h5diff があれば)
#        /timeseries と /freqdomain が一致すること。
#
# 完全一致ではなく相対誤差で見る。収束履歴は comm_average() の Allreduce に
# よる総和なので、加算順序がプロセス数で変わり最終桁が動く
# (ofd_cuda_mpi では実際に 3e-6 程度ずれる。ofd_mpi では偶然一致していた)。
# 一方、転送が壊れた場合の差は桁違いに大きいので、この許容幅でも検出できる。
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
#   TOL         … 相対許容誤差 (既定 1e-4)

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
TOL=${TOL:-1e-4}

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

# 数値を相対誤差で比較する (行数・列数が違えば即 NG)
compare_num() {
	awk -v tol="$3" '
	NR == FNR { ref[FNR] = $0; nref = FNR; next }
	{
		if (FNR > nref) { extra = 1; next }
		n = split(ref[FNR], a); m = split($0, b);
		if (n != m) { printf "  列数が違う (行 %d): %d vs %d\n", FNR, n, m; bad = 1; next }
		for (i = 1; i <= n; i++) {
			r = a[i] + 0; c = b[i] + 0;
			d = (r > c) ? r - c : c - r;
			s = (r < 0) ? -r : r;
			if (s < 1) s = 1;
			if (d / s > tol) {
				printf "  行 %d 列 %d : %s vs %s (相対 %.3g)\n", FNR, i, a[i], b[i], d / s;
				bad = 1;
			}
		}
		ncur = FNR;
	}
	END {
		if (extra || (ncur != nref)) { printf "  行数が違う : %d vs %d\n", nref, FNR; bad = 1 }
		exit bad ? 1 : 0
	}' "$1" "$2"
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

	if compare_num ref.txt "out_$np.txt" "$TOL" > cmp_$np.txt 2>&1; then
		echo "tcp -n $np : ofd.log 一致 (相対誤差 <= $TOL) -> OK"
	else
		echo "tcp -n $np : ofd.log が基準と異なる -> NG"
		head -10 cmp_$np.txt || true
		status=1
	fi

	# HDF5 も比較する (h5diff があるときだけ)
	if [ -f ref.h5 ] && [ -f time_series_data.h5 ] && command -v h5diff > /dev/null 2>&1; then
		for ds in /timeseries/E /timeseries/H /freqdomain/E /freqdomain/H; do
			if h5diff -p "$TOL" ref.h5 time_series_data.h5 "$ds" "$ds" > /dev/null 2>&1; then
				echo "           $ds 一致 -> OK"
			else
				echo "           $ds が基準と異なる -> NG"
				status=1
			fi
		done
	fi
done

# HDF5 を作れない状況でもデッドロックしないこと
#
# 集約 (comm_snapshot) は全ランクが参加する送受信なので、rank 0 が
# 「書けないから」と途中で抜けると他ランクが MPI_Send で待ち続けて
# ハングする。実際にこの経路でデッドロックしたことがある
# (症状は「1 プロセスなら正常、2 プロセス以上でハング」)。
#
# time_series_data.h5 と同名のディレクトリを置くと H5Fcreate だけが失敗し、
# ログなど他の出力はそのまま書けるので、この経路だけを再現できる。
echo "--- HDF5 を作れない場合 (デッドロックしないこと) ---"
rm -rf blocked && mkdir -p blocked
cp dipole.ofd blocked/
(
	cd blocked
	mkdir -p time_series_data.h5
	# タイムアウトが使えない環境では検査を飛ばす
	if command -v timeout > /dev/null 2>&1; then
		# shellcheck disable=SC2086
		if timeout 120 $MPIRUN $MPIRUN_OPTS $BTL_OPTS -n 2 "$SOLVER" $EXTRA -n 1 dipole.ofd > /dev/null 2>&1; then
			:
		elif [ $? -eq 124 ]; then
			echo "  ハングした (デッドロック) -> NG"
			exit 1
		fi
		if grep -q "normal end" ofd.log 2>/dev/null; then
			echo "  最後まで完走した -> OK"
		else
			echo "  normal end に達しなかった -> NG"
			exit 1
		fi
	else
		echo "  timeout が無いので省略"
	fi
) || status=1

if [ $status -ne 0 ]; then
	echo "*** MPI network-path check failed" >&2
	exit 1
fi

echo "MPI network-path check passed"
