#!/bin/sh
# tpa_decomp_check.sh — MPI 領域分割によらず結果が一致することの検証 (CI 用)
#
# data/sample/tpa_slab.ofd を複数の領域分割で実行し、ofd.log の
# "TPA: transmission = ..." が全ケースで完全一致することを判定する。
#
# 解析解との比較は tpa_slab_check.sh が行う。こちらが見るのは
# 「プロセス数・分割方向を変えても同じ答えになるか」だけで、
# ハロー交換の漏れ・入力の配り忘れ (comm_broadcast) といった
# 並列化固有のバグを捕まえるのが目的。基準値は 1 プロセス実行の結果。
#
# 使い方 : tpa_decomp_check.sh <ソルバー(絶対パス)> <作業ディレクトリ> [追加オプション...]
#   例 : tpa_decomp_check.sh /path/to/bin/ofd_mpi      /tmp/w
#        tpa_decomp_check.sh /path/to/bin/ofd_cuda_mpi /tmp/w -cpu
#
# 環境変数:
#   MPIRUN      … mpirun コマンド (既定 mpirun)
#   MPIRUN_OPTS … mpirun のオプション (既定 --oversubscribe)
#   DECOMP_LIST … 試す分割の一覧 (既定は下記 5 通り)

set -e

if [ $# -lt 2 ]; then
	echo "Usage: tpa_decomp_check.sh <solver> <workdir> [extra solver options...]" >&2
	exit 2
fi

SOLVER="$1"
WORK="$2"
shift 2
EXTRA="$*"

SRC="$(cd "$(dirname "$0")" && pwd)/tpa_slab.ofd"
MPIRUN="${MPIRUN:-mpirun}"
MPIRUN_OPTS="${MPIRUN_OPTS:---oversubscribe}"
# 1 プロセス / 各軸 1 方向のみ / 3 軸すべて の代表 5 通り
DECOMP_LIST="${DECOMP_LIST:-1 1 1;1 1 2;2 1 1;1 2 1;2 2 2}"

mkdir -p "$WORK"
cp "$SRC" "$WORK/tpa_run.ofd"

ref=""
status=0

# ";" 区切りの分割リストを 1 件ずつ処理する
# (パイプ + while にするとサブシェルになり ref/status を持ち出せないので
#  IFS を切り替えた for ループにする)
OLDIFS=$IFS
IFS=';'
for p in $DECOMP_LIST; do
	IFS=$OLDIFS
	[ -n "$p" ] || continue
	n=$(echo "$p" | awk '{print $1 * $2 * $3}')

	# $MPIRUN_OPTS / $EXTRA は複数語に分解させたいので意図的にクォートしない
	(cd "$WORK" && $MPIRUN $MPIRUN_OPTS -n "$n" "$SOLVER" $EXTRA -p $p -n 1 tpa_run.ofd > /dev/null)

	if ! grep -q "normal end" "$WORK/ofd.log"; then
		echo "*** -p $p ($n procs) : no 'normal end' in ofd.log" >&2
		exit 1
	fi
	t=$(grep "TPA: transmission" "$WORK/ofd.log" | tail -1 | awk '{print $4}')
	if [ -z "$t" ]; then
		echo "*** -p $p ($n procs) : no 'TPA: transmission' line in ofd.log" >&2
		exit 1
	fi

	if [ -z "$ref" ]; then
		ref="$t"
		echo "-p $p ($n procs) : T=$t (基準)"
	elif [ "$t" = "$ref" ]; then
		echo "-p $p ($n procs) : T=$t -> OK"
	else
		echo "-p $p ($n procs) : T=$t -> NG (基準 $ref と不一致)" >&2
		status=1
	fi
	IFS=';'
done
IFS=$OLDIFS

if [ "$status" -ne 0 ]; then
	echo "*** TPA decomposition check FAILED" >&2
else
	echo "TPA decomposition check passed"
fi
exit "$status"
