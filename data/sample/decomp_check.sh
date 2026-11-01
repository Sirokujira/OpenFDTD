#!/bin/sh
# decomp_check.sh — MPI の領域分割によらず結果が一致することの検証 (3 次元の問題、CI 用)
#
# tpa_decomp_check.sh は TPA の平板 (z 方向に進む平面波) だけを見るので、
# x / y 方向に一様な場しか通らない。x / y 分割のハロー交換を壊しても
# (Hz を受け取らない、送る面を 1 つずらす) 透過率は変わらず、検出できなかった
# (2026-11 の変異試験)。ここでは場が 3 次元に広がるサンプルを分割を変えて走らせ、
# 1 プロセスの結果と次が一致することを確かめる:
#   - ofd.log (プロセス数・スレッド数・メモリ量・日時・所要時間の行を除く。収束履歴・
#     インピーダンス表などが入る)
#   - ofd_post の *.log (遠方界・近傍界・給電・観測点の数値)
#
# 収束履歴の平均値 (<E> <H>) は各ランクの部分和を Allreduce で足すので、加算順で
# 揺らぐ。CPU / MPI は倍精度で 15 桁目 (tpa_half の -p 1 1 2)、CUDA 版は単精度で
# 足すので 7 桁目 (相対 1e-6) まで動く。収束履歴だけは相対 1e-5 で比べ、他の行は
# 完全一致を求める (平均値の Allreduce を壊すと桁違いにずれるので検出できる)。
#
# サンプル : dipole (給電・Mur) / sphere (平面波) / dipole_couple (PML・観測点) /
#   loop_load (負荷) / pbc (周期境界 x, y) / pml_reflection (PML・観測点 4 つ) /
#   tpa_half (tpa_slab の TPA 層を x / y とも半分 (1/4 の柱) にしたもの。場が x / y に
#   一様でなくなり Ey / Ez も立つので、TPA だけが使う E のハロー交換 mpi/comm_E.c と
#   成分の間の交換の順序 (mpi/solve.c) を通る。x だけ半分では Ey が対称性で 0 になり、
#   Ey の交換を省いても見えなかった)
#
# 使い方 : decomp_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数:
#   MPIRUN      … mpirun コマンド (既定 mpirun)
#   MPIRUN_OPTS … mpirun のオプション (既定 --oversubscribe)
#   OFD_ARGS    … ソルバーに渡す追加オプション (例: "-cpu")
#   OFD_POST    … 後処理 (既定: ソルバーと同じ場所の ofd_post)
#   DECOMP_LIST … 試す分割の一覧 (既定は下記。先頭が基準)

set -e

OFD="$1"
WORK="${2:-.}"
EXTRA="${OFD_ARGS:-}"
MPIRUN="${MPIRUN:-mpirun}"
MPIRUN_OPTS="${MPIRUN_OPTS:---oversubscribe}"
DECOMP_LIST="${DECOMP_LIST:-1 1 1;2 1 1;1 2 1;1 1 2;2 2 2}"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OFD" ]; then
	echo "Usage: decomp_check.sh <ofd_mpi> [workdir]" >&2
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

# ofd.log からプロセス数・スレッド数・メモリ量 (rank 0 の分)・日時・所要時間の行を除く
filter_log() {
	grep -v -E -e 'process=' -e 'thread=' -e 'Memory size' -e '^  part-' -e '^  total ' \
	        -e '^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) ' "$1" |
	awk 'NF == 3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9.e+-]+$/ && $3 ~ /^[0-9.e+-]+$/ { next } { print }'
}

# 収束履歴 (step <E> <H>) の行だけ取り出す
conv_log() {
	awk 'NF == 3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9.e+-]+$/ && $3 ~ /^[0-9.e+-]+$/ { print }' "$1"
}

# 収束履歴を相対 1e-5 で比べる (行数・step が同じで、各値の相対差が 1e-5 以下)
conv_same() {
	awk 'NR == FNR { a[FNR] = $0; n = FNR; next }
	     { split(a[FNR], r, " ")
	       if ((FNR > n) || (r[1] != $1)) bad = 1
	       for (c = 2; c <= 3; c++) {
	         d = r[c] - $c; if (d < 0) d = -d
	         m = (r[c] < 0) ? -r[c] : r[c]
	         if (d > 1e-5 * m + 1e-30) bad = 1
	       }
	       m2 = FNR }
	     END { exit (bad || (m2 != n)) ? 1 : 0 }' "$1" "$2"
}

status=0
for b in dipole sphere dipole_couple loop_load pbc pml_reflection tpa_half; do
	ref=""
	OLDIFS=$IFS
	IFS=';'
	for p in $DECOMP_LIST; do
		IFS=$OLDIFS
		[ -n "$p" ] || continue
		n=$(echo "$p" | awk '{print $1 * $2 * $3}')
		tag=$(echo "$p" | tr -d ' ')
		d="$WORK/$b/p$tag"
		mkdir -p "$d"
		if [ "$b" = "tpa_half" ]; then
			sed 's/^geometry = 3 1 0 0\.4e-6 0 0\.4e-6 /geometry = 3 1 0 0.2e-6 0 0.2e-6 /' "$SRCDIR/tpa_slab.ofd" > "$d/run.ofd"
			grep -q '^geometry = 3 1 0 0\.2e-6 0 0\.2e-6 ' "$d/run.ofd" || { echo "*** tpa_half : tpa_slab.ofd の TPA 層の行が見つからない" >&2; exit 1; }
		else
			cp "$SRCDIR/$b.ofd" "$d/run.ofd"
		fi
		# $MPIRUN_OPTS / $EXTRA / $p は複数語に分解させたいので意図的にクォートしない
		(cd "$d" && $MPIRUN $MPIRUN_OPTS -n "$n" "$OFD" $EXTRA -p $p -n 1 run.ofd > /dev/null 2>&1) || true
		if ! grep -q "normal end" "$d/ofd.log" 2>/dev/null; then
			echo "*** $b -p $p : solver failed" >&2
			status=1
			IFS=';'
			continue
		fi
		(cd "$d" && "$POST" -n 1 run.ofd > post.txt 2>&1) || true
		filter_log "$d/ofd.log" > "$d/ofd.filtered"
		conv_log "$d/ofd.log" > "$d/ofd.conv"
		if [ -z "$ref" ]; then
			ref="$d"
			IFS=';'
			continue
		fi
		bad=""
		cmp -s "$ref/ofd.filtered" "$d/ofd.filtered" || bad="$bad ofd.log"
		conv_same "$ref/ofd.conv" "$d/ofd.conv" || bad="$bad convergence"
		nlog=1
		for f in "$ref"/*.log; do
			g=$(basename "$f")
			[ "$g" = "ofd.log" ] && continue
			nlog=$((nlog + 1))
			cmp -s "$f" "$d/$g" || bad="$bad $g"
		done
		if [ -n "$bad" ]; then
			echo "$b -p $p : differ from -p $(basename "$ref" | sed 's/^p//') :$bad -> NG"
			{ diff "$ref/ofd.filtered" "$d/ofd.filtered"; diff "$ref/ofd.conv" "$d/ofd.conv"; } | head -8 || true
			status=1
		else
			echo "$b -p $p : $nlog logs identical -> OK"
		fi
		IFS=';'
	done
	IFS=$OLDIFS
done

if [ "$status" -ne 0 ]; then
	echo "*** decomposition check FAILED"
else
	echo "decomposition check passed"
fi
exit $status
