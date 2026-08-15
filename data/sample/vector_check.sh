#!/bin/sh
# vector_check.sh — ベクトル化更新経路 (-vector) と通常経路 (-no-vector) の一致検証 (CI 用)
#
# OpenFDTD の電磁界更新には **2 つの独立した実装**がある。
#
#   -no-vector (既定) … セル毎に材料 ID (iEx..iHz) を引き、材料テーブル
#                        C1/C2/D1/D2 の係数で更新する。PEC は if 分岐で
#                        特別扱いする (例 updateHx_p_no_vector: Hx = -fi)。
#   -vector           … setup_vector() が全セル分の係数を展開した配列
#                        K1Ex..K2Hz を先に作り、分岐なしで更新する。
#                        PEC は分岐ではなく K1 = K2 = 0 という係数で表現する
#                        (Hx = 0*Hx - 0*curl - 0*dfi - 1*fi = -fi と等価)。
#
# 対象は sol/update{Ex,Ey,Ez,Hx,Hy,Hz}.c の 12 関数 (_f_/_p_ × vector/no_vector)
# と sol/setup_vector.c、および CUDA 側の同名カーネル。**4 実装すべてが
# -vector を受け付ける** (sol / mpi / cuda / cuda_mpi の各 Main)。
#
# ここで捕まえたいのは、係数展開と分岐なし更新の等価性が崩れる類のバグ:
#   - K1/K2 (C1/C2, D1/D2) の取り違え。式の形が似ているので目視で気付けない
#     (Python 移植版で実際に K2Hx を書くべき所が K1Hx になっていた)
#   - setup_vector() の展開ループの範囲が更新ループの範囲と食い違う
#     (端の 1 面だけ係数が 0 のまま残る)
#   - PEC や分散性材料など、no_vector 側だけが持つ特別扱いの取りこぼし
# いずれもログは "normal end" で終わるのでスモークテストでは検出できない。
#
# 判定 : 同じ入力を -no-vector と -vector で実行し、
#        ofd.log の数値 (収束履歴・インピーダンス表・散乱断面積など) が
#        完全一致し、HDF5 の界データが相対誤差 TOL 以内で一致すること。
#
# ログは完全一致を要求する (実測でも全 27 サンプルで完全一致した)。
# HDF5 側だけ許容誤差を持たせているのは、2 つの経路で式の書き方が異なり
# コンパイラの FMA 縮約が変わるため、float32 で保存した界に 1 ULP
# (相対 6e-8 程度) の差が出るサンプルがあるから。既定の TOL = 1e-5 は
# その 100 倍以上の余裕がある一方、係数の取り違えのような実バグは
# 桁違いに大きな差を出すので取りこぼさない。
#
# 使い方 : vector_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数 (どれも省略時は既定動作):
#   OFD_LAUNCHER … 実行ファイルの前に置くコマンド
#                  例: OFD_LAUNCHER="mpirun --oversubscribe -n 2"
#   OFD_ARGS     … 実行ファイルに渡す追加オプション (例: "-cpu")
#   CASES        … 試すケース (既定は下を参照。SAMPLES でも可)
#   TOL          … HDF5 の相対許容誤差 (既定 1e-5)

set -e

OFD="$1"
WORK="${2:-.}"
LAUNCHER="${OFD_LAUNCHER:-}"
EXTRA="${OFD_ARGS:-}"
TOL=${TOL:-1e-5}

# 既定の実行ケース : 更新関数の分岐をひと通り踏む最小の組み合わせ。
# "<サンプル名>" か "<サンプル名>:pol<n>" (平面波の偏波を n に差し替えて実行)。
#
#   dipole                 … 給電 (NFeed) -> update*_f_* 側。PEC 細線
#   cube                   … 平面波 -> update*_p_* 側。PEC 直方体
#   debye                  … 平面波 + Debye 分散性材料 (setupDispersion)
#   square2d               … 平面波 + 周期境界 (pbc)
#   vector_lossy_feed      … 給電 + εr/μr が 1 でない材料
#   vector_lossy_pw:pol1/2 … 斜め入射の平面波 + 同じ材料 + PEC 板
#
# vector_lossy_* は必須。更新係数は C1 = εr/denom, C2 = 1/denom (H 側は
# D1 = μr/denom, D2 = 1/denom) なので、**真空 (εr = μr = 1, 無損失) では
# C1 = C2 = 1、PEC では C1 = C2 = 0** となり、C1 と C2 (D1 と D2) を
# 取り違えても結果が 1 ビットも変わらない。この 2 つを入れる前は、
# 既定ケースの dipole / cube / square2d がいずれも真空 + PEC のみ、
# debye も損失材料を geometry に配置していなかったため、係数の取り違えを
# 一切検出できなかった。
#
# 平面波を 2 偏波とも回すのも必須。sol/input2.c の設定では垂直偏波で
# 入射 H の z 成分が、水平偏波で入射 E の z 成分が、入射方向に依らず
# 厳密に 0 になるので、片方だけでは 6 成分の入射項を網羅できない。
# 詳細は各 .ofd の先頭コメント参照。
CASES=${CASES:-${SAMPLES:-"dipole cube debye square2d vector_lossy_feed vector_lossy_pw:pol1 vector_lossy_pw:pol2"}}

if [ -z "$OFD" ]; then
	echo "Usage: vector_check.sh <solver> [workdir]" >&2
	exit 2
fi

DATA=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)

# ログから経路によって当然変わる行を落とす
#   - 実行日時 (ctime の曜日始まり)
#   - "vector=on/off" を含む構成行
#   - メモリ使用量 (-vector は係数配列の分だけ増える)
#   - cpu time 以降 (所要時間)
filter_log() {
	awk '
	/^=== cpu time/       { exit }
	/^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) / { next }
	/vector=/             { next }
	/Memory size/         { next }
	                      { print }
	' "$1"
}

status=0
for c in $CASES; do
	# "<サンプル名>" または "<サンプル名>:pol<偏波>" (偏波は planewave の 3 番目)
	s=${c%%:*}
	pol=${c#*:}
	[ "$pol" = "$c" ] && pol=""

	src="$DATA/$s.ofd"
	if [ ! -f "$src" ]; then
		echo "*** not found : $src" >&2
		exit 2
	fi

	label="$s"
	tag="$s"
	if [ -n "$pol" ]; then
		label="$s (${pol})"
		tag="$s-$pol"
	fi

	for mode in no-vector vector; do
		d="$WORK/$tag-$mode"
		rm -rf "$d"
		mkdir -p "$d"
		if [ -n "$pol" ]; then
			# planewave = <θ> <φ> <偏波> の偏波だけ差し替える
			sed "s/^planewave = \\([^ ]*\\) \\([^ ]*\\) .*/planewave = \\1 \\2 ${pol#pol}/" "$src" > "$d/$s.ofd"
		else
			cp "$src" "$d/"
		fi
		# $LAUNCHER / $EXTRA は複数語に分解させたいので意図的にクォートしない
		(cd "$d" && $LAUNCHER "$OFD" $EXTRA "-$mode" -n 2 "$s.ofd" > /dev/null 2>&1)
		if ! grep -q "normal end" "$d/ofd.log" 2>/dev/null; then
			echo "$label ($mode) : normal end に達しなかった -> NG"
			status=1
			continue 2
		fi
	done

	a="$WORK/$tag-no-vector"
	b="$WORK/$tag-vector"
	s="$tag"

	filter_log "$a/ofd.log" > "$WORK/$s.a"
	filter_log "$b/ofd.log" > "$WORK/$s.b"
	if diff -q "$WORK/$s.a" "$WORK/$s.b" > /dev/null 2>&1; then
		echo "$label : ofd.log 完全一致 -> OK"
	else
		echo "$label : ofd.log が -vector と -no-vector で異なる -> NG"
		diff "$WORK/$s.a" "$WORK/$s.b" | head -10 || true
		status=1
	fi

	# HDF5 の界データも比較する (h5diff があるときだけ)
	if [ -f "$a/time_series_data.h5" ] && [ -f "$b/time_series_data.h5" ] \
	   && command -v h5diff > /dev/null 2>&1; then
		for ds in /timeseries/E /timeseries/H /freqdomain/E /freqdomain/H; do
			# 実装によっては書かれないデータセットがあるので、両方にある物だけ見る
			if ! h5ls "$a/time_series_data.h5$ds" > /dev/null 2>&1 \
			|| ! h5ls "$b/time_series_data.h5$ds" > /dev/null 2>&1; then
				continue
			fi
			if h5diff -p "$TOL" "$a/time_series_data.h5" "$b/time_series_data.h5" "$ds" "$ds" > /dev/null 2>&1; then
				echo "         $ds 一致 (相対 <= $TOL) -> OK"
			else
				echo "         $ds が -vector と -no-vector で異なる -> NG"
				h5diff -p "$TOL" "$a/time_series_data.h5" "$b/time_series_data.h5" "$ds" "$ds" 2>&1 | head -5 || true
				status=1
			fi
		done
	fi
done

if [ $status -ne 0 ]; then
	echo "*** vector / no-vector consistency check FAILED" >&2
	exit 1
fi

echo "vector / no-vector consistency check passed (tolerance $TOL)"
