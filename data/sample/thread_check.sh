#!/bin/sh
# thread_check.sh — OpenMP のスレッド数で結果が変わらないことの検証 (CI 用)
#
# 更新・境界・DFT などのセルごとの並列ループにデータ競合があると、結果が
# スレッド数や実行ごとに変わる。どの検証スクリプトもスレッド数を固定 (-n 2)
# して走らせるので、競合で値がずれても期待値の許容内なら見つからない。
# ここではサンプルを -n 1 と -n 4 で走らせ、次が一致することを確かめる:
#   - ofd.log (スレッド数の表示・日時・所要時間の行を除く。収束履歴・
#     インピーダンス表・TPA の透過率などが入る)
#   - ofd_post の *.log (遠方界・近傍界・給電・観測点の数値。近傍界 DFT と
#     給電の時系列から作るので、それらが一致すれば一致する)
#
# 2026-11 に 1 / 8 スレッドの ofd.out をバイト単位で比べ、違いは収束履歴
# (平均値 Eiter / Hiter) の区間だけだった。平均値は OpenMP の reduction で
# 足し合わせるので順序が変わり、最後の桁が揺らぐ (競合ではない)。ログには
# 小数 6 桁で出るので、その揺らぎはここでは見えない。
#
# サンプル : dipole (給電) / sphere (平面波) / debye (分散性材料) /
#   pml_reflection (PML) / tpa_slab (TPA) / dipole_couple (観測点) /
#   vector_lossy_feed (損失のある材料、-vector 経路も)
#
# 使い方 : thread_check.sh <ソルバー(絶対パス)> [作業ディレクトリ]
#
# 環境変数:
#   OFD_ARGS … ソルバーに渡す追加オプション (例: "-cpu")
#   OFD_POST … 後処理 (既定: ソルバーと同じ場所の ofd_post)

set -e

OFD="$1"
WORK="${2:-.}"
EXTRA="${OFD_ARGS:-}"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OFD" ]; then
	echo "Usage: thread_check.sh <ofd> [workdir]" >&2
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

# ofd.log からスレッド数・日時・所要時間の行を除く
filter_log() {
	grep -v -E -e 'thread=' -e '^  part-' -e '^  total ' \
	        -e '^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) ' "$1"
}

status=0
for run in dipole sphere debye pml_reflection tpa_slab dipole_couple vector_lossy_feed vector_lossy_feed:-vector; do
	b=${run%%:*}
	opt=""
	[ "$b" != "$run" ] && opt=${run#*:}
	tag="$b$opt"
	for n in 1 4; do
		d="$WORK/$tag/n$n"
		mkdir -p "$d"
		cp "$SRCDIR/$b.ofd" "$d/run.ofd"
		# $EXTRA / $opt は複数語に分解させたいので意図的にクォートしない
		(cd "$d" && "$OFD" $EXTRA $opt -n $n run.ofd > /dev/null 2>&1) || true
		grep -q "normal end" "$d/ofd.log" || { echo "*** $tag : solver failed (-n $n)" >&2; status=1; continue 2; }
		(cd "$d" && "$POST" -n 1 run.ofd > post.txt 2>&1) || true
		filter_log "$d/ofd.log" > "$d/ofd.filtered"
	done
	bad=""
	cmp -s "$WORK/$tag/n1/ofd.filtered" "$WORK/$tag/n4/ofd.filtered" || bad="$bad ofd.log"
	nlog=1
	for f in "$WORK/$tag/n1"/*.log; do
		g=$(basename "$f")
		[ "$g" = "ofd.log" ] && continue
		nlog=$((nlog + 1))
		cmp -s "$f" "$WORK/$tag/n4/$g" || bad="$bad $g"
	done
	if [ -n "$bad" ]; then
		echo "$tag : differ between -n 1 and -n 4 :$bad -> NG"
		diff "$WORK/$tag/n1/ofd.filtered" "$WORK/$tag/n4/ofd.filtered" | head -6
		status=1
	else
		echo "$tag : $nlog logs identical (-n 1 vs -n 4) -> OK"
	fi
done

if [ "$status" -ne 0 ]; then
	echo "*** thread check FAILED"
else
	echo "thread check passed"
fi
exit $status
