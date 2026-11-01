#!/bin/sh
# asan_check.sh — 全サンプルを AddressSanitizer 付きのソルバー / ofd_post で走らせる (CI 用)
#
# 配列外読み・解放後使用・リークは、結果がたまたま正しければ "normal end" で
# 終わるので、スモークテストや検証スクリプトでは見つからない。実際にこの方法で
#   - 給電・観測点・負荷を領域の端の節点に置くと配列外を読む (sol/input2.c)
#   - load の作業配列のリーク (sol/input_data.c)
#   - ofd_post の近傍界 1D 作図のリーク (post/plot2dNear1d0.c)
# を見つけている。ここでは data/sample/*.ofd をすべて反復 20 回に縮めて走らせ、
# ソルバーと ofd_post のどちらかで sanitizer が何か報告したら失敗にする。
#
# 使い方 : asan_check.sh <ASan 付き ofd> <ASan 付き ofd_post> [作業ディレクトリ]
#
# ASan 付きのビルドの作り方 (出力先を bin/ と分けること。CMakeLists は
# 実行ファイルを bin/ に置くので、そのままだと通常のバイナリを上書きする):
#   cmake -B build-asan -DCMAKE_BUILD_TYPE=Debug -DWITH_CUDA=OFF -DWITH_MPI=OFF \
#     -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
#     -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address" \
#     -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$PWD/build-asan/bin
#   cmake --build build-asan --target ofd ofd_post ofd_tests

set -e

OFD="$1"
POST="$2"
WORK="${3:-.}"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OFD" ] || [ -z "$POST" ]; then
	echo "Usage: asan_check.sh <ofd> <ofd_post> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
status=0
n=0
for f in "$SRCDIR"/*.ofd; do
	b=$(basename "$f" .ofd)
	d="$WORK/$b"
	mkdir -p "$d"
	sed 's/^solver = .*/solver = 20 10 0/' "$f" > "$d/run.ofd"
	# 終了コードは見ない (入力によっては意図して止まるものがある)。
	# sanitizer の報告だけで判定する
	(cd "$d" && "$OFD" -n 1 run.ofd > sol.txt 2>&1) || true
	(cd "$d" && "$POST" -n 1 run.ofd > post.txt 2>&1) || true
	for t in sol post; do
		if grep -q "ERROR: AddressSanitizer\|ERROR: LeakSanitizer" "$d/$t.txt"; then
			echo "*** $b ($t) :"
			grep -m1 "ERROR:" "$d/$t.txt"
			grep -m3 -E "^ +#[0-9]+ .* in " "$d/$t.txt"
			status=1
		fi
	done
	n=$((n + 1))
done

if [ "$status" -ne 0 ]; then
	echo "*** ASan check FAILED" >&2
else
	echo "ASan check passed ($n samples, ofd + ofd_post)"
fi
exit $status
