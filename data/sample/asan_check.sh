#!/bin/sh
# asan_check.sh — 全サンプルを AddressSanitizer 付きのソルバー / ofd_post で走らせる (CI 用)
#
# 配列外読み・解放後使用・リークは、結果がたまたま正しければ "normal end" で
# 終わるので、スモークテストや検証スクリプトでは見つからない。実際にこの方法で
#   - 給電・観測点・負荷を領域の端の節点に置くと配列外を読む (sol/input2.c)
#   - load の作業配列のリーク (sol/input_data.c)
#   - ofd_post の近傍界 1D 作図のリーク (post/plot2dNear1d0.c)
#   - frequency2 の無い入力で Debug ビルドが assert で止まる (sol/farfield.c)
#   - frequency1 の無い給電の入力で ofd.out の書き出しが落ちる (sol/outputChars.c)
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
	# HDF5 から読む経路 (ofd_post -h5、post/readhdf5.c)
	(cd "$d" && "$POST" -h5 -n 1 run.ofd > posth5.txt 2>&1) || true
	# vector 経路 (-vector) は係数配列 K1/K2 を別に確保するので、別に走らせる
	mkdir -p "$d/vec" && cp "$d/run.ofd" "$d/vec/"
	(cd "$d/vec" && "$OFD" -vector -n 1 run.ofd > ../solvec.txt 2>&1) || true
	for t in sol solvec post posth5; do
		if grep -q "ERROR: AddressSanitizer\|ERROR: LeakSanitizer" "$d/$t.txt"; then
			echo "*** $b ($t) :"
			grep -m1 "ERROR:" "$d/$t.txt"
			grep -m3 -E "^ +#[0-9]+ .* in " "$d/$t.txt"
			status=1
		fi
	done
	n=$((n + 1))
done

# frequency1 / frequency2 の片方が無い入力 (どちらも正当な入力。全サンプルが
# 両方を持つので dipole から消して作る)。以前は
#   - frequency2 が無いと alloc_farfield() の assert で Debug ビルドが計算の後に abort
#   - frequency1 が無いと Pin が求まらず writeout() が NULL を書こうとして SIGSEGV
for key in frequency1 frequency2; do
	d="$WORK/no_$key"
	mkdir -p "$d"
	sed -e 's/^solver = .*/solver = 20 10 0/' -e "/^$key/d" "$SRCDIR/dipole.ofd" > "$d/run.ofd"
	(cd "$d" && "$OFD" -n 1 run.ofd > sol.txt 2>&1) || true
	(cd "$d" && "$POST" -n 1 run.ofd > post.txt 2>&1) || true
	for t in sol post; do
		if grep -q "ERROR: AddressSanitizer\|ERROR: LeakSanitizer\|Assertion\|invalid file length" "$d/$t.txt" \
		   || ! grep -q "normal end" "$d/ofd.log" 2>/dev/null; then
			echo "*** no $key ($t) :"
			grep -m3 "ERROR:\|Assertion\|invalid" "$d/$t.txt"
			status=1
		fi
	done
done

if [ "$status" -ne 0 ]; then
	echo "*** ASan check FAILED" >&2
else
	echo "ASan check passed ($n samples + inputs without frequency1 / frequency2, ofd (no-vector / vector) + ofd_post (ofd.out / -h5))"
fi
exit $status
