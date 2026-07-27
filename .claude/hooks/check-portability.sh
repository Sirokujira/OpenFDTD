#!/bin/sh
# check-portability.sh — 編集後に走る移植性・ビルド整合チェック (PostToolUse hook)
#
# CI が検出できない失敗だけを、grep で確実に見つかる形で検査する。
# 引数・標準入力には依存せずリポジトリ全体を走査するので 3 OS どこでも動く。
# 違反があれば stderr に出して exit 2 (Claude にフィードバックされる)。

root="$(cd "$(dirname "$0")/../.." && pwd)"
status=0

# ── (1) 新規 sol/*.c のビルドリスト漏れ ─────────────────────────────
# CPU ビルドは file(GLOB sol/*.c) だが、CUDA (SOURCES2) と MPI (SOURCES3) は
# CMakeLists の手書きリスト。CI は CPU 版しかビルドしないので、新しい
# sol/*.c を足すと CUDA/MPI ビルドだけが静かに壊れる (リンクエラー)。
#
# 判定 : CMakeLists.txt に "sol/<name>" の記載が無く、cuda/<stem>.cu も
#        mpi/<stem>.c も無く、CPU 専用として文書化された除外でもないもの。
# 一部の実装だけ対応することを文書化済みのファイル (README の対応状況表を参照)
cpu_only=""
for f in "$root"/sol/*.c; do
	[ -e "$f" ] || continue
	n=$(basename "$f")
	stem=${n%.c}
	case " $cpu_only " in *" $n "*) continue ;; esac
	grep -q "sol/$n" "$root/CMakeLists.txt" 2>/dev/null && continue
	[ -f "$root/cuda/$stem.cu" ] && continue
	[ -f "$root/mpi/$stem.c" ] && continue
	echo "[hook] sol/$n が CMakeLists.txt のどのソースリストにも入っていません。" >&2
	echo "       CPU ビルドは glob なので通りますが、CUDA (SOURCES2) と" >&2
	echo "       MPI (SOURCES3) は手書きリストなのでリンクエラーになります。" >&2
	echo "       CI は CPU 版しかビルドしないため気付けません。" >&2
	echo "       対応 : 各リストに追加する / cuda/$stem.cu を実装する /" >&2
	echo "              CPU 専用なら README に明記しこのスクリプトの cpu_only に追加する。" >&2
	status=2
done

# ── (2) MSVC の OpenMP は 2.0 相当 (error C3015) ────────────────────
# #pragma omp parallel for の直後の for 文でインデックスを宣言できない。
hit=$(grep -n -A2 '^[[:space:]]*#pragma omp parallel for' \
	"$root"/sol/*.c "$root"/mpi/*.c "$root"/post/*.c 2>/dev/null \
	| grep -E 'for[[:space:]]*\([[:space:]]*(int|long|size_t|unsigned)[[:space:]]')
if [ -n "$hit" ]; then
	echo "[hook] MSVC C3015: '#pragma omp parallel for' の直後の for 文で" >&2
	echo "       ループ変数を宣言しています。'int i;' を前置してください。" >&2
	echo "$hit" >&2
	status=2
fi

# ── (3) C99 VLA 禁止 (MSVC C2057/C2466) ─────────────────────────────
# 配列サイズが定数でも大文字マクロでもない宣言を拾う。
hit=$(grep -nE '^[[:space:]]*(const[[:space:]]+)?(char|int|long|float|double|size_t)[[:space:]]+[a-z_][a-zA-Z0-9_]*\[[a-z_][a-zA-Z0-9_]*\][[:space:]]*;' \
	"$root"/sol/*.c "$root"/mpi/*.c "$root"/post/*.c 2>/dev/null)
if [ -n "$hit" ]; then
	echo "[hook] C99 VLA の可能性 (MSVC C2057/C2466)。malloc + 明示インデックスの" >&2
	echo "       フラット配列に置き換えてください (サイズが定数式なら無視して構いません)。" >&2
	echo "$hit" >&2
	status=2
fi

exit $status
