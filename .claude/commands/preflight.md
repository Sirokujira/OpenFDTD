---
description: push 前の一括点検 (クリーンビルド・警告掃討・移植性・3 実装の整合・回帰・TPA)
allowed-tools: Bash(cmake:*), Bash(gcc:*), Bash(./bin/ofd:*), Bash(./bin/ofd_post:*), Bash(sh data/sample/tpa_slab_check.sh:*), Bash(sh .claude/hooks/check-portability.sh:*), Bash(grep:*), Bash(ls:*), Bash(git status:*), Bash(git diff:*)
---

push 前の点検を順に実行し、結果を表 1 つにまとめて報告してください。
途中で落ちたら、そこで止めて原因を示してください。

1. **クリーンビルド (CPU)** —
   `rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF && cmake --build build`
   がエラー 0 件であること。
2. **警告掃討** — `for f in sol/*.c; do gcc -std=c11 -Wall -Wextra -fopenmp -Iinclude -c "$f" -o /dev/null; done`
   の出力を確認 (既存分より増えていないこと)。
3. **移植性・ビルド整合** — `sh .claude/hooks/check-portability.sh` が exit 0 であること。
   これは **CI が検出できない CUDA/MPI のソースリスト漏れ**を見る唯一の砦です。
4. **3 実装の整合** — `sol/solve.c` を触った場合、`mpi/solve.c` と
   `cuda/solve.cu` にも同じ修正が要らないか確認する。
5. **dipole 回帰** — `normal end` に達し、収束履歴とインピーダンス表が
   変更前と一致すること。`ofd_post` も通ること。
6. **TPA 検証** — `sh data/sample/tpa_slab_check.sh "$PWD/bin/ofd" /tmp/tpa-preflight`
   が 3 点とも OK であること。
7. **新キーの後方互換** — 入力キーを足した場合、そのキーを含む .ofd で
   `bin/ofd_post` が落ちないこと。
8. **差分の確認** — `git status --short` と `git diff --stat` で、意図しない
   ファイル (build/, bin/, *.log, *.out) が混ざっていないこと。

報告フォーマット:

| 点検 | 結果 |
|---|---|
| クリーンビルド (CPU) | OK / NG (理由) |
| ... | ... |

該当しない項目 (例: sol/solve.c を触っていない) は「該当なし」と書いてください。
すべて OK なら最後に「push 可」とだけ書いてください。コミットや push は
明示的に頼まれるまで実行しないこと。
