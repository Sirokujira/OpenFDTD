---
paths:
  - "CMakeLists.txt"
  - "sol/**/*.c"
  - "cuda/**"
  - "cuda_mpi/**"
  - "mpi/**/*.c"
  - "post/**/*.c"
---

# 3 実装 (CPU / MPI / CUDA) とビルドリストの整合

**このリポジトリで最も踏みやすい落とし穴。CI では検出できない。**

## ソースの選ばれ方が実装ごとに違う

| ターゲット | ソースの指定 | 影響 |
|---|---|---|
| `ofd` (CPU+OpenMP) | `file(GLOB sol/*.c)` | 新規ファイルが**自動で入る** |
| `ofd_cuda` | `set(SOURCES2 ...)` 手書き | 新規ファイルは**入らない** |
| `ofd_mpi` | `set(SOURCES3 ...)` 手書き + `mpi/*.c` の glob | 新規ファイルは**入らない** |
| `ofd_post` | `set(SOURCES3 ...)` 別定義 (5 ファイルのみ) + `post/*.c` | 同上 |

CI (`.github/workflows/ci.yml`) は 3 OS とも `WITH_CUDA=OFF -DWITH_MPI=OFF`
でしかビルドしない。つまり **`sol/` に新しい .c を足すと CPU ビルドだけ通り、
CUDA/MPI ビルドは誰も気付かないままリンクエラーになる**。

`.claude/hooks/check-portability.sh` が編集のたびにこれを検査する。

## 新しい sol/*.c を足すときの選択肢

1. **3 実装すべてで使う** → `SOURCES2` (CUDA) と `SOURCES3` (MPI) の両方に追加。
2. **GPU では別実装にする** → `cuda/<同名>.cu` を実装する
   (`updateEx.c` ↔ `cuda/updateEx.cu` のように対になっている)。
3. **一部の実装だけ対応する** → README の「実装ごとの対応状況」表に明記する。

## 現状の意図的な除外

- `sol/setupNear.c` — 呼び出し元がなく、削除済みグローバル `Fnorm` を参照する。
  CPU (REMOVE_ITEM) でも CUDA (SOURCES2) でも MPI (SOURCES3) でも除外済み。
- `sol/solve.c` — MPI は `mpi/solve.c`、CUDA は `cuda/solve.cu` と対になっている。
- `sol/updateTpa.c` — CPU と MPI は対応済み (SOURCES3 に登録)。
  **CUDA は未対応** (`cuda/updateTpa.cu` が無い)。

## MPI の既知の制約

- **2 プロセス以上でデッドロックする** (TPA と無関係な既存不具合)。
  `mpi/solve.c` の `if (commRank == 0)` ブロックが `H5Gcreate` / `H5Dcreate`
  という並列 HDF5 の集団操作を rank 0 だけで呼んでいるため。
  修正するときは create/close を全ランクで呼び、`H5Dwrite` だけを rank 0 に
  限る形にする必要がある (`ofd_post` が読む HDF5 の形式を壊さないこと)。
- MPI ビルドには**並列 HDF5** が必要 (`H5Pset_fapl_mpio` を使うため)。
  Ubuntu なら `libhdf5-openmpi-dev`、`HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi`
  を指定する。直列 HDF5 だとリンクエラーになる。
- `comm_X/Y/Z` は **H 面しか交換していない**。E の近傍を読む処理を追加する
  場合は `mpi/comm_E.c` の `comm_E_X/Y/Z()` を先に呼ぶこと。

## 既存機能に手を入れるとき

`sol/solve.c` を変更したら `mpi/solve.c` と `cuda/solve.cu` にも同じ修正が
要るか必ず確認する (同名の別実装が並んでいる)。
