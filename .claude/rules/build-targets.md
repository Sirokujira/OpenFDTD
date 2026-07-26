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
3. **CPU 専用にする** → README に対応状況を明記し、フックの `cpu_only`
   リストに追加する。現状の CPU 専用は `updateTpa.c` (二光子吸収) のみ。

## 現状の意図的な除外

- `sol/setupNear.c` — 呼び出し元がなく、削除済みグローバル `Fnorm` を参照する。
  CPU (REMOVE_ITEM) でも CUDA (SOURCES2) でも除外済み。
- `sol/solve.c` — MPI は `mpi/solve.c`、CUDA は `cuda/solve.cu` と対になっている。
- `sol/updateTpa.c` — TPA は CPU のみ (CUDA 未対応、README に明記)。

## 既存機能に手を入れるとき

`sol/solve.c` を変更したら `mpi/solve.c` と `cuda/solve.cu` にも同じ修正が
要るか必ず確認する (同名の別実装が並んでいる)。
