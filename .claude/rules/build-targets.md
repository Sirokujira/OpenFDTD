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
  CUDA は `cuda/updateTpa.cu` と対になっている (CUDA_SOURCES に登録、
  SOURCES2 には `sol/updateTpa.c` を入れないこと)。
  **`ofd_cuda_mpi` (CUDA+MPI) だけ未対応** — `CUDA_SOURCES2` に
  `updateTpa.cu` が無く、`MPI_SOURCES2` に `comm_E.c` も無い。

## MPI で踏んだ落とし穴 (すべて修正済み。同じ轍を踏まないこと)

- **入力データを配り忘れるとデッドロックする**。`mpi/comm.c` の
  `comm_broadcast()` は rank 0 が読んだ入力を他ランクへ配る。**新しい入力キーを
  足したらここにも足すこと**。配り忘れると rank != 0 でその機能が丸ごと無効に
  なり、機能内の集団操作 (通信・Allreduce) の呼び出しがランク間で食い違って
  デッドロックする。TPA (`NTpa` / `Tpa` / `WaveAmp` / `WaveOmega`) で実際に踏んだ。
  症状は「1 プロセスなら正常、2 プロセス以上でハング」。
- **rank 条件の中で HDF5 の集団操作を呼ばない**。並列 HDF5
  (`H5Pset_fapl_mpio`) では `H5Gcreate` / `H5Dcreate` / `H5Dclose` /
  `H5Gclose` / `H5Fclose` は全ランクが同じ順序で呼ぶ必要がある。
  以前は `if (io)` / `if (commRank == 0)` の中で呼んでいてデッドロックしていた。
  このファイルに書かれるのはすべて rank 0 の値なので、**直列ドライバで
  rank 0 だけが開く**形にして解決済み (`ofd_post` が読む構造は不変)。
  → **MPI ビルドに並列 HDF5 は不要**になった。
- `comm_X/Y/Z` は **H 面しか交換していない**。E の近傍を読む処理を追加する
  場合は `mpi/comm_E.c` の `comm_E_X/Y/Z()` を先に呼ぶこと。
- ログ出力の改行に全角の `¥n` を使わない (`\n`)。`ofd.log` が 1 行に潰れて
  検証スクリプトが値を取れなくなる。

## MPI の検証手順

```bash
apt-get install -y libhdf5-dev libopenmpi-dev openmpi-bin
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=ON
cmake --build build-mpi --target ofd_mpi -j4

# プロセス数・分割方向によらず結果が一致することを必ず確認する
for p in "1 1 1" "1 1 2" "2 1 1" "1 2 1" "2 2 2"; do
  mpirun --allow-run-as-root --oversubscribe -n $(echo $p | tr ' ' '*' | bc) \
    ./bin/ofd_mpi -p $p -n 1 tpa_slab.ofd >/dev/null
  grep -o "TPA: transmission = [0-9.]*" ofd.log | tail -1
done
```

## 既存機能に手を入れるとき

`sol/solve.c` を変更したら `mpi/solve.c` と `cuda/solve.cu` にも同じ修正が
要るか必ず確認する (同名の別実装が並んでいる)。
