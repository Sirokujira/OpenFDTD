# AGENTS.md — OpenFDTD

3 次元 FDTD 電磁界ソルバー (C)。OpenFDTD-X (GUI) から QProcess で起動される
処理カーネル。**CPU (OpenMP) / MPI / CUDA / CUDA+MPI の 4 実装**を持ち、
後処理は別バイナリ `ofd_post` が同じ `.ofd` を読んで行う。

> このファイルは Claude Code 用の `CLAUDE.md` / `.claude/rules/*.md` と同じ内容を
> 単独で読めるようまとめたもの。**片方だけ直すと食い違うので、規約を変えたら
> 両方直すこと。**

## ビルドとテスト

依存は HDF5 (と OpenMP)。素の環境では `find_package(HDF5)` で落ちるので先に入れる —
Linux: `apt-get install -y libhdf5-dev`、macOS: `brew install hdf5 libomp`、
Windows: vcpkg (CI と同じ `hdf5[core,zlib]:x64-windows-static-md`)。

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# 回帰 (dipole) : 収束履歴・インピーダンス表が変更前と一致すること
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
"$OLDPWD/bin/ofd" -n 2 dipole.ofd && grep "normal end" ofd.log

# 検証 (いずれも解析解・厳密解との比較)
sh data/sample/tpa_slab_check.sh        "$PWD/bin/ofd" /tmp/tpa-check      # TPA ±7%
sh data/sample/thermal_material_check.sh "$PWD/bin/ofd" /tmp/thermal-check # 熱解析の材料参照
sh data/sample/sphere_rcs_check.sh       "$PWD/bin/ofd" /tmp/rcs-check     # 平面波 RCS vs Mie
```

MPI / CUDA を触ったら該当ビルドも通すこと (下記「4 実装」参照)。

`.ofd` を変更・追加したら **`ofd_post` が読めることも確認する**
(`bin/ofd_post -n 2 <file>.ofd`)。同じファイルを別バイナリが読む契約になっている。

## ディレクトリ構成

| パス | 役割 |
|---|---|
| `sol/` | CPU (OpenMP) 実装。入力パース・setup・FDTD 更新・出力 |
| `mpi/` | MPI 版の差分 (`comm*.c`, `setup_mpi.c`, `solve.c`)。他は `sol/` を共用 |
| `cuda/`, `cuda_mpi/` | GPU カーネル (`.cu`)。`sol/` の一部と同名で対になる |
| `post/` | 後処理 `ofd_post` (遠方界・近傍界・作図) |
| `src/` | 各実装の `main` |
| `include/` | 共通ヘッダ (`ofd.h`, `complex.h`, `ofd_prototype.h` 等) |
| `data/sample/` | サンプル `.ofd` と検証スクリプト |

## 最も踏みやすい罠 : 4 実装とビルドリストの整合

**ソースの選ばれ方が実装ごとに違う。**

| ターゲット | ソースの指定 | 影響 |
|---|---|---|
| `ofd` (CPU+OpenMP) | `file(GLOB sol/*.c)` | 新規ファイルが**自動で入る** |
| `ofd_cuda` | `set(SOURCES2 ...)` + `CUDA_SOURCES` 手書き | 新規ファイルは**入らない** |
| `ofd_mpi` | `set(SOURCES3 ...)` 手書き + `mpi/*.c` の glob | 同上 |
| `ofd_cuda_mpi` | `SOURCES2` + `MPI_SOURCES2` + `CUDA_SOURCES2` + `cuda_mpi/*.cu` の glob | 同上 |
| `ofd_post` | `SOURCES3` 別定義 + `post/*.c` | 同上 |

`sol/` に新しい `.c` を足して `SOURCES2` / `SOURCES3` / `CUDA_SOURCES2` への追加を
忘れると、**CPU ビルドだけ通って他が静かにリンクエラーになる**。

新しい `sol/*.c` を足すときの選択肢 :

1. 4 実装すべてで使う → `SOURCES2` と `SOURCES3` の両方に追加。
2. GPU では別実装にする → `cuda/<同名>.cu` を実装する
   (`updateEx.c` ↔ `cuda/updateEx.cu` のように対になっている)。CUDA 版を作った場合は
   `sol/<同名>.c` を `SOURCES2` に入れないこと (二重定義になる)。
3. 一部の実装だけ対応する → README の「実装ごとの対応状況」表に明記する。

`.claude/hooks/check-portability.sh` がこれを機械的に検査する
(Codex からも `sh .claude/hooks/check-portability.sh` で直接叩ける)。

現状の意図的な除外 :

- `sol/setupNear.c` — 呼び出し元がなく削除済みグローバル `Fnorm` を参照。全ビルドから除外。
- `sol/solve.c` — MPI は `mpi/solve.c`、CUDA は `cuda/solve.cu`、CUDA+MPI は
  `cuda_mpi/solve.cu` と対になっている。
- `mpi/comm_E.c` — E ハロー交換の CPU 実装。CUDA+MPI では E が device 側にあり得るので
  `cuda_mpi/comm_cuda_E.cu` を使う (`MPI_SOURCES2` には入れない)。

**`sol/solve.c` を変更したら `mpi/solve.c` / `cuda/solve.cu` / `cuda_mpi/solve.cu` にも
同じ修正が要るか必ず確認する。**

## MPI で踏んだ落とし穴 (すべて修正済み。同じ轍を踏まないこと)

- **入力データを配り忘れるとデッドロックする**。`mpi/comm.c` の `comm_broadcast()` は
  rank 0 が読んだ入力を他ランクへ配る。**新しい入力キーを足したらここにも足すこと**。
  配り忘れると rank ≠ 0 でその機能が丸ごと無効になり、機能内の集団操作 (通信・Allreduce)
  の呼び出しがランク間で食い違ってデッドロックする。TPA で実際に踏んだ。
  症状は「1 プロセスなら正常、2 プロセス以上でハング」。
- **rank 条件の中で HDF5 の集団操作を呼ばない**。並列 HDF5 (`H5Pset_fapl_mpio`) では
  `H5Gcreate` / `H5Dcreate` / `H5Dclose` / `H5Gclose` / `H5Fclose` は全ランクが同じ順序で
  呼ぶ必要がある。現在は**直列ドライバで rank 0 だけが開く**形にしてある
  (書かれるのはすべて rank 0 の値なので並列ドライバを使う理由が無い)。
  → **MPI / CUDA+MPI ビルドに並列 HDF5 は不要**。
- `comm_X/Y/Z` は **H 面しか交換していない**。E の近傍を読む処理を追加する場合は
  `mpi/comm_E.c` の `comm_E_X/Y/Z()` (CUDA+MPI では `cuda_mpi/comm_cuda_E.cu` の
  `comm_cuda_E_X/Y/Z()`) を先に呼ぶこと。
- ログ出力の改行に**全角の `¥n` を使わない** (`\n`)。`ofd.log` が 1 行に潰れて検証
  スクリプトが値を取れなくなる。過去に mpi/cuda_mpi の solve で大量に混入していた。

MPI の検証 :

```bash
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=ON
cmake --build build-mpi -j"$(nproc)"
sh data/sample/tpa_decomp_check.sh "$PWD/bin/ofd_mpi" /tmp/mpi-decomp   # 領域分割不変性
OFD_LAUNCHER="mpirun --oversubscribe -n 2" OFD_ARGS="-p 1 1 2" \
  sh data/sample/tpa_slab_check.sh "$PWD/bin/ofd_mpi" /tmp/mpi-tpa       # 解析解
```

## 移植性の絶対規則 (Windows/MSVC で実際に踏んだもの)

1. **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + `a[i * n + j]` のフラット配列。
2. **OpenMP for のインデックスは事前宣言する** (MSVC C3015)。`#pragma omp parallel for`
   の直後に `for (int i = ...)` と書けない。`int i;` を前置する。既存の 74 箇所すべてこの形。
3. **float\* / double\* の取り違え禁止**。配列の実型と読み出しポインタ型の不一致は
   Windows で 0xC0000005 クラッシュする (glibc は偶然耐えるので Linux では発覚しない)。
   過去に `calculatePowerLoss` で踏んだ (DFT 配列 `cEx_r` 等は **float**)。
4. libm リンクは CMake の `MATH_LIB` 変数経由 (Windows に m.lib は無い)。置換するとき
   継続行の単独の `m` も引っかかるので、変更後に
   `grep -n "MATH_LIB\|[^_a-zA-Z]m)" CMakeLists.txt` で確認する。
5. MSVC フラグは CMakeLists の既存ブロックに従う
   (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`, `NOMINMAX`)。
6. **MSVC の既定スタックは 1 MB** (Linux は 8 MB)。FDTD ループには足りないので
   `ofd` に `/STACK:16777216` を付けている。新しい実行ファイルにも同様に必要。
7. 数学定数は `PI` / `EPS0` 等の既存マクロを使う (`M_PI` は MSVC では
   `_USE_MATH_DEFINES` 依存)。
8. `include/complex.h` は `d_complex_t` が**事前に typedef されている前提**で書かれて
   いる (`sol/utils.c` を参照)。単独で include しない。

## 入力キーを足すとき

- パースは `sol/input_data.c` に追加する。
- **既定値は「キー省略時に従来動作と完全一致」**にする (後方互換)。
  `input_data()` 冒頭の初期化ブロックに既定値を書く。
- **`ofd_post` が新キーを無害に無視できることを確認する**
  (`bin/ofd_post -n 2 <新キーを含む .ofd>` を実際に通す)。
- **MPI を使う機能なら `comm_broadcast()` にも足す** (上記の落とし穴)。
- README の入力キー説明に追記する。

## 平面波の波源と近傍界 DFT の正規化

近傍界 DFT は入射スペクトルで正規化される — `sol/setupDft.c` が波源波形を DFT した
`cFdft` で割る。**この正規化係数が 0 になると DFT が 0 除算で NaN になり、散乱断面積
などが静かに 0 になる**。ログは `normal end` で終わるのでスモークテストでは検出できない。

実際に踏んだ : `setup_planewave()` のパルス幅 `Planewave.ai` が 10⁴ 倍大きく、入射波が
全サンプル点で 0 になっていた。**全周波数帯の平面波サンプルが断面積 0 のまま長期間
放置されていた**。現在は `cFdft` が 0 なら `setupDft()` がエラーで停止し、
`data/sample/sphere_rcs_check.sh` が完全導体球の Mie 厳密解と比較している。
平面波の波源設定に手を入れたら必ずこの検証を通すこと。

波源のスペクトルが `frequency2` に乗らない構成では `waveamp` で CW に切り替えられる
(`tpa_slab.ofd` / `thermal_slab.ofd` がそうしている)。

## 検証ケースを足すとき (新機能には必須)

1. `data/sample/<name>.ofd` を作り、**先頭コメントに解析解の導出を書く** —
   使った公式、代入した数値、期待値、許容誤差の根拠。
2. 判定スクリプトは POSIX sh + awk/grep/sed のみ (`tpa_slab_check.sh` が雛形)。
   `OFD_LAUNCHER` (実行ファイルの前に置くコマンド) と `OFD_ARGS` (追加オプション) を
   見るようにしておくと、**同じスクリプト・同じ期待値のまま MPI (`mpirun -n 2`) や
   CUDA (`-cpu`) のバイナリにも掛けられる**。実装ごとに判定を書き分けないこと。
3. CI (`.github/workflows/ci.yml`) の関係するジョブすべてにステップを足す。
   **Windows ジョブは PowerShell なので sh スクリプトをそのまま呼べず、同じ判定を
   PowerShell で書き直している** — 期待値を変えたら両方直すこと (現状 TPA の期待値が
   `.sh` と `ci.yml` の PowerShell ブロックに二重管理されている)。
4. 期待値は**コードとは独立な出所**にする (解析解・文献値・別実装との比較)。
   コード自身の出力を期待値にすると回帰テストにしかならない。

## CI

`.github/workflows/ci.yml` の 5 ジョブ :

| ジョブ | ビルドするもの | 検証 |
|---|---|---|
| `build-cpu` / `build-macos` / `build-windows` | `ofd`, `ofd_post` | dipole スモーク + TPA 解析解 (Linux/macOS は熱解析、Linux は RCS も) |
| `build-mpi` | `ofd_mpi`, `ofd_post` | dipole の 1/2 プロセス一致 + 分割不変性 + TPA 解析解 (2 プロセス) |
| `build-cuda` | `ofd_cuda`, `ofd_cuda_mpi` | `-cpu` 実行で TPA 解析解 + 分割不変性 |

ランナーに GPU は無いので CUDA 版の検証は `-cpu` (GPU を使わない実行モード) で行う。
セル毎の演算は `__host__ __device__` の共通関数なので物理の正しさは判定できるが、
**カーネル起動構成と実機 GPU 実行は未検証**。

Windows は vcpkg (`hdf5[core,zlib]:x64-windows-static-md`)。szip は libaec の 429 で
落ちるため使わない。`build-cuda` は distro の `nvidia-cuda-toolkit` と gcc の組み合わせに
依存するので `ubuntu-24.04` に固定している。タグ `v*` push で Release にバイナリ添付。

## 既知の未対応・制限

- **熱解析レイヤは CPU 版 (`sol/solve.c`) のみ**。MPI / CUDA / CUDA+MPI には存在しない。
  熱拡散係数と初期温度はソース内の定数、境界は内点のみ更新の固定境界、セル幅を平均値で
  取るため不等間隔メッシュ非対応、温度配列を周波数ごとに持つ構造。
- `ofd_post` の coupling / cross_section の HDF5 展開が未実装 (正本は `ofd.out`)。
- `include/ofd_prototype.h` に定義の無い宣言が残っている
  (`dftNear1d` / `dftNear1dX/Y/Z` / `dftNear2d` / `dftNear2dX/Y/Z` / `rgbColor`)。
- `datalib/` (`.ofd` をプログラム生成する C API、5 ファイル) は CMake から参照されて
  おらずビルド対象外・テストなし。
- 遠方界の正規化そのものは精密には未検証 (Mie 比較で桁と傾向は確認済み)。
