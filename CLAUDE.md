# OpenFDTD

3 次元 FDTD 電磁界ソルバー (C)。OpenFDTD-X (GUI) から QProcess で起動される
処理カーネル。**CPU (OpenMP) / MPI / CUDA の 3 実装**を持ち、
後処理は別バイナリ `ofd_post` が同じ .ofd を読んで行う。

## ビルド / テスト

依存は HDF5 (と OpenMP)。素の環境では `find_package(HDF5)` で落ちるので
先に入れる — Linux: `apt-get install -y libhdf5-dev`、macOS: `brew install hdf5 libomp`、
Windows: vcpkg (CI と同じ `hdf5[core,zlib]:x64-windows-static-md`)。

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# 回帰 (dipole): ログの収束履歴・インピーダンス表が変更前と一致すること
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/ofd -n 2 dipole.ofd && grep "normal end" ofd.log

# TPA 検証 (解析解 T=1/(1+βI0L) と ±7%)
sh data/sample/tpa_slab_check.sh bin/ofd /tmp/tpa-check
```

`/check` (ビルド+回帰+TPA)、`/preflight` (push 前の一括点検)、
`/add-feature` (3 実装の整合と検証込みの機能追加) のスラッシュコマンドがある。

## ディレクトリ構成

| パス | 役割 |
|---|---|
| `sol/` | CPU (OpenMP) 実装。入力パース・setup・FDTD 更新・出力 |
| `mpi/` | MPI 版の差分 (`comm*.c`, `setup_mpi.c`, `solve.c`)。他は `sol/` を共用 |
| `cuda/`, `cuda_mpi/` | GPU カーネル (`.cu`)。`sol/` の一部と同名で対になる |
| `post/` | 後処理 `ofd_post` (遠方界・近傍界・作図) |
| `src/` | 各実装の `main` |
| `include/` | 共通ヘッダ (`ofd.h`, `complex.h`, `ofd_prototype.h` 等) |
| `data/sample/` | サンプル .ofd と検証スクリプト |

## 詳細な規則

作業対象のファイルに応じて `.claude/rules/` が自動で読み込まれる。
先に目を通しておくべきもの:

- `.claude/rules/build-targets.md` — **このリポジトリで最も踏みやすい罠**。
  CPU ビルドは `file(GLOB sol/*.c)` だが CUDA/MPI は手書きリストで、
  CI は CPU 版しかビルドしない。`sol/` に .c を足すと CUDA/MPI だけが
  静かに壊れる。編集のたびに `.claude/hooks/check-portability.sh` が検査する。
- `.claude/rules/portability.md` — MSVC で実際に踏んだ落とし穴
  (VLA 禁止 / OpenMP インデックス事前宣言 / float\*・double\* の取り違え /
  MATH_LIB / スタックサイズ)。
- `.claude/rules/features.md` — 入力キーの後方互換規則、`ofd_post` が
  新キーを無視できることの確認、検証ケースと CI (3 OS) の作り方。

## CI

`.github/workflows/ci.yml`: Linux / macOS (libomp) / Windows
(MSVC + Ninja + vcpkg `hdf5[core,zlib]:x64-windows-static-md` —
szip は libaec の 429 で落ちるため使わない)。
**3 OS とも `WITH_CUDA=OFF -DWITH_MPI=OFF`** なので、CUDA/MPI ビルドは
CI で検証されない。タグ `v*` push で Release にバイナリ添付。
