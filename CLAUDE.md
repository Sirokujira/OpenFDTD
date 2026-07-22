# OpenFDTD

3 次元 FDTD 電磁界ソルバー (C)。OpenFDTD-X (GUI) から QProcess で起動される
処理カーネル。CPU (OpenMP) / MPI / CUDA の 3 実装を持つ。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

# 回帰 (dipole): ログの収束履歴・インピーダンス表が変更前と一致すること
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/ofd -n 2 dipole.ofd && grep "normal end" ofd.log

# TPA 検証 (解析解 T=1/(1+βI0L) と ±7%)
bash data/sample/tpa_slab_check.sh bin/ofd /tmp/tpa-check
```

## 移植性の絶対規則 (Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
  フラット配列を使う。
- **float\*/double\* の取り違え禁止**: 配列の実型と読み出しポインタ型の
  不一致は Windows で 0xC0000005 クラッシュ (glibc は偶然耐える)。
  過去例: calculatePowerLoss。
- libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
  継続行の `m` も置換対象になるので注意。
- MSVC フラグは CMakeLists の既存ブロックに従う
  (`/utf-8`, `_USE_MATH_DEFINES`, `/STACK:16777216`)。
- 数学定数は `PI` / `EPS0` 等の既存マクロを使う。

## 機能追加の規則

- 入力キー追加は `sol/input_data.c` に、既定値は「キー省略時に従来動作と
  完全一致」になるよう初期化する (後方互換)。
- CPU 実装 (`sol/`) にだけ追加した機能は README に CUDA/MPI 対応状況を
  明記する。現状: `tpa` (二光子吸収) は CPU のみ、CUDA 未対応。
- 新機能には data/sample/ の検証ケース + CI スモーク (3 OS) を必ず付ける。
- `ofd_post` が新キーを無害に無視できることを確認する。

## CI

`.github/workflows/ci.yml`: Linux / macOS (libomp) / Windows
(MSVC + Ninja + vcpkg `hdf5[core,zlib]:x64-windows-static-md` —
szip は libaec の 429 で落ちるため使わない)。タグ `v*` push で
Release にバイナリ添付。
