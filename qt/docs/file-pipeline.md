# OpenFDTD-X — ファイル処理パイプライン設計

> 本家 OpenFDTD の `.ofd` テキスト形式と CLI バイナリ群 (`ofd` / `ofd_mpi` / `ofd_cuda` / `ofd_cuda_mpi`) を**そのまま再利用**しつつ、光・音響・水中・tidy3d 拡張ドメインを上位ラッパーで扱う構造。

## 1. データフロー全体図

```
                ┌─────────────────────────────────────────┐
                │   OpenFDTD-X GUI (Qt6 Widgets, C++)    │
                │   ─────────────────────────────────    │
                │   Project (in-memory data model)        │
                └─────────────────────────────────────────┘
                              │
                              ▼ OfdIO::saveToOfd()
       ┌──────────────────────┴───────────────────────┐
       ▼                                              ▼
  ┌─────────────┐                            ┌──────────────────┐
  │ patch.ofd   │ (OpenFDTD 互換テキスト)    │ patch.ofdx       │ (JSON サイドカー)
  │ UTF-8, KB   │                            │ 拡張ドメインメタ │
  │             │                            │ - 光: BPF/Ring   │
  │             │                            │ - 音響: RT60等   │
  │             │                            │ - 水中: SSP      │
  │             │                            │ - Raycast 設定   │
  └─────────────┘                            └──────────────────┘
       │
       │ QProcess
       ▼
  ┌──────────────────────────────────────────────────────┐
  │  Runner — 実行エンジン (4種から動的選択)             │
  │  ofd / ofd_mpi / ofd_cuda / ofd_cuda_mpi            │
  │                                                      │
  │   [分割モード] $ ofd -solver patch.ofd               │
  │       → patch.out (中間バイナリ, GB) + ofd.log       │
  │   [後処理]    $ ofd -post  patch.ofd                 │
  │       → ev.ev2, ev.ev3 (MB)                          │
  │   [一括モード] $ ofd patch.ofd                       │
  │       → ev.ev2, ev.ev3 直接生成                      │
  └──────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────┐  ┌─────────────┐  ┌──────────────┐
  │ ev.ev2      │  │ ev.ev3      │  │ patch.h5     │ (拡張)
  │ 2D図形      │  │ 3D図形      │  │ HDF5 時系列  │
  └─────────────┘  └─────────────┘  └──────────────┘
       │                 │                 │
       ▼                 ▼                 ▼
  ┌──────────────────────────────────────────────────┐
  │ GUI 内蔵ビューワー (ev2d.exe / ev3d.exe 機能を   │
  │ Qt6 Charts + QOpenGLWidget で再実装)             │
  └──────────────────────────────────────────────────┘
```

## 2. .ofd テキスト形式（本家準拠）

最小例（ダイポールアンテナ）:
```
OpenFDTD 2 7
xmesh = -0.075 30 0.075
ymesh = -0.075 30 0.075
zmesh = -0.075 10 -0.025 11 0.025 10 0.075
material = 2.0 0.0 1.0 0.0    # εr σ μr σm
geometry = 1 1 0 0 0 0 -0.025 0.025   # mat shape coord...
feed = Z 0 0 0 1 0 50
frequency1 = 2e9 3e9 10       # 給電点周波数 (Hz)
frequency2 = 3e9 3e9 0        # 遠方界周波数
solver = 1000 100 1e-3        # 最大反復 出力間隔 収束判定
abc = 1 5 2 1e-5              # PML層数 次数 反射係数
end
```

特徴:
- 1行1コマンド、`key = value...`
- UTF-8（BOM無し）、複数行可
- `material` / `geometry` / `feed` / `point` は複数行繰り返し
- 最終行は `end`

## 3. .ofdx 拡張サイドカー（JSON）

`.ofd` は本家互換のまま残し、拡張情報は別ファイル `.ofdx` (JSON) に分離。本家 OpenFDTD は `.ofdx` を無視するため**完全に下位互換**。

```json
{
  "schemaVersion": "1.0",
  "domain": "optical",
  "linkedOfd": "patch.ofd",

  "optical": {
    "mode": "bpf",
    "wavelength": { "min_nm": 1500, "max_nm": 1600, "div": 201 },
    "bpf": { "band_nm": [1540, 1560], "Q": 10000, "IL_dB": 0.5 },
    "dispersion": [
      { "matId": 5, "model": "drude", "wp": 1.37e16, "gamma": 1.07e14 }
    ],
    "raycast": null
  },

  "acoustic": null,
  "underwater": null,

  "raycast": {
    "enabled": false,
    "nRays": 1000000,
    "maxBounces": 12,
    "minEnergy_dB": -60,
    "sampling": "qmc",
    "fresnel": true,
    "dispersion": true
  },

  "tidy3d": {
    "exportTarget": "td.Simulation",
    "resolution": "medium",
    "autoPml": true
  },

  "voxelization": {
    "delta_m": 5.0e-4,
    "surfaceMethod": "conformal",
    "PVF": true,
    "octree": { "enabled": false, "maxLevel": 3 }
  },

  "import3D": [
    { "file": "antenna_housing.stl", "scale": 1.0, "matId": 2,
      "offset_mm": [0,0,0], "rot_deg": [0,0,0] }
  ]
}
```

## 4. 出力ファイル一覧

| 拡張子 | 内容 | サイズ | 用途 |
|---|---|---|---|
| `.ofd`  | 入力データ (本家互換, UTF-8テキスト) | KB | エディタ編集可、git管理可 |
| `.ofdx` | 拡張データ (JSON サイドカー) | KB | OpenFDTD-X 専用 |
| `.out`  | 計算結果中間バイナリ (本家) | GB | 分割モード時のみ、ポスト用 |
| `.log`  | 標準出力ログ | KB-MB | 進捗確認 |
| `.ev2`  | 2D図形データ (本家独自) | MB | ev2d.exe 用 |
| `.ev3`  | 3D図形データ (本家独自) | MB | ev3d.exe 用 |
| `.h5`   | HDF5 時系列・場分布 (拡張) | MB-GB | Python/Matlab連携 |
| `.s2p`  | Touchstone Sパラメータ (拡張) | KB | RF/光回路設計 |
| `.py`   | tidy3d Python 出力 | KB | クラウド送信 |
| `.wav`  | オーラリゼーション音声 (音響) | MB | バイノーラル試聴 |
| `.htm`  | HTML図形 (本家、引数 `-ev` 省略時) | MB | ブラウザで開ける |

## 5. 実装ファイルの責務マッピング

| ファイル | 担当 |
|---|---|
| `core/Project.h`           | データモデル (in-memory) |
| `io/OfdIO.{h,cpp}`         | `.ofd` 読み書き (本家互換パーサ) |
| `io/OfdxIO.{h,cpp}`        | `.ofdx` JSON 読み書き |
| `io/H5Writer.{h,cpp}`      | HDF5 時系列エクスポート |
| `io/Touchstone.{h,cpp}`    | `.s2p` 形式 |
| `io/Tidy3dExporter.{h,cpp}`| tidy3d `.py` 生成 |
| `io/StlImporter.{h,cpp}`   | STL/OBJ/STEP → mesh |
| `io/Voxelizer.{h,cpp}`     | mesh → Yee格子 (conformal/PVF) |
| `io/EvReader.{h,cpp}`      | `.ev2/.ev3` パーサ (Qt Charts へ転送) |
| `kernel/Runner.{h,cpp}`    | QProcess で本家バイナリを実行 |
| `kernel/ProcessKind.h`     | ofd/ofd_mpi/ofd_cuda/ofd_cuda_mpi 選択ロジック |

## 6. 操作モードの GUI への対応

本家の 3 つの利用パターン (講習資料 p.40, p.56-57) を Tweaks 化:

```
[ ] スパコン連携モード (CLI のみ)
[●] PC GUI モード
    Pipeline:
      (A) 一括   ofd patch.ofd     ← 1回実行
      (B) 分割   ofd -solver  →  ofd -post   ← ポスト試行錯誤可
      (C) 計算のみ ofd -solver  →  PCで -post  ← FOCUS連携
```

## 7. 既存ライブラリの再利用

本家 OpenFDTD 配布物の `datalib/` (C言語データ作成ライブラリ) は、`ofd_xsection()`, `ofd_geometry()`, `ofd_feed()` 等の関数で `.ofd` を生成します。本ラッパー GUI も同等の API を C++ 側で持ち、必要なら `datalib` を内部リンクして使う設計も可能。

```cpp
// 例: Project → .ofd 出力時、datalib を呼ぶ実装にも切替可
extern "C" {
  void ofd_init(void);
  void ofd_xsection(int n, ...);
  void ofd_geometry(int mat, int shape, ...);
  void ofd_outdata(const char *path);
}
```

これにより本家コードに完全準拠した `.ofd` 出力が保証されます。
