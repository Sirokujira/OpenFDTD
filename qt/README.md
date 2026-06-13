# OpenFDTD-X — Qt6 マルチドメインGUI

OpenFDTD 本家の `.ofd` カーネルをそのまま再利用しつつ、**電磁 / 光 / 室内音響 /
水中音響** の各ドメインを統一 UI で扱う Qt6 Widgets デスクトップアプリ。
`claude.ai/design` で作成された HTML/CSS/JS モックを実コードへ起こしたもの。

光ドメインでは FDTD に加えて **RCWA (OpenRCWA / `orcwa`)** と
**BPM (OpenBPM / `obpm`)** の姉妹カーネルを切り替えて実行でき、tidy3d
(Flexcompute 社の光FDTDクラウド) への Python スクリプト書き出しにも対応する。

## ビルド

### 必要環境
- Qt 6.2+ (Widgets モジュールのみ必須)
- CMake 3.21+
- C++17 コンパイラ

```bash
cd qt
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### オプション
| オプション | 既定 | 説明 |
|---|---|---|
| `-DUSE_HDF5=ON`   | OFF | HDF5 時系列/プロジェクト出力 (`io/H5Writer`) |
| `-DUSE_LIBIGL=ON` | OFF | STL→Yee格子ボクセル化 (`io/Voxelizer`, `docs/libigl-integration.md`) |
| `-DBUILD_TESTS=ON`| ON  | `.ofd` ラウンドトリップ自己テスト (`ofdx_selftest`) |

### 実行
```bash
# 既存プロジェクトを開く
./build/openfdtd_x ../data/sample/dipole.ofd

# 起動ドメイン・言語を指定 (ja|en|both)
./build/openfdtd_x --domain optical --lang both

# 本家カーネルの場所 (なければ PATH を探索)
export OPENFDTD_HOME=/path/to/OpenFDTD   # ofd, ofd_mpi, ofd_cuda ...
export OPENRCWA_HOME=/path/to/OpenRCWA   # orcwa, orcwa_post ...
export OPENBPM_HOME=/path/to/OpenBPM     # obpm, obpm_post ...
```

### 自己テスト
`data/sample` と `data/benchmark` の全 `.ofd` をロード → シリアライズ → 再パースし、
構造が完全一致することを確認する (本家フォーマット互換性の保証)。
```bash
./build/ofdx_selftest ../data/sample      # 24 files, 1496 checks, 0 failures
./build/ofdx_selftest ../data/benchmark   # 8 files,  528 checks, 0 failures
```

## アーキテクチャ

データモデル中心。タブはすべて `Project` の View であり、編集すると
`Project::changed()` / `loaded()` シグナルでビューポート・ツリー・ステータスバーが
自動更新される。

```
src/
├── main.cpp                  アプリエントリ (CLI: file / --lang / --domain / --screenshot)
├── MainWindow.{h,cpp}        メインシェル (メニュー/ツールバー/ドック/ステータスバー)
├── DomainBar.{h,cpp}         電磁/光/室内音響/水中 切替タブ
├── RightDock.{h,cpp}         プロジェクトツリー + 実行ログ
├── I18n.{h,cpp}              日英バイリンガル翻訳テーブル
│
├── core/                     データモデル (値オブジェクト + Project)
│   ├── Project.{h,cpp}         全体モデル + CFL/メモリ推定 + 永続化
│   ├── Domain.h               物理ドメイン enum + アクセントカラー
│   ├── Material.h             material 行 (通常/分散性) + load 行
│   ├── Geometry.h             geometry 行 (本家 shape コード)
│   ├── Source.h               feed / planewave / point 行
│   ├── MeshAxis.h             xmesh/ymesh/zmesh (座標+分割)
│   └── PostOpts.h             plot* ポスト処理キー一式
│
├── io/
│   ├── OfdIO.{h,cpp}          .ofd テキスト読み書き (本家完全互換) + .ofdx JSON サイドカー
│   ├── Touchstone.{h,cpp}     S パラメータ .s1p/.s2p 出力
│   ├── Tidy3dExporter.{h,cpp} 光プロジェクト → tidy3d Python スクリプト生成
│   ├── H5Writer.{h,cpp}       HDF5 出力 (USE_HDF5 時)
│   ├── StlImporter.{h,cpp}    STL (バイナリ/ASCII) 取込
│   └── EvReader.h             ev2/ev3 ネイティブパーサ骨格
│
├── kernel/
│   └── Runner.{h,cpp}         QProcess で ofd/orcwa/obpm (+_mpi/_cuda) と *_post を起動
│
├── widgets/
│   ├── SectionBox.{h,cpp}     見出し付きグループボックス
│   ├── UnitNav.{h,cpp}        ユニット番号ナビ (◀ n/総数 ▶)
│   ├── Viewport3D.{h,cpp}     QPainter 製 3Dワイヤフレームビュー (OpenGL不要)
│   ├── PlotPanel.{h,cpp}      波源波形プレビュー + 収束履歴プロット
│   ├── LogConsole.{h,cpp}     等幅ログペイン
│   └── EvViewer.{h,cpp}       ev2/ev3 図形を3戦略 (HTML/外部exe/ネイティブ) で表示
│
└── tabs/                     左ドック各タブ (本家章立て + ドメイン拡張)
    ├── GeneralTab    全般 (solver/abc/pbc/frequency)
    ├── MeshTab       メッシュ (xmesh/ymesh/zmesh)
    ├── MaterialTab   物性値・集中定数 (material/load)
    ├── GeometryTab   物体形状 (geometry + STL取込)
    ├── SourceTab     波源・観測点 (feed/planewave/point)
    ├── Post1Tab      ポスト処理(1) 周波数特性
    ├── Post2Tab      ポスト処理(2) 遠方界・近傍界
    ├── OpticalTab    光解析 (FDTD/RCWA/BPM/FMM 切替) ※光ドメイン時のみ
    ├── AcousticTab   室内音響 (RT60/C80/STI/可聴化)  ※音響ドメイン時のみ
    ├── UnderwaterTab 水中音響 (SSP/SOFAR/ソナー)     ※水中ドメイン時のみ
    └── Tidy3dTab     tidy3d クラウド連携             ※光ドメイン時のみ
```

## 設計判断

1. **本家 `.ofd` 形式に完全互換** — `OfdIO` は `sol/input_data.c` と
   `post/post_data.c` が解釈する全キーを 1:1 で読み書きする。GUI が
   モデル化しないキーは `Project::extraLines()` に保持し、保存時にそのまま
   書き戻すので手編集ファイルがラウンドトリップで壊れない。

2. **拡張ドメインは `.ofdx` (JSON) に分離** — 光・音響・水中・tidy3d の設定は
   `.ofd` と同じ basename のサイドカー JSON に保存。本家カーネルは無視するので
   **下位互換 100%**。

3. **既存カーネルを subprocess として再利用** — `Runner` が CPU/GPU × MPI の
   4種バイナリと姉妹ソルバー (FDTD=`ofd`, RCWA=`orcwa`, BPM=`obpm`) を
   `QProcess` で起動。OpenMP は `OMP_NUM_THREADS`、MPI は `mpiexec -n N`。
   solver→post の2段実行 (一括モード) にも対応。

4. **tidy3d は光FDTD専用のクラウドバックエンド** — 独立した物理ドメインではなく
   光ドメインのタブとして配置。クラウド送信ボタンは光ドメイン選択時のみ有効。

5. **ビューポートは OpenGL 非依存** — `Viewport3D` は `QPainter` の正射影
   ワイヤフレーム。ヘッドレス/リモート環境でも動作する。

## デザインリファレンス
元の HTML/CSS/JS モック (`claude.ai/design` バンドル) と各画面が対応する。
`docs/` に `.ofd` 処理パイプライン・ev2/ev3 形式・libigl 統合の設計メモを収録。
