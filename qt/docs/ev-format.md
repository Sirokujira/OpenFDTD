# ev2/ev3 図形データ処理 — 3つの戦略

## 1. ev2/ev3 形式の正体

調査の結果、ev2/ev3 は **OpenFDTD専用ではない汎用描画コマンドライブラリ** であることが判明:

> <cite index="10-1">関数名はすべて"ev3d_"で始まります。実数はすべてdoubleです。座標単位は任意です。</cite>

呼び出しパターンは PostScript / DisplayPDF 系の典型例:

> <cite index="10-3">ev3d_init (必須) ev3d_newPage (必須) (描画関数) ev3d_newPage (描画関数) ・・・ ev3d_file (オプション) ev3d_output (必須)</cite>

`ev3d_file()` は **出力形式を切替可能**:
> <cite index="10-1">type : 出力形式、0=HTML形式、1=ev3形式 fn : ファイル名(拡張子はそれぞれ.htmと.ev3を推奨) binary : 1のときバイナリファイルを出力する</cite>

つまり ev2/ev3 = **描画コマンドのページ単位ダンプ** であり、対応する viewer (`ev2d.exe` / `ev3d_otk.exe`) はそれを replay する。

3D版は **OpenTK (C# OpenGL)** で実装されているため、ソースもC#。OpenFDTDのカーネル内で `ev3d_*()` 関数を呼んで生成 → exe が再描画、という分離構造。

## 2. 取り扱い戦略（推奨順）

### 戦略 A: HTML出力モードを使う（最推奨・即実装可）

本家には `-html` オプションがあり、`.ev2/.ev3` の代わりに `ev2d.htm` / `ev3d.htm` を出力できる:

> <cite index="3-1">引数"-html"をつけるとev.ev2,ev.ev3の代わりにev2d.htm,ev3d.htmが出力されます。これらはHTMLファイルなのでブラウザで開いてください。</cite>

Linuxでは既にこれが標準。`QWebEngineView` に埋め込めば GUI 内に統合表示できる。

**メリット:** 実装数行で完了、Linux/macOS対応、本家と完全同じ結果
**デメリット:** WebEngine ≒ 70MB依存追加、相互作用 (選択・回転) はJS実装に依存

```cpp
// 実装イメージ — 数行で完了
auto *view = new QWebEngineView();
view->load(QUrl::fromLocalFile("ev2d.htm"));
m_centerStack->addWidget(view);
```

### 戦略 B: ev2d.exe / ev3d_otk.exe を内部Process起動

最も簡単。本家viewerを別ウィンドウで起動して `.ev2/.ev3` をパス指定:

```cpp
QProcess::startDetached("ev2d.exe", { QDir(workDir).filePath("ev.ev2") });
QProcess::startDetached("ev3d_otk.exe", { QDir(workDir).filePath("ev.ev3") });
```

**メリット:** 即動作、本家機能100%
**デメリット:** Windows限定、別ウィンドウなので統合感が薄い

### 戦略 C: ev2/ev3 ネイティブパーサ + Qt再描画（将来拡張）

ev2/ev3 のバイナリ仕様を逆解析し、Qt の `QPainter` / `QOpenGLWidget` で直接描画。

**メリット:** GUI完全統合、Linux/macOS/組込で同等動作、可視化のフルカスタマイズ可
**デメリット:** 仕様解析が必要 (本家GitHub: `Sirokujira/OpenFDTD` の `post/` 内 `ev2d_*.c` / `ev3d_*.c` を参照)

## 3. 実装方針

```
┌─────────────────────────────────────────────┐
│ PlotPanel (Qt6, 中央タブ「結果表示」)       │
│ ─────────────────────────────────────────── │
│ 表示モード切替: [HTML] [Process] [Native]   │
└─────────────────────────────────────────────┘
        │
        ├─ HTML mode    → QWebEngineView (ev2d.htm/ev3d.htm)
        ├─ Process mode → QProcess::startDetached("ev2d.exe")
        └─ Native mode  → EvViewer2D / EvViewer3D
                          └─ EvReader が .ev2/.ev3 をパース
                             → QPainter / QOpenGLWidget で再描画
```

ユーザーの環境・好みに応じて切替可能にする。

## 4. ev描画コマンド (3D) の主なAPI

`ev3d_*` 関数群（資料 chap6 より）から推測される描画プリミティブ:

| カテゴリ | 関数 (推定) | Qt等価 |
|---|---|---|
| 初期化 | `ev3d_init`, `ev3d_newPage` | `QPainter::begin`, page break |
| 出力 | `ev3d_file`, `ev3d_output` | `QPainter::end` |
| 線 | `ev3d_line(x1,y1,z1,x2,y2,z2)` | `drawLine` (3D投影) |
| ポリゴン | `ev3d_polygon(...)` | `drawPolygon` |
| 三角形メッシュ | `ev3d_triangle(...)` | OpenGL `GL_TRIANGLES` |
| 等高線 | `ev3d_contour(...)` | カラーマップ + ライン |
| ベクトル場 | `ev3d_vector(...)` | 矢印描画 |
| 軸 | `ev3d_axis(...)` | XYZ軸 |
| 文字 | `ev3d_text(...)` | `drawText` |
| 色設定 | `ev3d_color(r,g,b)` | `QPen`/`QBrush` |
| カメラ | `ev3d_view(theta,phi)` | 視点変換 |

`ev2d_*` (2D) は同様で、Smith chart / 極座標 / 直交座標プロットなどに対応:

| 描画用途 | 関数（推定） |
|---|---|
| 直交座標プロット | `ev2d_plot(x[], y[], n)` |
| 対数軸 | `ev2d_logx`, `ev2d_logy` |
| 極座標 | `ev2d_polar(...)` |
| スミスチャート | `ev2d_smith(...)` |
| カラーマップ | `ev2d_contour2d(...)` |
| 凡例 | `ev2d_legend(...)` |

## 5. .ev2 / .ev3 ファイルの内部構造（推定）

OpenFDTDのソースを読まずに構造を推定:

```
┌─────────────────────────────────┐
│ Header (magic, version)         │
├─────────────────────────────────┤
│ Page 1 metadata (title, axes)   │
│ Drawing commands                │
│  ├─ COLOR r g b                 │
│  ├─ LINE x1 y1 [z1] x2 y2 [z2]  │
│  ├─ POLY n x1 y1 ... xn yn      │
│  ├─ TEXT x y "string"           │
│  └─ ...                         │
├─────────────────────────────────┤
│ Page 2 metadata                 │
│ Drawing commands                │
├─────────────────────────────────┤
│ ... (animation = many pages)    │
└─────────────────────────────────┘
```

バイナリ形式の場合は `binary=1` のとき各値が `double`/`int` で並ぶだけ。
テキスト形式 (`binary=0`) は人間が読める ASCII コマンド列。

## 6. OpenFDTD-X での対応プラン（推奨ロードマップ）

| Phase | 実装内容 | 工数目安 |
|---|---|---|
| **1. MVP** | 戦略B: `ev2d.exe` / `ev3d_otk.exe` を `QProcess::startDetached` | 0.5日 |
| **2. 統合** | 戦略A: `-html` 出力 + `QWebEngineView` 埋込 | 2日 |
| **3. ネイティブ** | 戦略C: `EvReader` パーサ + `QPainter`/`QOpenGLWidget` viewer | 1〜2週 (仕様解析含む) |
| **4. 拡張連携** | HDF5 出力との二重化 (Python/Matlab連携用) | 1週 |

Phase 1 で動くものを出し、Phase 2 で Linux/macOS 対応、Phase 3 で完全GUI統合を目指す。
