---
paths:
  - "sol/input_data.c"
  - "post/post_data.c"
  - "data/sample/**"
  - ".github/workflows/*.yml"
  - "README.md"
---

# 機能追加と検証の規則

## 入力キーを足すとき

- パースは `sol/input_data.c` に追加する。
- **既定値は「キー省略時に従来動作と完全一致」**になるよう初期化する
  (後方互換)。`input_data()` 冒頭の初期化ブロックに既定値を書く。
- **`ofd_post` が新キーを無害に無視できることを確認する**。`ofd_post` は
  同じ .ofd を読む別バイナリで、知らないキーで落ちてはいけない。
  `bin/ofd_post -n 2 <新キーを含む .ofd>` が通ることを実際に確かめる。
- README の入力キー説明に追記する。

## 実装を足すとき

`.claude/rules/build-targets.md` を必ず読むこと。CPU (`sol/`) にだけ追加した
機能は **README に CUDA/MPI の対応状況を明記する**。
現状: `tpa` (二光子吸収) と `waveamp` (CW 波源) は CPU / MPI / CUDA /
CUDA+MPI の 4 実装すべて対応済み。熱解析レイヤは CPU (`sol/solve.c`) のみ。

## 検証ケース (新機能には必須)

1. `data/sample/<name>.ofd` を作る。**先頭コメントに解析解の導出を書く** —
   使った公式、代入した数値、期待値、許容誤差の根拠。
   例: `data/sample/tpa_slab.ofd` は解析解 `T = 1/(1 + β I0 L)` の導出を
   コメントに持ち、`tpa_slab_check.sh` が振幅 3 点を掃引して ±7% で比較する。
2. 判定スクリプトは POSIX sh + awk/grep/sed のみで書く
   (`data/sample/tpa_slab_check.sh` が雛形)。
3. CI (`.github/workflows/ci.yml`) の関係するジョブすべてにステップを足す
   (`build-cpu` / `build-macos` / `build-windows` / `build-mpi` / `build-cuda`)。
   **Windows ジョブは PowerShell なので sh スクリプトをそのまま呼べず、
   同じ判定を PowerShell で書き直している** — 期待値を変えたら両方直すこと
   (現状 TPA の期待値が Linux/macOS の .sh と Windows の .ps1 相当ブロックに
   二重管理されている)。
   判定スクリプトは `OFD_LAUNCHER` (実行ファイルの前に置くコマンド) と
   `OFD_ARGS` (追加オプション) を見るので、**同じスクリプト・同じ期待値のまま
   MPI (`mpirun -n 2`) や CUDA (`-cpu`) のバイナリにも掛けられる**。
   実装ごとに判定を書き分けないこと。
4. 期待値は**コードとは独立な出所**にする。解析解・文献値・別実装との比較。
   コード自身の出力を期待値にすると回帰テストにしかならない。

## 回帰の基準

dipole サンプルのログ (収束履歴・インピーダンス表) が変更前と一致すること。

```bash
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/ && cd /tmp/smoke
"$OLDPWD/bin/ofd" -n 2 dipole.ofd && grep "normal end" ofd.log
```
