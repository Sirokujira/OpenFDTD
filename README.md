# OpenFDTD

3次元 FDTD (時間領域有限差分) 法による電磁界ソルバー。
[本家 OpenFDTD](http://www.e-em.co.jp/OpenFDTD/) をベースに、
HDF5 出力・熱解析レイヤ・Touchstone 出力などを拡張したフォークです。
GUI フロントエンド [OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X)
から subprocess として起動されます (テキストの `.ofd` が両者の契約)。

## 処理部の構成

| ディレクトリ | 役割 |
|---|---|
| `src/sol_Main.c` | ソルバー `ofd` のエントリ。入力読込 → セットアップ → `solve()` → 出力 |
| `src/post_Main.c` | ポストプロセッサ `ofd_post` のエントリ。`post_data()` → `readout()` → `post()` |
| `sol/` | 処理本体: Yee 格子更新 (`updateEx..updateHz`)、境界条件 (Mur/PML/PBC)、給電/平面波、近傍界 DFT (`dftNear3d`)、遠方界 (`farfield`)、S パラメータ/入力インピーダンス/結合度/散乱断面積の出力 |
| `sol/solve.c` | メイン時間ループ。平均収束判定に加え、**HDF5 時系列出力と熱解析レイヤ (実験的)** を含む |
| `post/` | 描画・出力: ev2/ev3 図形 (`ev2d`/`ev3d`)、周波数特性・遠方界・近傍界プロット、HDF5 読込 (`readhdf5`) |
| `mpi/`, `cuda/` | MPI 並列版 (`ofd_mpi`) と CUDA 版 (`ofd_cuda`) の差し替えカーネル |
| `include/` | 共有ヘッダ (`ofd.h` のグローバル格子/配列、`ofd_prototype.h`) |

### 入出力

- 入力: `.ofd` テキスト (`sol/input_data.c` が解釈。mesh/material/geometry/
  feed/planewave/point/abc/pbc/frequency1/2/solver/tpa/waveamp +
  `plot*` ポストキー)
- 出力:
  - `ofd.log` — 実行ログ (収束履歴、`=== normal end ===` で正常終了)
  - `ofd.out` — ポスト処理用バイナリ (ポストの正本)
  - `time_series_data.h5` — HDF5。`/metadata` (S パラメータ・Zin・結合度・
    Surface 等) と、出力ステップ毎の `/dataNNNNNN` (E/H/表面電流/発熱密度)
  - `ofd_post` 実行後: `ev.ev2` / `ev.ev3` (図形)、`far1d.log` /
    `far2d.log` / `near2d.log` など

### 熱解析レイヤ (実験的)

`sol/solve.c` は DFT 済み電磁界から発熱密度 P_loss を計算し、
3次元熱拡散 (`updateTemperature`) で温度分布を更新して HDF5 に書き出します。
現状は次の制約があります:

- 材料 ID が先頭材料に固定 (セル毎の材料参照は未対応)
- 出力ステップ毎の `/dataNNNNNN` は大容量 (dipole サンプルで ~50MB)。
  不要な場合は該当ブロックを無効化してください

### 二光子吸収 (TPA) 非線形材料

メタマテリアル装荷 Si 導波路の光活性化関数
(Honda, Shoji, Amemiya, Opt. Lett. **49**, 5811 (2024)。Si: β = 424 cm/GW)
を扱うための強度依存吸収 dI/dz = −β I² のモデルです (`sol/updateTpa.c`)。

入力キー (どちらも省略時は従来と完全に同一挙動):

```
tpa = <material_id> <β [cm/GW]>   # 複数行可。id は material の番号 (0=真空, 1=PEC, 2〜=ユーザー材料)
waveamp = <E0 [V/m]>              # 平面波を CW 正弦波 (振幅 E0) にする。角周波数は frequency1 の先頭値
```

物理式: E 更新後に各 E 成分へ一次減衰 `E *= 1/(1+γ)`,
`γ = Δt σ_NL/(2 ε0 εr)`, `σ_NL = (4/3) β εr ε0² c² |E|²` (瞬時値) を適用します
(→ `γ = (2/3) β ε0 c² Δt |E|²`, εr に依らない)。係数は CW 定常のサイクル平均が
dI/dz = −β I² (I = ½ n ε0 c E0²) を再現するよう調和平衡で導出しています
(導出全文は `sol/updateTpa.c` 冒頭コメント)。|E|² は Yee 格子上の
colocated 近似 (成分位置に他成分を 4 点平均) で評価します。
平面波入射 (散乱界定式化) では全電界 E_scat + E_inc に減衰を適用します。

検証 (`data/sample/tpa_slab.ofd` + `data/sample/tpa_slab_check.sh`):
真空 / λ/4 AR 層 / Si 相当 TPA スラブ (εr=12.1, L=2µm, β=2000 cm/GW 検証値) /
AR 層 / 真空 の擬似 1D 構成 (x,y は PBC)。AR 層により線形透過率 ≈ 1 かつ
スラブ内は前進波のみになるので、解析解は Fresnel/Fabry-Perot 補正なしで
`T = 1/(1 + β I0 L)`, `I0 = ½ ε0 c E0²`。実行すると `ofd.log` に
`TPA: transmission = <T> (I0=<I0> W/m^2)` が出力され、CI では振幅 3 点
(E0 = 5e7/1e8/1.5e8 V/m) で解析解との差 ±7% を判定します
(実測誤差は ±1.6% 程度。残差は数値分散による AR 層のわずかな不整合・
Mur-1st の残留反射・スラブ端 1 セルの material 判定に起因)。

```sh
sh data/sample/tpa_slab_check.sh /path/to/bin/ofd /tmp/tpa
```

## ビルド

必要環境: C99 コンパイラ / CMake 3.18+ / libhdf5

```sh
# CPU 版 (ofd, ofd_post → bin/)
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j

# オプション
#   -DWITH_MPI=ON   … ofd_mpi (要 MPI)
#   -DWITH_CUDA=ON  … ofd_cuda (要 CUDA Toolkit)
```

macOS では Homebrew の libomp/hdf5 を使用します
(`.github/workflows/ci.yml` の build-macos ジョブ参照)。

## 実行

```sh
cd /tmp && cp /path/to/data/sample/dipole.ofd .
/path/to/bin/ofd -n 4 dipole.ofd        # -n = OpenMP スレッド数
grep "normal end" ofd.log
/path/to/bin/ofd_post -n 4 dipole.ofd   # ev.ev2/ev.ev3 等を生成
```

## CI / Release

- push / PR ごとに Linux (gcc) / macOS (AppleClang) / Windows (MSVC) で
  CPU ビルド + dipole サンプルのスモーク実行 (`normal end` 判定) +
  TPA スラブ検証 (解析解 ±7% 判定)
- ビルド成果物は artifact (`ofd-linux-x64` / `ofd-macos-arm64`) に保存
- `v*` タグを push すると GitHub Release に `ofd-<platform>.tar.gz` が
  自動添付されます (OpenFDTD-X や nightly 統合テストの取得元)

## 姉妹リポジトリ

| リポジトリ | 手法 | バイナリ |
|---|---|---|
| [OpenRCWA](https://github.com/Sirokujira/OpenRCWA) | 周期構造 RCWA | `orcwa` |
| [OpenBPM](https://github.com/Sirokujira/OpenBPM) | 導波路 BPM | `obpm` |
| [OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X) | Qt6 GUI | `openfdtd_x` |

## Reference

- OpenFDTD — http://www.e-em.co.jp/OpenFDTD/
