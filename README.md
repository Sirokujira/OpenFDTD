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
  feed/planewave/point/abc/pbc/frequency1/2/solver/tpa/waveamp/hdf5 +
  `plot*` ポストキー)
- 出力:
  - `ofd.log` — 実行ログ (収束履歴、`=== normal end ===` で正常終了)
  - `ofd.out` — ポスト処理用バイナリ (ポストの正本)
  - `time_series_data.h5` — HDF5。表示用の時系列・分布データ
    (→ [HDF5 出力](#hdf5-出力-画面表示用))
  - `ofd_post` 実行後: `ev.ev2` / `ev.ev3` (図形)、`far1d.log` /
    `far2d.log` / `near2d.log` など

### HDF5 出力 (画面表示用)

`time_series_data.h5` は GUI ([OpenFDTD-X](https://github.com/Sirokujira/OpenFDTD-X))
での可視化を想定した出力です。書き込みは `sol/outputHdf5.c` に集約しています
(API とファイル構成の詳細は `include/ofd_hdf5.h` のコメント)。

| グループ | 内容 | 用途 |
|---|---|---|
| `/geometry` | `Xn/Yn/Zn` (節点座標), `Xc/Yc/Zc`, `Gline` (形状の線分) | 格子・構造の描画 |
| `/timeseries` | `E`/`H` {nsnap, Nx+1, Ny+1, Nz+1, 3} の**瞬時値**, `time`/`time_H`/`itime` | 時間領域アニメーション |
| `/freqdomain` | `E`/`H` {NFreq2, Nx+1, Ny+1, Nz+1, 3, 2} の複素振幅, `freq` | 振幅・位相の分布図 |
| `/loss` | `P_loss` {NFreq2, Nx+1, Ny+1, Nz+1} [W/m³] | 発熱・損失分布 |
| `/convergence` | `iter`/`E`/`H` | 収束履歴のグラフ |
| `/metadata` | `Nx..Nz`, `Dt`, `Freq1/2`, `VFeed`/`IFeed`, `Surface`, `input_impedance` 他 | 解析条件・Zin 等 |

- 電磁界配列は **(i,j,k,成分) の自然な 3 次元形状**で、PML/Mur の冗長領域を
  除いた物理領域のみ。ソルバー内部の 1 次元添字 (`NA(i,j,k)`) を表示側が
  知る必要はありません。
- 値は **float32 + gzip 圧縮 + チャンク化**。`/timeseries` は
  `H5S_UNLIMITED` の追記形式で、`Solver.nout` ステップごとに 1 枚増えます。
- 成分は Yee 格子上の生の値です。節点 (i,j,k) の `Ex` は
  (i,j,k)-(i+1,j,k) 稜線上の値なので、必要なら表示側で補間してください
  (周波数領域の節点補間は `sol/nearfield_c.c` の `NodeE_c`/`NodeH_c` が実装)。
- `E` と `H` は leapfrog により半ステップずれます。`time` が `E` の時刻、
  `time_H = time - Dt/2` が `H` の時刻です。

#### 出力の制御 (`hdf5` キー)

```
hdf5 = <output> [interval]
  output   : 0 = HDF5 を出力しない / 1 = 出力する (既定 1)
  interval : 瞬時値スナップショットの間隔 [ステップ]
             (既定 0 = solver の nout に従う)
```

省略時は従来と完全に同じ挙動です。表示用の時系列は容量が大きいので、
不要なら `hdf5 = 0`、粗くてよければ `interval` を大きく取ります
(dipole で `hdf5 = 1 200` にすると 12 枚 → 3 枚、5.2MB → 1.6MB)。
`interval` は瞬時値を出せるのが `solver` の `nout` ごとなので、
その倍数に切り上げられます。

#### 実装ごとの対応状況

| 実装 | `/timeseries` (瞬時値) | その他のグループ |
|---|---|---|
| CPU (`ofd`) | 対応 | 対応 |
| MPI (`ofd_mpi`) | 対応 | 対応 |
| CUDA (`ofd_cuda`) | 対応 | 対応 |
| CUDA+MPI (`ofd_cuda_mpi`) | 対応 | 対応 |

MPI 版は各ランクが部分領域しか持たないため、`mpi/comm.c` の
`comm_snapshot()` が時間ループ内で全域を rank 0 に集めてから書きます
(`/freqdomain` は従来どおり `comm_near3d()` の集約後)。
`comm_snapshot()` は全ランクが参加する集団操作なので、**rank 条件の中で
呼んではいけません** (デッドロックします)。

dipole サンプルで、`/timeseries` `/freqdomain` `/geometry` のいずれも
**1 / 2 / 4 プロセスの MPI 実行が CPU 版とビット一致**することを
`h5diff` で確認しています。

### 熱解析レイヤ (実験的)

`sol/solve.c` は DFT 済み電磁界から発熱密度 P_loss を計算し、
3次元熱拡散 (`updateTemperature`) で温度分布を更新して HDF5 に書き出します。

発熱密度は各セル・各成分について

```
p = (1/2) σ_e |E|² + (1/2) σ_m |H|²   [W/m³]
```

で、σ_e は導電率 (material の 3 番目の値)、σ_m は磁気導電率 (同 5 番目) です。
**材料は成分ごとの材料 ID 配列 (`iEx`..`iHz`) からセル毎に引きます**。
実行すると `ofd.log` に、セル体積で重み付けした総和が周波数ごとに出力されます:

```
Thermal: dissipated[0] = <値> (f=<周波数> Hz)
```

近傍界 DFT は入射スペクトルで正規化されている (`sol/setupDft.c`) ため、
この値は入射振幅 1 あたりの相対量です。検証は
`data/sample/thermal_slab.ofd` + `data/sample/thermal_material_check.sh` で、
σ を変えたときの 4 つの性質 (σ>0 で正 / 同一材料なら material id に依らず
一致 / σ=0 で 0 / σ 2 倍で 2 倍) を判定します。

```sh
sh data/sample/thermal_material_check.sh /path/to/bin/ofd /tmp/thermal
```

現状は次の制約があります:

- **CPU 版 (`ofd`) のみ**。`ofd_mpi` / `ofd_cuda` / `ofd_cuda_mpi` には
  熱解析レイヤ自体が存在しません
- 熱拡散係数 α と初期温度がソース内の定数 (入力キーが無い)。
  温度の境界条件も内点のみ更新の固定境界
- セル幅を平均値 `(Xn[Nx]-Xn[0])/Nx` で取るため不等間隔メッシュ非対応
- 温度配列を周波数ごとに持つ構造 (物理的には重ね合わせ後に 1 つのはず)

### 平面波励振の検証 (完全導体球 RCS)

平面波入射の散乱断面積を、完全導体球の Mie 厳密解と比較します
(`data/sample/sphere_rcs.ofd` + `data/sample/sphere_rcs_check.sh`)。

```sh
sh data/sample/sphere_rcs_check.sh /path/to/bin/ofd /tmp/rcs
```

半径 a = 0.05 m、ka = 3.0 の完全導体球に対し、Mie 級数から
後方 σ = 4.0901e-03 m²、前方 σ = 8.4797e-02 m² が期待値です。
球を直交格子で階段近似するため厳密一致はせず、実測比は Δ=λ/25 で
後方 1.13 / 前方 1.37 (Δ=λ/20 では 1.41 / 1.45 で、細かくすると Mie 解へ
近づく)。判定は 0.625〜1.6 倍と広めに取っています。

この検証を入れたのは、**平面波の波源設定が壊れると断面積が 0 や NaN になる
一方でログは `normal end` で終わる**ためです (実際に、ガウス微分パルスの幅
パラメータが 10⁴ 倍大きく入射波が事実上存在しない状態が長期間続き、
全周波数帯の平面波サンプルの断面積が恒等的に 0 になっていました)。
近傍界 DFT の正規化係数が 0 になる場合は `sol/setupDft.c` がエラーで停止します。

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

# MPI / CUDA のバイナリにも同じスクリプト・同じ期待値を掛けられる
OFD_LAUNCHER="mpirun --oversubscribe -n 2" OFD_ARGS="-p 1 1 2" \
  sh data/sample/tpa_slab_check.sh /path/to/bin/ofd_mpi /tmp/tpa-mpi
OFD_ARGS="-cpu" sh data/sample/tpa_slab_check.sh /path/to/bin/ofd_cuda /tmp/tpa-cuda

# 領域分割の取り方によらず結果が一致することの検証 (MPI 版)
sh data/sample/tpa_decomp_check.sh /path/to/bin/ofd_mpi /tmp/tpa-decomp
```

#### 実装ごとの対応状況

| 実装 | TPA | 備考 |
|---|---|---|
| CPU (`ofd`) | 対応 | CI (3 OS) で解析解 ±7% を判定 |
| MPI (`ofd_mpi`) | 対応 | 7 通りの領域分割で CPU 版と完全一致。CI (`build-mpi`) で分割不変性と解析解を判定 |
| CUDA (`ofd_cuda`) | 対応 | `-cpu` 実行で解析解 3 点合格。CI (`build-cuda`) で判定。GPU カーネルは nvcc でコンパイル検証のみ (実機 GPU 未検証) |
| CUDA+MPI (`ofd_cuda_mpi`) | 対応 | `-cpu` 実行で 5 通りの領域分割が `ofd_cuda` と完全一致 (0.646838)、解析解 3 点合格 |

MPI 版では `updateTpa` が `|E|²` の colocated 近似のために隣接セルの E 成分を
読むため、領域分割時は E のハローを交換してから適用します (`mpi/comm_E.c`)。
この FDTD の領域分割では E 更新が読むのは H の近傍だけなので、既存の
`comm_X/Y/Z` は H 面しか交換していませんでした — TPA が E ハローを必要とする
最初の処理です。

検証結果 (`data/sample/tpa_slab.ofd`、透過率は CPU 版と完全一致):

| 実行 | 透過率 |
|---|---|
| CPU (`ofd`) | 0.646848 |
| MPI `-p 1 1 1` / `1 1 2` / `1 1 4` | 0.646848 |
| MPI `-p 2 1 1` / `1 2 1` / `2 2 1` / `2 2 2` | 0.646848 |

MPI 2 プロセスでの解析解検証も CPU 版と同じ誤差 (−0.30% / −0.98% / +1.57%) で
3 点とも合格します。

CUDA 版では `waveamp` の CW 波源も同時に対応しました。`include/finc_cuda.h` は
ガウス微分パルスしか実装しておらず、CW 波源 (TPA と同時に追加された機能) が
GPU 側へ反映されていなかったためです。CW パラメータは `param_t` の
`waveAmp` / `waveOmega` 経由でデバイスへ渡します。

CUDA 版の検証は `ofd_cuda -cpu` (GPU を使わない実行モード) で行っています。
セル毎の演算は `__host__ __device__` の共通関数なので GPU カーネルと同じ式を
通りますが、**カーネル起動構成と実機 GPU での実行は未検証**です
(この環境に GPU が無いため。nvcc によるコンパイル・リンクは通っています)。

| 実行 | 透過率 |
|---|---|
| CPU (`ofd`) | 0.646848 |
| CUDA (`ofd_cuda -cpu`) | 0.646838 (差 0.0015%、`real_t`=float の精度差) |

MPI 版で TPA を使うには**並列 MPI ビルド用の HDF5 は不要**です。

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

- push / PR ごとに次の 5 ジョブを実行する
  - `build-cpu` / `build-macos` / `build-windows` — Linux (gcc) / macOS
    (AppleClang) / Windows (MSVC) で CPU ビルド + dipole サンプルの
    スモーク実行 (`normal end` 判定) + TPA スラブ検証 (解析解 ±7% 判定)。
    Linux / macOS では熱解析のセル毎材料検証、Linux ではさらに
    平面波散乱の Mie 検証 (完全導体球 RCS) を実行
  - `build-mpi` — `ofd_mpi` をビルドし、dipole の 1/2 プロセス一致、
    TPA の領域分割不変性 (5 通り)、2 プロセスでの解析解を判定
  - `build-cuda` — `ofd_cuda` と `ofd_cuda_mpi` を nvcc でビルドし、
    `-cpu` 実行で解析解と領域分割不変性を判定
    (ランナーに GPU が無いためカーネル起動構成は未検証)
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
