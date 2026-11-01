# C 実装のユニットテスト

`sol/` と `include/` に実装されたロジックのうち、**ソルバーを走らせずに
単体で検証できる純粋な関数**を対象にしたテスト。CTest から実行する。

既存の検証との関係:

| 層 | 何を見るか | 実体 |
|---|---|---|
| ユニットテスト (ここ) | 波形・DFT・形状判定・複素数演算の**式** | `tests/`, CTest |
| スモークテスト | ソルバーが最後まで走ること | `bin/ofd -n 2 dipole.ofd` |
| 検証テスト | 物理が解析解と合うこと | `data/sample/*_check.sh` (TPA / 熱 / Mie) |

スモーク・検証は「壊れたら結果が変わる」ことしか言えないが、
ユニットテストは**どの式が壊れたか**を指す。

## 実行方法

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc)"

ctest --test-dir build --output-on-failure   # 全件
./bin/ofd_tests                              # 直接実行 (全件)
./bin/ofd_tests finc                         # 1 件だけ
./bin/ofd_tests --list                       # テスト名の一覧
```

`-DWITH_TESTS=OFF` でテストのビルドを外せる (既定は ON)。

## カバーしている実装

| テスト名 | 対象 | 検証内容 |
|---|---|---|
| `complex` | `include/complex.h` | 四則・絶対値・偏角・sin/cos/sqrt/inv、0 除算が (0,0) を返すこと、float 版との整合 |
| `finc` | `include/finc.h` | ガウス微分パルスと CW 正弦波の両モード、伝搬遅延、`(a·t)²≥16` の打ち切り、ランプの連続性、`dfi = dt·d(fi)/dt` |
| `utils` | `sol/utils.c` | `nearest` / `getspan` / `calcdft` / `tokenize` |
| `vfeed` | `sol/vfeed.c` | 給電波形のピーク振幅 1 への正規化、奇対称性、遅延の平行移動 |
| `ingeometry` | `sol/ingeometry.c` | 直方体・楕円体・円柱・三角柱・角錐台・円錐台の内外判定、未知コードは常に外 |
| `geomlines` | `sol/geomlines.c` | 形状ごとの線分数、**mode=1 が mode=0 の見積もりを超えて書き込まないこと** (バッファオーバーラン防止) |
| `setupid` | `sol/setupId.c` | 形状→Yee 格子の焼き付け。6 成分が正しい実位置 (節点/セル中心) を参照すること (**半セルずれの検出**)、shape=1 の getspan 高速経路と汎用経路の一致、ドメイン端に接する物体の焼き付け漏れ、PEC/分散性材料の表面補正 (`setupId_surface`) |
| `outputzin` | `sol/outputZin.c`, `sol/outputChars.c` | 入力インピーダンス `Zin = DFT(V)/DFT(I)` と 電流 DFT の半ステップシフト、反射 `Ref = 20log10\|Γ\|`、入力電力 `Pin = Re(V conj(I))` と不整合損の補正 `/(1-\|Γ\|²)`、派生量 `zinDerived` (アドミタンス [mS]・VSWR)。**ログと HDF5 が同じ値を出すこと**。`outputChars()` の振り分け: frequency1 の無い給電でも Pin (遠方界の利得用) が求まること |
| `outputspara` | `sol/outputSpara.c` | S パラメータ。伝送線路の 3 点観測から進行波 V+ と 反射波 V- を分離する計算 (`S11 = V-/V+`, `Sn1 = Vn/V+`) を、無損失・有損失の両方で **合成した時系列から元の S を復元できるか**で検証。無反射で `S11=0`、全反射で `|S11|=1`、振幅倍率に対する不変性。dB・位相への変換 `dbDeg` (`sol/utils.c`。`EPS2` による下限 −240 dB。結合度と共通)。**ログと HDF5 が同じ値を出すこと** |
| `outputcoupling` | `sol/outputCoupling.c` | 結合度 `C = DFT(Vpoint)/DFT(Vfeed)`。2 周波数を重ねた信号と給電・観測点 2 本ずつで、**周波数・給電・観測点の添字の取り違え**まで検出 |
| `outputcross` | `sol/outputCross.c` | 散乱断面積。前方にだけ放射する Huygens 源で **後方 0・前方 4π [m²]** (物理)、前方の方向を単位ベクトル `−r̂` から独立に求めて対蹠点 `(180−θ, 180+φ)` を検証、平面波の正規化 `ffctr² = π` (λ=1 m) |
| `dispersion` | `sol/setupDispersion.c`, `sol/dispersionE[xyz].c` | 分散性材料 (`material = 2 einf ae be ce`、`ε(ω) = ε∞ + ae/(jω) + be/(ce+jω)`) の再帰的畳み込み。係数 `xi0`/`dxi0`/`ke` を `χ(t)` の数値積分と比較、**`ce = 0` (導電率のみ) で 0/0 にならないこと**、curl H = 0 の一様媒質での緩和を解析解 (Debye: 時定数 `ε∞/(ce εs)` で `E0 ε∞/εs` へ、導電率: `exp(−ae t/ε∞)`) と比較、格子点リストの範囲 (自成分方向だけ `< Max`) と材料ごとの係数、平面波入射で畳み込みに入射界 (各成分の Yee 位置・偏波成分) を足すこと |
| `pml` | `sol/setupPml.c`, `sol/setupPml[EH]?.c`, `sol/pml[EH]?.c` | PML 吸収境界。損失分布から σmax を取り出し、**理論反射率 `R(0) = exp(−2∫σ dx)` が入力の R0 に戻ること**、E 節点と H セル中心を深さ順に並べて狭義単調増加 (半セルずれの検出)・両側で一致・PML 内のセル幅は境界セル幅の延長、材料係数 `1/εr`・`1/μr`・`1/√(εμ)`、格子点リストの個数 (外周 PEC 壁の接線 E は除く)・重複なし・境界材料は最も近い内部の点、6 成分の分割界の更新を Maxwell の回転 (符号・後退/前進差分・双対/主セル幅) と整合損失 `g/√(εμ)` から計算した値と比較。不等間隔・`Nx≠Ny≠Nz`・非真空の材料で行う。物理での確認は `data/sample/pml_reflection_check.sh` |
| `mur` | `sol/murH.c`, `sol/setupMurH?.c` (`factorMur`) | Mur 1 次吸収境界 (既定の `abc`)。外向きの厳密解 `F(t + s/v)` を内側の点に与え続けて外側が追従すること、**誤差が ω を半分にして約 1/8 (1 ステップ分の局所誤差 O(h³)) になる収束次数**、`v dt = d` で特性線に厳密一致。面ごとにその面の接線 E だけを誘電体/磁性体にして、位相速度 `c/√(εμ)` の材料・境界セル幅 (不等間隔で面ごとに別)・内側の点・点の個数 (辺を含まない) を検査、PEC は `g = −1` |
| `pbc` | `sol/pbc[xyz].c` | 周期境界。`solve.c` と同じ順に呼んだ後、**辺・角を含むすべてのゴーストが周期の写し (`−1 → N−1`, `N → 0`) に一致**し、それ以外 (内部・法線成分・周期でない方向のゴースト) は書き換わらないこと。単一/2 方向/全方向の 7 通り |
| `feedload` | `sol/efeed.c`, `sol/eload.c`, `sol/input2.c` (`setup_feed` / `setup_load` / `getindex`) | 給電と集中定数負荷。辺を囲む H を**線電流の右ねじ (â × r̂) で与えて電流を独立に決め**、Thevenin 電源 `V = Vs − R I`・`E = −V/d`・`t = (n+1)dt`・パルス後は PEC の辺だけ保持 (3 方向 × 空隙/線上)。R/C 負荷は辺を満たす材料 `σ = (d/A)/R`・`εr = C(d/A)/ε0`。L 負荷は辺の容量 `ε0 A/d` との **LC 共振 `ω0 = 1/√(L C)`** の周期・振幅、一定電流駆動の応答 (1 ステップ目は厳密)。辺の長さ (双対セル幅) は**境界の節点でも配列内**を読むこと |
| `update` | `sol/update[EH][xyz].c` | 電磁界の更新 (給電モード/平面波モード × vector/no-vector)。各点の増分を **Maxwell の積分形** (Ampère `ΔE = c dt ∮H·dl/A`、Faraday `ΔH = −c dt ∮E·dl/A`) から求めて比較 (周回の辺の位置・向き・長さは Yee 格子の配置から組み立てる)、材料係数 `C1/C2`・`D1/D2` の掛かり方、**散乱界の式を「全電界は材料の式、入射界は真空の式」から導いて比較** (真空で入射項が消え、PEC で全電界 0)、更新範囲外を書き換えないこと。不等間隔・セルごとにばらばらの材料・斜め入射。`vector_check.sh` は 2 経路の一致を、こちらは各経路が正しいことを見る |
| `observe` | `sol/setupDft.c`, `sol/dftNear3d.c`, `sol/vpoint.c`, `sol/average.c` | 周波数領域変換と観測。正規化 `cFdft` を**ガウス微分パルスの解析的フーリエ変換**と比較 (給電・平面波)、DFT 係数の時刻 (H は `(n+1/2)dt`、E は `(n+1)dt`)、近傍界 DFT の累積 (実部/虚部・周波数オフセット・H のゴースト 1 層を含む範囲)、観測点の電圧 `V = −E d` と平面波の全電界、収束判定の平均 (一様場・市松模様で 0・一次関数でセル中心値) |
| `setup` | `sol/setup.c` (`setup()` を通しで) | 時間刻みが **3 次元 Courant 条件** `c dt = 1/√(Σ 1/d_min²)` (最小セルを各軸の末尾に置き探索範囲を検査)、パルス幅 `Tw = 1.27/f0` と f0 の決め方 (frequency1 → frequency2 → 1/(20dt))、材料係数 C1/C2・D1/D2 が**陰的損失の Ampère / Faraday 差分式を恒等式として満たす**こと (`εr ΔE + (σdt/ε0) E(n+1) = X`、分散性は ε∞ と ξ0)・導体中の減衰率 ≈ exp(−σdt/ε)、メッシュ係数 (双対/主セル幅)、入力で与えた dt/Tw を上書きしないこと |
| `input` | `sol/input_data.c` (+ `sol/input2.c`) | 入力の解釈。**キー省略時の既定値** (後方互換)、メッシュの区間等分割、周波数の点列、solver/HDF5 間隔の切り上げ、材料 (type 1/2/**旧形式**)、形状、PML と PBC の併用で Mur、平面波の f0 (frequency2 優先、無ければ frequency1)、**CRLF の行末**、コメント・end、長い title の切り詰め、入力の誤り 25 通りが 1 を返して**落ちない**こと (nout = 0、frequency2 の無い平面波を含む) |
| `tpa` | `sol/updateTpa.c` | 二光子吸収の非線形減衰。係数を**エネルギー収支から独立に導いて比較** (CW 波の 1 周期平均の減少率 = v β I dt、電界と磁界のエネルギーが等しいことから `⟨E²−E′²⟩/(2⟨E²⟩) = ½βε0c²E0²dt`、εr に依らない)、他成分の 4 点平均 (互いに違う値で取りこぼしを検査)、平面波は全電界を減衰させて散乱界に戻すこと、β 表、β=0 の点を書き換えないこと。通しの物理検証は `tpa_slab_check.sh` |
| `size` | `sol/setupSize.c` | MPI の領域分割と配列レイアウト。全ランクの担当範囲が領域を**隙間・重なり無く覆い**、各ランクが空でなく最後のランクが端数を持つこと、ランク番号と (Ipx,Ipy,Ipz) の対応、余白 (Mur 1 層 / PML L 層) を含む索引 NA() が隅で 0 / NN−1 の 1 対 1、**分割数がセル数を超えると 1 を返す**こと |
| `vector` | `sol/setup_vector.c` | vector 経路の係数展開。同じ形状で setup() を no-vector / vector の両方で通し、**全点の実効係数 (K1/K2 と C1/C2[id]) が一致**し、分散性の点のリストの数も一致すること。境界の材料は後に焼いた形状が勝つので、PEC と分散性の箱を先に置き後の誘電体と各方向で面を共有させ、PEC 優先・分散性優先 (周囲 4 つの H の各分岐)・添字 1 の境界を踏む。`vector_check.sh` は計算結果の通しの一致、こちらは食い違う点を指す |
| `outfile` | `sol/writeout.c` ↔ `post/readout.c` | ソルバーと ofd_post の間の **ofd.out の契約**。全項目に値を入れて書き、グローバルを消して ofd_post の readout で読み戻し、全項目が一致・ファイルを読み切ることを確かめる。int64_t の NGline / NSurface (ファイル上は 4 バイト) は上位ビットにごみを入れてから読み、正しく置き換わること |
| `alloc` | `sol/memallocfree?.c`, `sol/initfield.c` | 配列の確保と解放。構成 (Mur/PML × vector 有無 × 近傍界の有無) ごとに確保し、各配列を**使う側の範囲いっぱいまで書き込み**、解放後に全ポインタを NULL にする。CI の ASan ステップで、確保不足は heap-buffer-overflow、解放漏れは (どこからも指されなくなるので) LeakSanitizer で検出される。書き込んだ後に `initfield()` が時間発展の状態 (E/H、Mur、PML 分割成分、収束履歴、給電・観測点の時系列、近傍界 DFT の和、`Niter`) をすべて 0 に戻すことも確かめる |
| `monitor` | `sol/monitor.c`, `sol/cputime.c` | ofd.log の行の形 (検証スクリプトが grep する `ABC = PML`・`normal end` など)、メモリ量・出力ファイルの見積りが**実際の確保量の ±5% で正** (1 周波数 × 1300³ セルで int の桁あふれが起きないこと)。`cputime()` が秒の境目 2 回の間で約 1 秒を返し、減らず、10 ms 未満で刻むこと |
| `plot3d` | `sol/plot3dGeom.c`, `sol/ev3d.c` | 形状の 3D 表示 (plot3dgeom)。geom.ev3 を読み、形状の線分がすべて描かれ**材料 id で色が分かれる** (PEC 黒 / その他マゼンタ)、格子線が 2(N+1)×3 本、給電・インダクタ・観測点 (差分用 2 点込み) がその辺を各色で、タイトル、HTML 版の出力 |
| `hdf5` | `sol/outputHdf5.c` | time_series_data.h5 の**値** (構成は `hdf5_layout_check.sh`)。一意な値を入れた場・近傍界 DFT・損失を書いて HDF5 API で読み戻し、(i,j,k,成分[,実部/虚部]) への詰め替えと余白を含まないこと、time_H = time − Dt/2、interval の間引き、収束履歴の iter = n·nout、metadata・Title、hdf5 = 0 で何も作らないこと |
| `post` | `post/calcFar1d.c`, `post/calcFar2d.c`, `post/calcNear1d.c`, `post/calcNear2d.c`, `post/planewave.c` | ofd_post の計算。遠方界は x 方向微小電流の放射電力 1 − (sinθ cosφ)² で各面 (X/Y/Z/V/H) の角度の約束・far1dnorm・far2d の格子を確かめる。近傍界は線 (X:(y,z), Y:(z,x), Z:(x,y)) と面の (i,j,k) の割り当てと並び。平面波の入射界は中心で Ei、進行方向 λ/4 で −j 倍、垂直方向で位相不変 |
| `farfield` | `sol/farfield.c`, `sol/nearfield_c.c` | 微小ダイポールの放射パターン (電気/磁気/Huygens)、等価電磁流 `n̂×H` / `-n̂×E` の 6 成分すべて、位相因子 `exp(+jk r̂·r')` の符号、偏波分解 (直線/円/楕円と `\|E\|²=長軸²+短軸²=RHCP²+LHCP²`)、節点補間が一次関数を厳密に再現すること、面の幾何 (法線・位置・面素・PBC)、正規化係数 `farfactor` 。frequency2 の無い入力 (NFreq2 = 0) でも面の幾何が作れること |

## テストを足すとき

1. `tests/test_<name>.c` を作り、`void test_<name>(void)` を定義する。
2. `tests/test_main.c` の `tests[]` テーブルに登録する。
3. `CMakeLists.txt` の `TEST_SOURCES` と `add_test` の `foreach` に追加する。

**`sol/` に .c を足さないこと。** CPU ビルドは `file(GLOB sol/*.c)` だが
CUDA (`SOURCES2`) と MPI (`SOURCES3`) は手書きリストなので、
CPU だけ通って CUDA/MPI が静かにリンクエラーになる
(`.claude/rules/build-targets.md`)。テストは必ず `tests/` に置く。

グローバル状態を必要としない翻訳単位 (`sol/utils.c` `sol/vfeed.c`
`sol/ingeometry.c` `sol/geomlines.c` — いずれも `ofd.h` を include
していない) はそのままテストできる。グローバルを読むロジックも、
`tests/test_main.c` が `MAIN` を定義して実体を用意しているので、
テスト側でグローバルを組み立てれば対象にできる —
`finc` (`WaveAmp`/`WaveOmega` を読む) と `setupid` (格子・材料・
ID 配列一式をテスト内で構築する。`tests/test_setupid.c` の
`setup_grid()` が雛形) がその例。ただし後者は他のテストと同じ
プロセスで走るので、確保したものは後片付けし、書き換えた
グローバルは元に戻すこと。

## 期待値の作り方

`.claude/rules/features.md` の方針に従い、**コードとは独立な出所**から取る。
このテスト群では定義式・数学的な恒等式・物理的な不変量を使っている
(実装式をそのまま写した期待値は回帰テストにしかならない)。

例:
- `d_sin(z)` は `sin x cosh y + i cos x sinh y` と比較する
- `finc` の `dfi` は `fi` の中心差分と比較する
- `vfeed` のピークは正規化の定義から `t = tw + td + (tw/4)/√2` で 1
- `calcdft` は線形性と、周波数 0 で総和になることで検証する
