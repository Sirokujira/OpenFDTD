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
| `outputzin` | `sol/outputZin.c` | 入力インピーダンス `Zin = DFT(V)/DFT(I)` と 電流 DFT の半ステップシフト、反射 `Ref = 20log10\|Γ\|`、入力電力 `Pin = Re(V conj(I))` と不整合損の補正 `/(1-\|Γ\|²)`、派生量 `zinDerived` (アドミタンス [mS]・VSWR)。**ログと HDF5 が同じ値を出すこと** |
| `outputspara` | `sol/outputSpara.c` | S パラメータ。伝送線路の 3 点観測から進行波 V+ と 反射波 V- を分離する計算 (`S11 = V-/V+`, `Sn1 = Vn/V+`) を、無損失・有損失の両方で **合成した時系列から元の S を復元できるか**で検証。無反射で `S11=0`、全反射で `|S11|=1`、振幅倍率に対する不変性。dB・位相への変換 `dbDeg` (`sol/utils.c`。`EPS2` による下限 −240 dB。結合度と共通)。**ログと HDF5 が同じ値を出すこと** |
| `outputcoupling` | `sol/outputCoupling.c` | 結合度 `C = DFT(Vpoint)/DFT(Vfeed)`。2 周波数を重ねた信号と給電・観測点 2 本ずつで、**周波数・給電・観測点の添字の取り違え**まで検出 |
| `outputcross` | `sol/outputCross.c` | 散乱断面積。前方にだけ放射する Huygens 源で **後方 0・前方 4π [m²]** (物理)、前方の方向を単位ベクトル `−r̂` から独立に求めて対蹠点 `(180−θ, 180+φ)` を検証、平面波の正規化 `ffctr² = π` (λ=1 m) |
| `dispersion` | `sol/setupDispersion.c`, `sol/dispersionE[xyz].c` | 分散性材料 (`material = 2 einf ae be ce`、`ε(ω) = ε∞ + ae/(jω) + be/(ce+jω)`) の再帰的畳み込み。係数 `xi0`/`dxi0`/`ke` を `χ(t)` の数値積分と比較、**`ce = 0` (導電率のみ) で 0/0 にならないこと**、curl H = 0 の一様媒質での緩和を解析解 (Debye: 時定数 `ε∞/(ce εs)` で `E0 ε∞/εs` へ、導電率: `exp(−ae t/ε∞)`) と比較、格子点リストの範囲 (自成分方向だけ `< Max`) と材料ごとの係数、平面波入射で畳み込みに入射界 (各成分の Yee 位置・偏波成分) を足すこと |
| `pml` | `sol/setupPml.c`, `sol/setupPml[EH]?.c`, `sol/pml[EH]?.c` | PML 吸収境界。損失分布から σmax を取り出し、**理論反射率 `R(0) = exp(−2∫σ dx)` が入力の R0 に戻ること**、E 節点と H セル中心を深さ順に並べて狭義単調増加 (半セルずれの検出)・両側で一致・PML 内のセル幅は境界セル幅の延長、材料係数 `1/εr`・`1/μr`・`1/√(εμ)`、格子点リストの個数 (外周 PEC 壁の接線 E は除く)・重複なし・境界材料は最も近い内部の点、6 成分の分割界の更新を Maxwell の回転 (符号・後退/前進差分・双対/主セル幅) と整合損失 `g/√(εμ)` から計算した値と比較。不等間隔・`Nx≠Ny≠Nz`・非真空の材料で行う。物理での確認は `data/sample/pml_reflection_check.sh` |
| `mur` | `sol/murH.c`, `sol/setupMurH?.c` (`factorMur`) | Mur 1 次吸収境界 (既定の `abc`)。外向きの厳密解 `F(t + s/v)` を内側の点に与え続けて外側が追従すること、**誤差が ω を半分にして約 1/8 (1 ステップ分の局所誤差 O(h³)) になる収束次数**、`v dt = d` で特性線に厳密一致。面ごとにその面の接線 E だけを誘電体/磁性体にして、位相速度 `c/√(εμ)` の材料・境界セル幅 (不等間隔で面ごとに別)・内側の点・点の個数 (辺を含まない) を検査、PEC は `g = −1` |
| `pbc` | `sol/pbc[xyz].c` | 周期境界。`solve.c` と同じ順に呼んだ後、**辺・角を含むすべてのゴーストが周期の写し (`−1 → N−1`, `N → 0`) に一致**し、それ以外 (内部・法線成分・周期でない方向のゴースト) は書き換わらないこと。単一/2 方向/全方向の 7 通り |
| `farfield` | `sol/farfield.c`, `sol/nearfield_c.c` | 微小ダイポールの放射パターン (電気/磁気/Huygens)、等価電磁流 `n̂×H` / `-n̂×E` の 6 成分すべて、位相因子 `exp(+jk r̂·r')` の符号、偏波分解 (直線/円/楕円と `\|E\|²=長軸²+短軸²=RHCP²+LHCP²`)、節点補間が一次関数を厳密に再現すること、面の幾何 (法線・位置・面素・PBC)、正規化係数 `farfactor` |

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
