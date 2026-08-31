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
