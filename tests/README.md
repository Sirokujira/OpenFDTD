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

## テストを足すとき

1. `tests/test_<name>.c` を作り、`void test_<name>(void)` を定義する。
2. `tests/test_main.c` の `tests[]` テーブルに登録する。
3. `CMakeLists.txt` の `TEST_SOURCES` と `add_test` の `foreach` に追加する。

**`sol/` に .c を足さないこと。** CPU ビルドは `file(GLOB sol/*.c)` だが
CUDA (`SOURCES2`) と MPI (`SOURCES3`) は手書きリストなので、
CPU だけ通って CUDA/MPI が静かにリンクエラーになる
(`.claude/rules/build-targets.md`)。テストは必ず `tests/` に置く。

対象にできるのは **ofd.h のグローバル状態を必要としない翻訳単位**に限る。
現状 `sol/utils.c` `sol/vfeed.c` `sol/ingeometry.c` `sol/geomlines.c` が
それに当たる (いずれも `ofd.h` を include していない)。
グローバルを読むロジック (`finc` の `WaveAmp`/`WaveOmega` など) は、
`tests/test_main.c` が `MAIN` を定義して実体を用意しているので利用できる。

## 期待値の作り方

`.claude/rules/features.md` の方針に従い、**コードとは独立な出所**から取る。
このテスト群では定義式・数学的な恒等式・物理的な不変量を使っている
(実装式をそのまま写した期待値は回帰テストにしかならない)。

例:
- `d_sin(z)` は `sin x cosh y + i cos x sinh y` と比較する
- `finc` の `dfi` は `fi` の中心差分と比較する
- `vfeed` のピークは正規化の定義から `t = tw + td + (tw/4)/√2` で 1
- `calcdft` は線形性と、周波数 0 で総和になることで検証する
