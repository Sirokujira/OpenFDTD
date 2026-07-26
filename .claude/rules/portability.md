---
paths:
  - "sol/**/*.c"
  - "mpi/**/*.c"
  - "post/**/*.c"
  - "include/**/*.h"
  - "CMakeLists.txt"
---

# 移植性の絶対規則 (Windows CI で実際に踏んだもの)

Linux/macOS では通るが Windows (MSVC) で落ちるものだけを挙げる。
`.claude/hooks/check-portability.sh` が編集のたびに (1)(2) を自動検査する。

1. **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
   フラット配列を使う。`a[i * n + j]` の形にする。
2. **OpenMP for のインデックスは事前宣言する** (MSVC C3015)。MSVC の OpenMP は
   2.0 相当で `#pragma omp parallel for` の直後に `for (int i = ...)` と
   書けない。`int i;` を前置して `for (i = ...)` にする。
   既存の 74 箇所はすべてこの形になっている。
3. **float\*/double\* の取り違え禁止**。配列の実型と読み出しポインタ型の
   不一致は Windows で 0xC0000005 クラッシュする (glibc は偶然耐えるので
   Linux では発覚しない)。過去に `calculatePowerLoss` で踏んだ。
   `float` の配列を `double *` で読むコードを書いていないか、配列の確保箇所と
   読み出し箇所の両方を見て確認する。
4. libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
   置換するとき**継続行の単独の `m` も引っかかる**ので、変更後に
   `grep -n "MATH_LIB\|[^_a-zA-Z]m)" CMakeLists.txt` で確認する。
5. MSVC フラグは CMakeLists の既存ブロックに従う
   (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`, `NOMINMAX`)。
6. **MSVC の既定スタックは 1 MB** (Linux は 8 MB)。FDTD ループには足りないので
   `ofd` ターゲットに `/STACK:16777216` を付けている。新しい実行ファイルを
   足すときも同様に必要。
7. 数学定数は `PI` / `EPS0` 等の既存マクロを使う (`M_PI` は MSVC では
   `_USE_MATH_DEFINES` 依存)。
8. `include/complex.h` は `d_complex_t` が**事前に typedef されている前提**で
   書かれている (`sol/utils.c` を参照)。単独で include しない。

## 検査

```bash
# 3 OS 共通の警告掃討
for f in sol/*.c; do gcc -std=c11 -Wall -Wextra -fopenmp -Iinclude -c "$f" -o /dev/null; done
```
