---
description: CPU 版をビルドして dipole 回帰と TPA 検証を走らせる
allowed-tools: Bash(cmake:*), Bash(./bin/ofd:*), Bash(./bin/ofd_post:*), Bash(sh data/sample/tpa_slab_check.sh:*), Bash(grep:*), Bash(ls:*)
---

CI と同じ範囲をローカルで再現し、結果を報告してください。

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF -DWITH_MPI=OFF
cmake --build build -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"

# 回帰 (dipole) : ログの収束履歴・インピーダンス表が変更前と一致すること
mkdir -p /tmp/smoke && cp data/sample/dipole.ofd /tmp/smoke/
(cd /tmp/smoke && "$OLDPWD/bin/ofd" -n 2 dipole.ofd && grep -q "normal end" ofd.log && "$OLDPWD/bin/ofd_post" -n 2 dipole.ofd)

# TPA 検証 (解析解 T = 1/(1 + β I0 L) と ±7%)
sh data/sample/tpa_slab_check.sh "$PWD/bin/ofd" /tmp/tpa-check
```

報告は簡潔に:

- ビルド成否 (警告があればファイル名・行・種別のみ)
- dipole が `normal end` に達し `ofd.out` が非空か
- TPA の 3 点それぞれの誤差と OK/NG
- NG があれば原因の当たりを 1〜2 行で

**注意**: これは CPU 版のみです。CUDA/MPI ビルドは CI でも検証されないので、
`sol/` にファイルを足した場合は `.claude/rules/build-targets.md` に従って
ソースリストへの登録を別途確認してください (フックが自動検査します)。

長い出力をそのまま貼らないこと。
