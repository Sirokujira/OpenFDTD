# OpenFDTD Python ソルバー ユニットテスト

`python/sol` 以下に実装された FDTD ロジックの単体テスト。
期待値は C 版リファレンス (`sol/*.c`) および物理的な不変量から独立に導出している。

## 実行方法

```sh
cd python
pip install numpy numba pytest
python -m pytest
```

## カバーしているロジック

| テストファイル | 対象 | 検証内容 |
| --- | --- | --- |
| `test_material_factor.py` | `sol.material.factor` | 空気/PEC/通常媒質/分散性媒質の係数、安定性の範囲 |
| `test_getspan.py` | `sol.material._getSpan` | 区間探索 (格子点・境界外・クランプ) |
| `test_dft.py` | `sol.dft.calc` | DFT の参照実装一致・線形性・周波数応答 |
| `test_planewave.py` | `sol.planewave.f` | 入射平面波と時間微分、伝搬遅延、打ち切り |
| `test_feed_volt.py` | `sol.feed.volt` | 給電波形の正規化・対称性 |
| `test_update_equivalence.py` | `sol.update{Ex..Hz}` | vector 版 / no-vector 版の数学的等価性 |
| `test_update_reference.py` | `sol.updateEx/Hx` | Yee 格子ステンシルの参照実装一致 |
| `test_average.py` | `sol.average.calcA` | 平均電磁界の和、加法性・非負性 |
