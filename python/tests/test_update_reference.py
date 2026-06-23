# -*- coding: utf-8 -*-
"""
no-vector 更新カーネルを純 Python の参照実装と照合するテスト

numba カーネルが正しい Yee 格子ステンシル(隣接オフセット・符号)を
実装していることを、独立に書き下した素朴なループ実装で検証する。
参照式は C 版 (sol/updateEx.c, sol/updateHx.c) と一致している。
"""

import numpy as np

import sol.updateEx as updateEx
import sol.updateHx as updateHx

from test_update_equivalence import make_index, NMAT


def ref_Ex(Ex, Hy, Hz, imat, C1E, C2E, RYn, RZn, bounds):
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0 = bounds
    out = Ex.copy()
    for i in range(iMin, iMax):
        for j in range(jMin, jMax + 1):
            for k in range(kMin, kMax + 1):
                n = Ni * i + Nj * j + Nk * k + N0
                m = imat[n]
                out[n] = C1E[m] * Ex[n] + C2E[m] * (
                    RYn[j] * (Hz[n] - Hz[n - Nj]) - RZn[k] * (Hy[n] - Hy[n - Nk]))
    return out


def ref_Hx(Hx, Ey, Ez, imat, C1H, C2H, RYc, RZc, bounds):
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0 = bounds
    out = Hx.copy()
    for i in range(iMin, iMax + 1):
        for j in range(jMin, jMax):
            for k in range(kMin, kMax):
                n = Ni * i + Nj * j + Nk * k + N0
                m = imat[n]
                out[n] = C1H[m] * Hx[n] - C2H[m] * (
                    RYc[j] * (Ez[n + Nj] - Ez[n]) - RZc[k] * (Ey[n + Nk] - Ey[n]))
    return out


def test_Ex_no_vector_matches_reference():
    idx = make_index(3, 4, 5)
    NN, bounds = idx[-1], idx[:-1]
    rng = np.random.default_rng(101)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1E = rng.uniform(0.1, 1.0, NMAT)
    C2E = rng.uniform(0.1, 1.0, NMAT)
    Ex = rng.standard_normal(NN)
    Hy = rng.standard_normal(NN)
    Hz = rng.standard_normal(NN)
    RY = rng.uniform(0.5, 1.5, NN)
    RZ = rng.uniform(0.5, 1.5, NN)

    got = Ex.copy()
    updateEx._Ex_f_no_vector(got, Hy, Hz, imat, C1E, C2E, RY, RZ, *bounds)
    exp = ref_Ex(Ex, Hy, Hz, imat, C1E, C2E, RY, RZ, bounds)

    assert np.allclose(got, exp, rtol=1e-12, atol=1e-14)


def test_Hx_no_vector_matches_reference():
    idx = make_index(3, 4, 5)
    NN, bounds = idx[-1], idx[:-1]
    rng = np.random.default_rng(202)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1H = rng.uniform(0.1, 1.0, NMAT)
    C2H = rng.uniform(0.1, 1.0, NMAT)
    Hx = rng.standard_normal(NN)
    Ey = rng.standard_normal(NN)
    Ez = rng.standard_normal(NN)
    RY = rng.uniform(0.5, 1.5, NN)
    RZ = rng.uniform(0.5, 1.5, NN)

    got = Hx.copy()
    updateHx._Hx_f_no_vector(got, Ey, Ez, imat, C1H, C2H, RY, RZ, *bounds)
    exp = ref_Hx(Hx, Ey, Ez, imat, C1H, C2H, RY, RZ, bounds)

    assert np.allclose(got, exp, rtol=1e-12, atol=1e-14)


def test_vacuum_update_is_pure_curl():
    """空気(C1=C2=1)では Ex の更新が純粋な回転(curl H)項になる。"""
    idx = make_index(3, 3, 3)
    NN, bounds = idx[-1], idx[:-1]
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0 = bounds
    rng = np.random.default_rng(303)
    imat = np.zeros(NN, dtype='i4')      # 全セル空気
    C1E = np.array([1.0, 0.0, 0.5, 0.5])
    C2E = np.array([1.0, 0.0, 0.5, 0.5])
    Ex = rng.standard_normal(NN)
    Hy = rng.standard_normal(NN)
    Hz = rng.standard_normal(NN)
    RY = rng.uniform(0.5, 1.5, NN)
    RZ = rng.uniform(0.5, 1.5, NN)

    got = Ex.copy()
    updateEx._Ex_f_no_vector(got, Hy, Hz, imat, C1E, C2E, RY, RZ, *bounds)

    # 1 セルを抜き出して手計算と照合
    i, j, k = 1, 2, 1
    n = Ni * i + Nj * j + Nk * k + N0
    expected = Ex[n] + (RY[j] * (Hz[n] - Hz[n - Nj]) - RZ[k] * (Hy[n] - Hy[n - Nk]))
    assert got[n] == np.float64(expected) or abs(got[n] - expected) < 1e-12
