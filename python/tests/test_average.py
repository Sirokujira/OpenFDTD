# -*- coding: utf-8 -*-
"""
sol.average.calcA のテスト

計算領域内の電界・磁界の節点平均の絶対値和 (sume, sumh) を計算する。
純 Python の参照ループと照合し、加法性・非負性も確認する。
"""

import numpy as np

import sol.average as average

from test_update_equivalence import make_index


def ref_calcA(Ex, Ey, Ez, Hx, Hy, Hz, bounds):
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0 = bounds
    sume = 0.0
    sumh = 0.0
    for i in range(iMin, iMax):
        for j in range(jMin, jMax):
            for k in range(kMin, kMax):
                n = Ni * i + Nj * j + Nk * k + N0
                ni = Ni * (i + 1) + Nj * j + Nk * k + N0
                nj = Ni * i + Nj * (j + 1) + Nk * k + N0
                nk = Ni * i + Nj * j + Nk * (k + 1) + N0
                njk = Ni * i + Nj * (j + 1) + Nk * (k + 1) + N0
                nki = Ni * (i + 1) + Nj * j + Nk * (k + 1) + N0
                nij = Ni * (i + 1) + Nj * (j + 1) + Nk * k + N0
                sume += (abs(Ex[n] + Ex[nj] + Ex[nk] + Ex[njk])
                         + abs(Ey[n] + Ey[nk] + Ey[ni] + Ey[nki])
                         + abs(Ez[n] + Ez[ni] + Ez[nj] + Ez[nij]))
                sumh += (abs(Hx[n] + Hx[ni])
                         + abs(Hy[n] + Hy[nj])
                         + abs(Hz[n] + Hz[nk]))
    return sume, sumh


def _fields(NN, seed):
    rng = np.random.default_rng(seed)
    return [rng.standard_normal(NN) for _ in range(6)]


def test_matches_reference_loop():
    idx = make_index(4, 4, 4)
    NN, bounds = idx[-1], idx[:-1]
    Ex, Ey, Ez, Hx, Hy, Hz = _fields(NN, 11)

    sume, sumh = average.calcA(Ex, Ey, Ez, Hx, Hy, Hz, *bounds)
    rsume, rsumh = ref_calcA(Ex, Ey, Ez, Hx, Hy, Hz, bounds)

    assert sume == np.float64(rsume) or abs(sume - rsume) < 1e-9 * (1 + abs(rsume))
    assert sumh == np.float64(rsumh) or abs(sumh - rsumh) < 1e-9 * (1 + abs(rsumh))


def test_non_negative():
    idx = make_index(3, 3, 3)
    NN, bounds = idx[-1], idx[:-1]
    Ex, Ey, Ez, Hx, Hy, Hz = _fields(NN, 22)
    sume, sumh = average.calcA(Ex, Ey, Ez, Hx, Hy, Hz, *bounds)
    assert sume >= 0.0
    assert sumh >= 0.0


def test_zero_fields_give_zero():
    idx = make_index(3, 3, 3)
    NN, bounds = idx[-1], idx[:-1]
    z = [np.zeros(NN) for _ in range(6)]
    sume, sumh = average.calcA(*z, *bounds)
    assert sume == 0.0
    assert sumh == 0.0


def test_additive_over_i_split():
    # i 方向に領域を 2 分割した和が全体の和に等しい(加法性)
    idx = make_index(6, 3, 3)
    NN, bounds = idx[-1], idx[:-1]
    iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0 = bounds
    Ex, Ey, Ez, Hx, Hy, Hz = _fields(NN, 33)

    full = average.calcA(Ex, Ey, Ez, Hx, Hy, Hz, *bounds)
    isplit = 3
    b1 = (iMin, isplit, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    b2 = (isplit, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0)
    p1 = average.calcA(Ex, Ey, Ez, Hx, Hy, Hz, *b1)
    p2 = average.calcA(Ex, Ey, Ez, Hx, Hy, Hz, *b2)

    assert full[0] == np.float64(p1[0] + p2[0]) or abs(full[0] - (p1[0] + p2[0])) < 1e-9 * (1 + abs(full[0]))
    assert full[1] == np.float64(p1[1] + p2[1]) or abs(full[1] - (p1[1] + p2[1])) < 1e-9 * (1 + abs(full[1]))
