# -*- coding: utf-8 -*-
"""
電磁界更新カーネルの vector 版 / no-vector 版 等価性テスト

OpenFDTD の各成分更新には 2 つの実装がある:
  - no-vector 版 : セルごとに物性値番号 m=i{E,H}[n] を引き、係数 C1/C2[m] を使う
  - vector 版    : あらかじめセルごとに展開した係数 K1/K2[n] を使う

K1{comp}[n] = C1[i{comp}[n]], K2{comp}[n] = C2[i{comp}[n]] とおけば、
両者は数学的に完全に一致しなければならない (給電点モード source==0)。
この不変量により、片方だけに紛れ込んだ係数の取り違えを検出できる。
"""

import numpy as np
import pytest

import sol.updateEx as updateEx
import sol.updateEy as updateEy
import sol.updateEz as updateEz
import sol.updateHx as updateHx
import sol.updateHy as updateHy
import sol.updateHz as updateHz


def make_index(Nx, Ny, Nz):
    """単一領域(MPIなし)・Mur ABC を想定した 1 次元配列インデックス係数。"""
    lx = ly = lz = 1
    iMin = jMin = kMin = 0
    iMax, jMax, kMax = Nx, Ny, Nz
    Nk = 1
    Nj = (kMax - kMin + (2 * lz) + 1)
    Ni = (jMax - jMin + (2 * ly) + 1) * Nj
    N0 = -((iMin - lx) * Ni + (jMin - ly) * Nj + (kMin - lz) * Nk)
    NN = (iMax + lx) * Ni + (jMax + ly) * Nj + (kMax + lz) * Nk + N0 + 1
    return iMin, iMax, jMin, jMax, kMin, kMax, Ni, Nj, Nk, N0, NN


NMAT = 4


def _coeffs(rng, NN, imat):
    """物性値係数とセル展開係数を生成する。"""
    C1 = rng.uniform(0.1, 1.0, NMAT)
    C2 = rng.uniform(0.1, 1.0, NMAT)
    K1 = C1[imat].astype('f8')
    K2 = C2[imat].astype('f8')
    return C1, C2, K1, K2


def _common(seed=0, dims=(3, 4, 5)):
    Nx, Ny, Nz = dims
    idx = make_index(Nx, Ny, Nz)
    NN = idx[-1]
    bounds = idx[:-1]  # iMin..N0
    rng = np.random.default_rng(seed)
    return rng, NN, bounds


def test_Ex_vector_equals_no_vector():
    rng, NN, bounds = _common(1)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1E, C2E, K1, K2 = _coeffs(rng, NN, imat)
    Ex = rng.standard_normal(NN)
    Hy = rng.standard_normal(NN)
    Hz = rng.standard_normal(NN)
    RY = rng.uniform(0.5, 1.5, NN)
    RZ = rng.uniform(0.5, 1.5, NN)

    vec = Ex.copy()
    updateEx._Ex_f_vector(vec, Hy, Hz, K1, K2, RY, RZ, *bounds)
    nov = Ex.copy()
    updateEx._Ex_f_no_vector(nov, Hy, Hz, imat, C1E, C2E, RY, RZ, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)


def test_Ey_vector_equals_no_vector():
    rng, NN, bounds = _common(2)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1E, C2E, K1, K2 = _coeffs(rng, NN, imat)
    Ey = rng.standard_normal(NN)
    Hz = rng.standard_normal(NN)
    Hx = rng.standard_normal(NN)
    RZ = rng.uniform(0.5, 1.5, NN)
    RX = rng.uniform(0.5, 1.5, NN)

    vec = Ey.copy()
    updateEy._Ey_f_vector(vec, Hz, Hx, K1, K2, RZ, RX, *bounds)
    nov = Ey.copy()
    updateEy._Ey_f_no_vector(nov, Hz, Hx, imat, C1E, C2E, RZ, RX, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)


def test_Ez_vector_equals_no_vector():
    rng, NN, bounds = _common(3)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1E, C2E, K1, K2 = _coeffs(rng, NN, imat)
    Ez = rng.standard_normal(NN)
    Hx = rng.standard_normal(NN)
    Hy = rng.standard_normal(NN)
    RX = rng.uniform(0.5, 1.5, NN)
    RY = rng.uniform(0.5, 1.5, NN)

    vec = Ez.copy()
    updateEz._Ez_f_vector(vec, Hx, Hy, K1, K2, RX, RY, *bounds)
    nov = Ez.copy()
    updateEz._Ez_f_no_vector(nov, Hx, Hy, imat, C1E, C2E, RX, RY, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)


def test_Hx_vector_equals_no_vector():
    rng, NN, bounds = _common(4)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1H, C2H, K1, K2 = _coeffs(rng, NN, imat)
    Hx = rng.standard_normal(NN)
    Ey = rng.standard_normal(NN)
    Ez = rng.standard_normal(NN)
    RY = rng.uniform(0.5, 1.5, NN)
    RZ = rng.uniform(0.5, 1.5, NN)

    vec = Hx.copy()
    updateHx._Hx_f_vector(vec, Ey, Ez, K1, K2, RY, RZ, *bounds)
    nov = Hx.copy()
    updateHx._Hx_f_no_vector(nov, Ey, Ez, imat, C1H, C2H, RY, RZ, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)


def test_Hy_vector_equals_no_vector():
    rng, NN, bounds = _common(5)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1H, C2H, K1, K2 = _coeffs(rng, NN, imat)
    Hy = rng.standard_normal(NN)
    Ez = rng.standard_normal(NN)
    Ex = rng.standard_normal(NN)
    RZ = rng.uniform(0.5, 1.5, NN)
    RX = rng.uniform(0.5, 1.5, NN)

    vec = Hy.copy()
    updateHy._Hy_f_vector(vec, Ez, Ex, K1, K2, RZ, RX, *bounds)
    nov = Hy.copy()
    updateHy._Hy_f_no_vector(nov, Ez, Ex, imat, C1H, C2H, RZ, RX, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)


def test_Hz_vector_equals_no_vector():
    rng, NN, bounds = _common(6)
    imat = rng.integers(0, NMAT, NN).astype('i4')
    C1H, C2H, K1, K2 = _coeffs(rng, NN, imat)
    Hz = rng.standard_normal(NN)
    Ex = rng.standard_normal(NN)
    Ey = rng.standard_normal(NN)
    RX = rng.uniform(0.5, 1.5, NN)
    RY = rng.uniform(0.5, 1.5, NN)

    vec = Hz.copy()
    updateHz._Hz_f_vector(vec, Ex, Ey, K1, K2, RX, RY, *bounds)
    nov = Hz.copy()
    updateHz._Hz_f_no_vector(nov, Ex, Ey, imat, C1H, C2H, RX, RY, *bounds)

    assert np.allclose(vec, nov, rtol=1e-12, atol=1e-14)
