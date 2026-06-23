# -*- coding: utf-8 -*-
"""
sol.material.factor の物性値係数計算ロジックのテスト

期待値は C 版リファレンス (sol/setup.c) と同じ式から独立に導出している。
  - 通常媒質 (type==1):
        edenom = epsr + esgm * ETA0 * C * dt
        C1E = epsr / edenom,   C2E = 1 / edenom
  - 分散性媒質 (type==2):
        ke    = exp(-ce*dt)
        xi0   = ae*dt + (be/ce)*(1-ke)
        denom = einf + xi0
        C1E = einf / denom,    C2E = 1 / denom
  - 透磁率 (全媒質共通):
        mdenom = amur + msgm / ETA0 * C * dt
        C1H = amur / mdenom,   C2H = 1 / mdenom
"""

import math

import numpy as np
import pytest

import sol.input_data as input_data
import sol.material as material


def base_parm():
    P = input_data.const()
    P['dt'] = 1.0e-12
    P['f_dtype'] = 'f8'
    return P


def test_air_and_pec_coefficients():
    """媒質 0=空気 は (1,1)、媒質 1=PEC は (0,0)。"""
    P = base_parm()
    iMaterial = np.array([0, 1], dtype='i4')
    fMaterial = np.zeros((2, 8), dtype='f8')

    C1E, C2E, C1H, C2H = material.factor(P, iMaterial, fMaterial)

    assert (C1E[0], C2E[0], C1H[0], C2H[0]) == (1, 1, 1, 1)
    assert (C1E[1], C2E[1], C1H[1], C2H[1]) == (0, 0, 0, 0)


def test_lossless_dielectric_c1_is_one():
    """無損失誘電体 (esgm=0) では C1E=1、C2E=1/epsr という物理的不変量。"""
    P = base_parm()
    epsr = 4.0
    iMaterial = np.array([0, 1, 1], dtype='i4')
    fMaterial = np.zeros((3, 8), dtype='f8')
    fMaterial[2, 0] = epsr   # epsr
    fMaterial[2, 1] = 0.0    # esgm
    fMaterial[2, 2] = 1.0    # amur
    fMaterial[2, 3] = 0.0    # msgm

    C1E, C2E, C1H, C2H = material.factor(P, iMaterial, fMaterial)

    assert C1E[2] == pytest.approx(1.0)
    assert C2E[2] == pytest.approx(1.0 / epsr)
    # 無損失・非磁性 → 透磁率係数も (1, 1)
    assert C1H[2] == pytest.approx(1.0)
    assert C2H[2] == pytest.approx(1.0)


def test_lossy_dielectric_matches_reference_formula():
    """損失あり通常媒質の C1E/C2E/C1H/C2H が定義式と一致する。"""
    P = base_parm()
    epsr, esgm, amur, msgm = 2.5, 0.01, 1.3, 0.002
    iMaterial = np.array([0, 1, 1], dtype='i4')
    fMaterial = np.zeros((3, 8), dtype='f8')
    fMaterial[2, 0:4] = [epsr, esgm, amur, msgm]

    C1E, C2E, C1H, C2H = material.factor(P, iMaterial, fMaterial)

    edenom = epsr + (esgm * P['ETA0'] * P['C'] * P['dt'])
    mdenom = amur + (msgm / P['ETA0'] * P['C'] * P['dt'])

    assert C1E[2] == pytest.approx(epsr / edenom)
    assert C2E[2] == pytest.approx(1.0 / edenom)
    assert C1H[2] == pytest.approx(amur / mdenom)
    assert C2H[2] == pytest.approx(1.0 / mdenom)


def test_dispersive_medium_matches_reference_formula():
    """分散性媒質 (type==2) の係数が定義式と一致する。"""
    P = base_parm()
    einf, ae, be, ce = 1.5, 1.0e10, 2.0e20, 1.0e11
    amur, msgm = 1.0, 0.0
    iMaterial = np.array([0, 1, 2], dtype='i4')
    fMaterial = np.zeros((3, 8), dtype='f8')
    fMaterial[2, 2] = amur
    fMaterial[2, 3] = msgm
    fMaterial[2, 4] = einf
    fMaterial[2, 5] = ae
    fMaterial[2, 6] = be
    fMaterial[2, 7] = ce

    C1E, C2E, C1H, C2H = material.factor(P, iMaterial, fMaterial)

    ke = math.exp(-ce * P['dt'])
    xi0 = (ae * P['dt']) + (be / ce) * (1 - ke)
    denom = einf + xi0

    assert C1E[2] == pytest.approx(einf / denom)
    assert C2E[2] == pytest.approx(1.0 / denom)


def test_coefficients_are_bounded():
    """受動媒質では係数は (0, 1] の範囲に収まる(安定性の必要条件)。"""
    P = base_parm()
    iMaterial = np.array([0, 1, 1, 1], dtype='i4')
    fMaterial = np.zeros((4, 8), dtype='f8')
    fMaterial[2, 0:4] = [1.0, 0.05, 1.0, 0.0]
    fMaterial[3, 0:4] = [9.0, 0.10, 2.0, 0.01]

    C1E, C2E, C1H, C2H = material.factor(P, iMaterial, fMaterial)

    for m in (2, 3):
        assert 0.0 < C1E[m] <= 1.0
        assert 0.0 < C2E[m] <= 1.0
        assert 0.0 < C1H[m] <= 1.0
        assert 0.0 < C2H[m] <= 1.0
