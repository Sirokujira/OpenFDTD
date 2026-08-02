# -*- coding: utf-8 -*-
"""
sol.dft.calc のテスト

calc(ntime, f, freq, dt, shift) は離散フーリエ変換
    csum = sum_{n=0}^{ntime-1} exp(-i * omega * (n + shift) * dt) * f[n]
    (omega = 2*pi*freq)
を計算する汎用関数。期待値はソース中のコメントにある参照ループから独立に算出する。
"""

import cmath
import math

import numpy as np
import pytest

import sol.dft as dft


def reference_dft(ntime, f, freq, dt, shift):
    omega = 2 * math.pi * freq
    csum = complex(0, 0)
    for n in range(ntime):
        ot = omega * (n + shift) * dt
        csum += cmath.exp(complex(0, -ot)) * f[n]
    return csum


def test_matches_reference_loop():
    rng = np.random.default_rng(12345)
    f = rng.standard_normal(64)
    dt = 2.0e-12
    freq = 3.0e9
    for shift in (0.0, -0.5, 1.0):
        got = dft.calc(len(f), f, freq, dt, shift)
        exp = reference_dft(len(f), f, freq, dt, shift)
        assert got == pytest.approx(exp, rel=1e-10, abs=1e-12)


def test_zero_frequency_is_plain_sum():
    f = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    got = dft.calc(len(f), f, 0.0, 1.0e-12, 0.0)
    assert got == pytest.approx(complex(np.sum(f), 0.0))


def test_linearity():
    rng = np.random.default_rng(7)
    a = rng.standard_normal(32)
    b = rng.standard_normal(32)
    dt, freq, shift = 1.0e-12, 5.0e9, -0.5
    ca = dft.calc(32, a, freq, dt, shift)
    cb = dft.calc(32, b, freq, dt, shift)
    cab = dft.calc(32, 2.0 * a + 3.0 * b, freq, dt, shift)
    assert cab == pytest.approx(2.0 * ca + 3.0 * cb, rel=1e-10)


def test_partial_window():
    # ntime < len(f) のとき先頭 ntime 点のみを使う
    f = np.arange(10).astype('f8')
    ntime = 4
    got = dft.calc(ntime, f, 1.0e9, 1.0e-12, 0.0)
    exp = reference_dft(ntime, f, 1.0e9, 1.0e-12, 0.0)
    assert got == pytest.approx(exp, rel=1e-10)


def test_single_tone_recovered():
    # 純粋な複素正弦波 exp(+i w0 t) を入力すると、その周波数で大きなピークが立つ
    dt = 1.0e-12
    ntime = 256
    freq0 = 4.0e9
    omega0 = 2 * math.pi * freq0
    n = np.arange(ntime)
    f = np.cos(omega0 * n * dt)
    amp_on = abs(dft.calc(ntime, f, freq0, dt, 0.0))
    amp_off = abs(dft.calc(ntime, f, freq0 * 3.0, dt, 0.0))
    assert amp_on > amp_off
