# -*- coding: utf-8 -*-
"""
sol.feed.volt のテスト

給電点の電圧波形(微分ガウスパルス):
    arg = (t - tw - td) / (tw / 4)
    V   = sqrt(2) * exp(0.5) * arg * exp(-arg^2)

前置係数 sqrt(2)*exp(0.5) はピーク振幅が 1 になるよう正規化されている。
"""

import math

import numpy as np
import pytest

import sol.feed as feed


def reference_volt(t, tw, td):
    arg = (t - tw - td) / (tw / 4.0)
    return math.sqrt(2.0) * math.exp(0.5) * arg * math.exp(-arg ** 2)


def test_matches_reference_formula():
    tw, td = 1.0e-10, 2.0e-11
    rng = np.random.default_rng(99)
    for _ in range(20):
        t = rng.uniform(0, 5e-10)
        assert feed.volt(t, tw, td) == pytest.approx(reference_volt(t, tw, td), rel=1e-12)


def test_zero_at_center():
    # arg = 0 (t = tw + td) で 0
    tw, td = 1.0e-10, 2.0e-11
    assert feed.volt(tw + td, tw, td) == pytest.approx(0.0)


def test_peak_amplitude_is_normalized_to_one():
    # ピークは arg = 1/sqrt(2) すなわち t = tw + td + (tw/4)/sqrt(2) で振幅 1
    tw, td = 1.0e-10, 0.0
    t_peak = tw + td + (tw / 4.0) / math.sqrt(2.0)
    assert feed.volt(t_peak, tw, td) == pytest.approx(1.0, rel=1e-12)


def test_peak_is_global_maximum():
    tw, td = 1.0e-10, 0.0
    t_peak = tw + td + (tw / 4.0) / math.sqrt(2.0)
    peak = feed.volt(t_peak, tw, td)
    ts = np.linspace(0, 4e-10, 2001)
    vals = np.array([feed.volt(t, tw, td) for t in ts])
    assert peak >= vals.max() - 1e-9


def test_odd_symmetry_about_center():
    # arg の奇関数成分: V(center + d) = -V(center - d)
    tw, td = 1.0e-10, 1.0e-11
    center = tw + td
    d = 1.5e-11
    assert feed.volt(center + d, tw, td) == pytest.approx(-feed.volt(center - d, tw, td), rel=1e-12)
