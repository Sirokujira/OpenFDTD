# -*- coding: utf-8 -*-
"""
sol.planewave.f のテスト

入射平面波の電界 fi とその時間微分項 dfi を計算する。
    tau = t - ((x-p9)*p6 + (y-p10)*p7 + (z-p11)*p8) / p14   # 遅延時間
    at  = p12 * tau
    ex  = exp(-at^2)  (at^2 >= 16 では 0 に打ち切り)
    fi  = at * ex * f0
    dfi = p13 * p12 * (1 - 2*at^2) * ex * f0

ここで p6..p8=伝搬方向単位ベクトル, p9..p11=基準点, p12=波形係数 a,
p13=dt, p14=光速 C。
"""

import math

import numpy as np
import pytest

import sol.planewave as planewave


def make_p(ai=1.0e10, dt=1.0e-12, c=2.99792458e8):
    p = np.zeros(15, dtype='f8')
    # 伝搬方向 +x
    p[6], p[7], p[8] = 1.0, 0.0, 0.0
    # 基準点
    p[9], p[10], p[11] = 0.0, 0.0, 0.0
    p[12] = ai
    p[13] = dt
    p[14] = c
    return p


def reference_f(x, y, z, t, f0, p):
    tau = t - ((x - p[9]) * p[6] + (y - p[10]) * p[7] + (z - p[11]) * p[8]) / p[14]
    at = p[12] * tau
    ex = math.exp(-at ** 2) if at ** 2 < 16 else 0.0
    fi = at * ex * f0
    dfi = p[13] * p[12] * (1 - 2 * at ** 2) * ex * f0
    return fi, dfi


def test_matches_reference_formula():
    p = make_p()
    f0 = 2.0
    rng = np.random.default_rng(1)
    for _ in range(20):
        x, y, z = rng.uniform(-1, 1, 3)
        t = rng.uniform(0, 5e-10)
        fi, dfi = planewave.f(x, y, z, t, f0, p)
        rfi, rdfi = reference_f(x, y, z, t, f0, p)
        assert fi == pytest.approx(rfi, rel=1e-10, abs=1e-15)
        assert dfi == pytest.approx(rdfi, rel=1e-10, abs=1e-15)


def test_zero_at_retarded_origin():
    # 遅延時間 tau = 0 では at=0 -> fi=0, ex=1 -> dfi = p13*p12*f0
    p = make_p()
    f0 = 3.0
    fi, dfi = planewave.f(0.0, 0.0, 0.0, 0.0, f0, p)
    assert fi == pytest.approx(0.0)
    assert dfi == pytest.approx(p[13] * p[12] * f0)


def test_far_tail_is_cut_off():
    # at^2 >= 16 では波形は 0 に打ち切られる
    p = make_p(ai=1.0e10)
    # at = p12 * t = 16 となる t を十分超える
    t_big = 100.0 / p[12]
    fi, dfi = planewave.f(0.0, 0.0, 0.0, t_big, 1.0, p)
    assert fi == 0.0
    assert dfi == 0.0


def test_propagation_delay_shifts_waveform():
    # +x 伝搬なので、x だけ離れた点では tau = t - x/c。
    # 同じ波形値が遅れて現れることを確認する。
    p = make_p()
    f0 = 1.0
    t0 = 2.0e-11
    x = 0.05
    fi_origin, _ = planewave.f(0.0, 0.0, 0.0, t0, f0, p)
    fi_shifted, _ = planewave.f(x, 0.0, 0.0, t0 + x / p[14], f0, p)
    assert fi_shifted == pytest.approx(fi_origin, rel=1e-9, abs=1e-15)


def test_dfi_is_dt_times_time_derivative():
    # dfi = dt * d(fi)/dt を中心差分で確認(物理的な意味の検証)
    p = make_p()
    f0 = 1.5
    x = y = z = 0.0
    t = 1.3e-11
    h = 1.0e-15
    fi_p, _ = planewave.f(x, y, z, t + h, f0, p)
    fi_m, _ = planewave.f(x, y, z, t - h, f0, p)
    num_deriv = (fi_p - fi_m) / (2 * h)
    _, dfi = planewave.f(x, y, z, t, f0, p)
    assert dfi == pytest.approx(p[13] * num_deriv, rel=1e-5)
