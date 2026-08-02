# -*- coding: utf-8 -*-
"""
sol.material._getSpan のテスト

仕様 (ソース中のコメント):
    p[i1] <= p[n1] <= p1 <= p2 <= p[n2] <= p[i2]
すなわち [p1, p2] を覆う格子点インデックス区間 [n1, n2] を返す。
区間が格子の外にある場合は n1 > n2 (空区間) を返す。
"""

import numpy as np
import pytest

import sol.material as material

EPS = 1.0e-6


@pytest.fixture
def grid():
    # p = [0, 1, 2, 3, 4, 5, 6, 7]
    return np.arange(8).astype('f8')


def test_node_aligned(grid):
    n1, n2 = material._getSpan(grid, len(grid), 0, 6, 2.0, 5.0, EPS)
    assert (n1, n2) == (2, 5)


def test_interior_between_nodes(grid):
    n1, n2 = material._getSpan(grid, len(grid), 0, 6, 2.3, 5.7, EPS)
    assert (n1, n2) == (2, 6)
    # 不変量: 区間 [p[n1], p[n2]] が [p1, p2] を包含する
    assert grid[n1] - EPS <= 2.3
    assert 5.7 <= grid[n2] + EPS


def test_below_range_returns_empty(grid):
    # [p1, p2] が格子下端より小さい -> n1 > n2 (空区間)
    n1, n2 = material._getSpan(grid, len(grid), 0, 6, -5.0, -3.0, EPS)
    assert n1 > n2


def test_above_range_returns_empty(grid):
    # [p1, p2] が格子上端より大きい -> n1 > n2 (空区間)
    n1, n2 = material._getSpan(grid, len(grid), 0, 6, 10.0, 12.0, EPS)
    assert n1 > n2


def test_single_point_range(grid):
    # i1 == i2 のとき n1 == n2 を即座に返す
    n1, n2 = material._getSpan(grid, len(grid), 3, 3, 0.0, 100.0, EPS)
    assert (n1, n2) == (3, 3)


def test_reversed_p1_p2_is_normalized(grid):
    # p1 > p2 で渡しても内部で min/max により正規化される
    a = material._getSpan(grid, len(grid), 0, 6, 5.0, 2.0, EPS)
    b = material._getSpan(grid, len(grid), 0, 6, 2.0, 5.0, EPS)
    assert a == b


def test_indices_clamped_to_grid(grid):
    # i1, i2 が範囲外でも [0, n-1] にクランプされ、結果が格子内に収まる
    n1, n2 = material._getSpan(grid, len(grid), -10, 100, 2.0, 5.0, EPS)
    assert (n1, n2) == (2, 5)
