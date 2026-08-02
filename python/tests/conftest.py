# -*- coding: utf-8 -*-
"""
pytest 共通設定

OpenFDTD の Python ソルバーは ``import sol.xxx`` / ``import sol_cuda.xxx`` の
形式でモジュールを参照する。テストをどのディレクトリから起動しても
これらが解決できるよう、``python/`` ディレクトリを sys.path 先頭に追加する。
"""

import os
import sys

PYTHON_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)
