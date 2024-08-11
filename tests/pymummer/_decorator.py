# -*- coding: utf-8 -*-
"""
 * @Date: 2023-10-22 20:59:23
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 10:29:44
 * @FilePath: /genome/tests/genome/_decorator.py
 * @Description:
"""
# """

import sys
import pytest


MARK_LIMIT_RESOURCE = "-m not limit_resource" not in " ".join(sys.argv)
pytest_mark_resource = pytest.mark.skipif(
    MARK_LIMIT_RESOURCE, reason="no mark '-m \"not limit_resource\"'"
)
