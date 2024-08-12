# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-08-12 13:12:52
 * @FilePath: /pymummer/tests/setup.py
 * @Description:
"""

import os
from pathlib import Path

repo_path = Path(__file__).parent.parent

os.chdir(repo_path)


def get_version(file: str | Path):
    "first try get version via versionneer, then try to get version from changelog.md"
    try:
        from tests._version import get_versions

        version = get_versions()["version"]
        if version:
            return version
    except ImportError:
        pass
    with open(file) as f:
        for line in f:
            if line.startswith("## changelog"):
                break
        for line in f:
            if line.strip():
                v = line.strip().rsplit(maxsplit=1)
                break
        else:
            v = ["", "0+unknown"]
        version = v[1].rstrip(":")
    return version


if __name__ == "__main__":
    from setuptools import setup, find_packages

    package_name = "pymummer"
    package_description = "handle mummer delta file"

    setup(
        name=package_name,
        version=get_version(repo_path / "README.md"),
        author="hwrn.aou",
        author_email="hwrn.aou@sjtu.edu.cn",
        description=package_description,
        packages=find_packages(include=[f"{package_name}*"]),
        include_package_data=True,
    )
