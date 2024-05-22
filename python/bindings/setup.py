# Copyright (c) 2024 TriLite
# This file is part of the TriLite project, a C++23 library for triangular mesh
# processing. Distributed under the MIT License. The full license text can be
# found at: https://github.com/MeshLite/TriLite/blob/main/LICENSE
# This notice must remain intact in all copies or substantial portions of the
# file.

from setuptools import setup, Extension
import pybind11
import os

# Specify the compiler and linker
os.environ["CC"] = "g++-14"
os.environ["CXX"] = "g++-14"

ext_modules = [
    Extension(
        "trilite",
        ["trilite_module.cpp"],
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(user=True),
            "/usr/include/eigen3",
        ],
        language="c++",
        extra_compile_args=["-std=c++23"],
    ),
]

setup(
    name="trilite",
    version="0.1",
    author="MeshLite",
    author_email="meshlite.developers@gmail.com",
    description="Python binding for the TriLite library",
    ext_modules=ext_modules,
    install_requires=["pybind11", "numpy"],
)
