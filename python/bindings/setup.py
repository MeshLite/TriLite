# MIT License
#
# Copyright (c) 2024 TriLite https://github.com/MeshLite/TriLite
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
