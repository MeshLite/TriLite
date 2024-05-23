# MIT License
#
# Copyright (c) 2024 TriLite https:#github.com/MeshLite/TriLite
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

import unittest
import trilite as TL
import os
import sys
import numpy as np
from math import sqrt, log

if len(sys.argv) < 2:
    print("Usage: python test_trimesh_functions.py /path/to/your/dataset")
    sys.exit(1)

dataset_dir = sys.argv[1]


class TestMeshProcessing(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dataset_dir = dataset_dir
        cls.files = os.listdir(cls.dataset_dir)

    def test_mean_distance_to_barycenter(self):
        v0 = np.array([0.0, 0.0, 0.0])
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.5, 0.5 * sqrt(3), 0.0])
        mesh1 = TL.Trimesh([v0, v1, v2])
        bary = (v0 + v1 + v2) / 3.0
        mesh2 = TL.Trimesh(
            [bary, bary + [0.0, 0.001, 1.0], bary + [0.0, -0.001, 1.0]]
        )
        self.assertAlmostEqual(
            TL.Distance.AsymetricMeanEuclidean(mesh1, mesh2, 1e-2),
            (2 * sqrt(3) + log(sqrt(3) + 2)) / 18.0,
            delta=1e-3,
        )
        self.assertAlmostEqual(
            TL.Distance.AsymetricMeanEuclidean(mesh2, mesh1, 1e-2),
            2.0 / 3.0,
            delta=1e-3,
        )
        tree = TL.Distance.Tree(mesh1)
        self.assertAlmostEqual(tree.Distance([0.0, 0.0, 1.0]), 1.0)
        self.assertTrue(
            np.array_equal(
                tree.ClosestPoint([0.0, 0.0, 1.0]),
                np.array([0.0, 0.0, 0.0]),
            )
        )

    def test_mean_distance_to_mesh(self):
        ref_mesh = None
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                if ref_mesh == None:
                    ref_mesh = mesh
                TL.Distance.MeanEuclidean(mesh, ref_mesh, 1.0)


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
