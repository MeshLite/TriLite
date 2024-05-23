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

if len(sys.argv) < 2:
    print("Usage: python test_trimesh_functions.py /path/to/your/dataset")
    sys.exit(1)

dataset_dir = sys.argv[1]


class TestMeshProcessing(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dataset_dir = dataset_dir
        cls.files = os.listdir(cls.dataset_dir)

    def test_decimation(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.ReadMeshFile(filepath)
                target_face_count = mesh.NumFaces() // 2
                TL.DecimateMesh(mesh, target_face_count)
                self.assertLessEqual(
                    mesh.NumFaces(),
                    target_face_count,
                    "Mesh faces should be decimated",
                )

    def test_hole_filling(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.ReadMeshFile(filepath)
                mesh.DisconnectFacesUntilManifold()
                for h in range(mesh.NumHalfedges()):
                    self.assertTrue(mesh.EdgeIsManifold(h))
                for v in range(mesh.NumVertices()):
                    self.assertTrue(mesh.VIsManifold(v))

                TL.FillMeshHoles(mesh, 0)

                for h in range(mesh.NumHalfedges()):
                    self.assertTrue(mesh.HOpposite(h) < mesh.NumHalfedges())

    def test_taubin_smoothing(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.ReadMeshFile(filepath)
                TL.TaubinSmoothing(mesh, 1)


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
