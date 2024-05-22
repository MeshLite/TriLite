# Copyright (c) 2024 TriLite
# This file is part of the TriLite project, a C++23 library for triangular mesh
# processing. Distributed under the MIT License. The full license text can be
# found at: https://github.com/MeshLite/TriLite/blob/main/LICENSE
# This notice must remain intact in all copies or substantial portions of the
# file.

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
