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

import unittest
import trilite as TL
import os
import sys

if len(sys.argv) < 2:
    print("Usage: python test_registration.py /path/to/your/dataset")
    sys.exit(1)

dataset_dir = sys.argv[1]


class TestMeshRegistration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dataset_dir = dataset_dir
        cls.files = os.listdir(cls.dataset_dir)

    def test_find_closest_points(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                indices, distances = TL.Registration.FindClosestPoints(
                    mesh, mesh
                )

    def test_find_best_rotation(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                rotation_matrix = TL.Registration.FindBestRotation(mesh, mesh)

    def test_find_best_rigid_transformation(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                rotation_matrix, translation_vector = (
                    TL.Registration.FindBestRigidTransformation(mesh, mesh)
                )

    def test_find_best_similarity_transformation(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                rotation_matrix, translation_vector, scaling_factor = (
                    TL.Registration.FindBestSimilarityTransformation(
                        mesh, mesh
                    )
                )

    def test_icp(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                rotation_matrix, translation_vector = TL.Registration.ICP(
                    mesh, mesh
                )

    def test_find_best_rigid_registration(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)
                mesh = TL.IO.ReadMeshFile(filepath)
                while mesh.NumFaces() > 1000:
                    TL.Processing.Simplify(mesh, 0.01, False, False)
                rotation_matrix, translation_vector = (
                    TL.Registration.RigidRegistrationHeuristics(mesh, mesh)
                )


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
