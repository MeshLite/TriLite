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

if len(sys.argv) < 3:
    print("Usage: python test_trimesh_functions.py /input/dataset /output/dir")
    sys.exit(1)

dataset_dir = sys.argv[1]
output_dir = sys.argv[2]


class TestTrimeshFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dataset_dir = dataset_dir
        cls.output_dir = output_dir
        cls.files = os.listdir(cls.dataset_dir)
        os.makedirs(cls.output_dir, exist_ok=True)

    def test_all_models(self):
        for filename in os.listdir(self.__class__.dataset_dir):
            with self.subTest(filename=filename):
                filepath = os.path.join(self.__class__.dataset_dir, filename)

                mesh = TL.ReadMeshFile(filepath)

                structure = [[] for _ in range(mesh.NumVertices())]
                for v in range(mesh.NumVertices()):
                    structure[v] = mesh.VStartings(v)

                for extension in [".stl", ".obj", ".off", ".ply"]:
                    for binary_mode in [True, False]:
                        self.__class__.binary_mode = binary_mode
                        self.__class__.cur_path = os.path.join(
                            self.__class__.output_dir, f"test{extension}"
                        )
                        TL.WriteMeshFile(
                            mesh, self.__class__.cur_path, binary_mode
                        )
                        rounded_mesh = TL.ReadMeshFile(self.__class__.cur_path)
                        # The mesh topology must remain exactly the same
                        # as the initial mesh, for all file extensions
                        self.assertEqual(
                            rounded_mesh.NumVertices(), len(structure)
                        )
                        for v in range(rounded_mesh.NumVertices()):
                            self.assertEqual(
                                rounded_mesh.VStartings(v), structure[v]
                            )


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
