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

if len(sys.argv) < 3:
    print("Usage: python test_trimesh_functions.py /input/dataset /output/dir")
    sys.exit(1)

dataset_dir = sys.argv[1]
output_dir = sys.argv[2]


class TestIO(unittest.TestCase):

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

                mesh = TL.IO.ReadMeshFile(filepath)
                TL.Processing.Simplify(mesh, 0.05, True)
                structure = [[] for _ in range(mesh.NumVertices())]
                for v in range(mesh.NumVertices()):
                    structure[v] = mesh.VStartings(v)

                for extension in [".stl", ".obj", ".off", ".ply"]:
                    for binary_mode in [True, False]:
                        self.__class__.binary_mode = binary_mode
                        self.__class__.cur_path = os.path.join(
                            self.__class__.output_dir, f"test{extension}"
                        )
                        TL.IO.WriteMeshFile(
                            mesh, self.__class__.cur_path, binary_mode
                        )
                        rounded_mesh = TL.IO.ReadMeshFile(
                            self.__class__.cur_path
                        )
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
