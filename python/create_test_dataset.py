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

import os
import numpy as np
import trilite as TL
import random
import sys

if len(sys.argv) < 2:
    print("Usage: python create_test_dataset.py /path/to/output/dataset")
    sys.exit(1)


def create_random_mesh(num_vertices, num_faces):
    mesh = TL.Trimesh()
    vertices = []
    for _ in range(num_vertices):
        vertices.append(np.array([random.random() for _ in range(3)]))
    for _ in range(num_faces):
        indices = random.sample(range(num_vertices), 3)
        mesh.AddFace([vertices[indices[i]] for i in range(3)])
    return mesh


def main():
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    num_meshes = 100
    input_dir = sys.argv[1]
    for i in range(num_meshes):
        num_vertices = random.randint(3, 100)
        num_faces = random.randint(0, 1000)
        mesh = create_random_mesh(num_vertices, num_faces)
        filename = os.path.join(input_dir, f"mesh_{i}.stl")
        TL.IO.WriteMeshFile(mesh, filename, random.choice([False, True]))
    mesh = TL.Trimesh()
    filename = os.path.join(input_dir, f"empty_mesh.stl")
    TL.IO.WriteMeshFile(mesh, filename)


if __name__ == "__main__":
    main()
