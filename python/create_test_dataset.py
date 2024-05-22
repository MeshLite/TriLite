# Copyright (c) 2024 TriLite
# This file is part of the TriLite project, a C++23 library for triangular mesh
# processing. Distributed under the MIT License. The full license text can be
# found at: https://github.com/MeshLite/TriLite/blob/main/LICENSE
# This notice must remain intact in all copies or substantial portions of the
# file.

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
        TL.WriteMeshFile(mesh, filename, random.choice([False, True]))


if __name__ == "__main__":
    main()
