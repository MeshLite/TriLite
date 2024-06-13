![TriLite Mesh Repairing](https://raw.githubusercontent.com/MeshLite/TriLite/images/scan_mesh_repair.png)

*Figure : The TriLite library provides tools to repair a [scanned mesh](http://redwood-data.org/3dscan/models.html?i=1833), resulting in a watertight mesh.*

TriLite is a header-only C++ library designed for efficient manipulation of triangle meshes. It offers a variety of functionalities including mesh reading/writing, distance calculations, and mesh processing algorithms, all without relying on a garbage collector. This ensures that the mesh topology remains consistent across various operations, regardless of the format used for writing and reading back the mesh.

Features
========

- **Mesh I/O**: Supports multiple mesh file formats including STL, OBJ, OFF, and PLY.
- **Distance Calculations**: Includes functions to compute Hausdorff and mean Euclidean distances between meshes.
- **Mesh Processing**: Provides tools for mesh decimation, hole filling, smoothing, and repairing for 3D printing.
- **Consistency**: Ensures the same mesh topology is maintained when writing and reading back meshes.

No Garbage Collector
====================

A unique aspect of TriLite is its design choice to avoid using a garbage collector. This means:
- **Compact Mesh Structure**: There can't be at any time any isolated vertex or any face marked as deleted (as with a garbage collector). If a face is deleted or a vertex becoming isolated, it is swapped with the last element and popped back to remove it definitively from the structure in O(L) (with L the largest vertex valence). Any order of operations always leads to a compact mesh structure.
- **Consistent Mesh Topology**: When writing a mesh to a file and reading it back, the topology remains unchanged. For example, when writing an STL file, points that coincide are slightly adjusted to prevent merging when read back, preserving the original mesh structure.

Installation
============

### C++ Usage
To integrate TriLite into your C++ project, include the header files and compile with the Eigen flag.

### Python Usage
For Python, you can install the latest released version via pip:

```sh
sudo apt-get install -y g++-12 libeigen3-dev
pip install trilite
```

Alternatively, you can build, install, and test the module using the provided Makefile:

```sh
make
```

Prerequisites
=============

- g++-12 or later
- Eigen3
- Python 3.x
- pybind11
- numpy

Usage
=====

### C++ Example

```cpp
#include "TriLite/Modules/Distance.hpp"
#include "TriLite/Modules/IO.hpp"
#include "TriLite/Modules/Processing.hpp"

int main() {
  // Reading a mesh
  TL::Trimesh mesh1 = TL::IO::ReadMeshFile("input.obj");

  // Copying the mesh
  TL::Trimesh mesh2(mesh1);

  // Simplify the copied mesh
  TL::Processing::Simplify(mesh2);

  // Writing the simplified mesh back to a file
  TL::IO::WriteMeshFile(mesh2, "output.stl");

  // Calculating Hausdorff distance between the meshes
  double distance = TL::Distance::Hausdorff(mesh1, mesh2, 1e-6);
  std::cout << "Hausdorff distance: " << distance << std::endl;

  return 0;
}
```

### Python Example

```python
import trilite as TL

# Reading a mesh
mesh1 = TL.IO.ReadMeshFile("input.obj")

# Copying the mesh
mesh2 = TL.Trimesh(mesh1)

# Simplify the copied mesh
TL.Processing.Simplify(mesh2)

# Writing the simplified mesh back to a file
TL.IO.WriteMeshFile(mesh2, "output.stl")

# Calculating Hausdorff distance between the meshes
distance = TL.Distance.Hausdorff(mesh1, mesh2, 1e-6)
print(f"Hausdorff distance: {distance}")
```

Testing
=======

TriLite is rigorously tested over the entire [Thingi10k dataset](https://ten-thousand-models.appspot.com/). The Python tests for these are located in the `python/tests` directory. You can run these tests using:

```sh
make DATASET=/path/to/Thingi10k/raw_meshes
```

Coding Standards
================

- The C++ code follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
- The Python code follows [PEP8](https://peps.python.org/pep-0008/).

License
=======

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Contributing
============

Contributions are welcome! Please fork the repository and submit a pull request.

Acknowledgements
================

TriLite is developed and maintained by the MeshLite team. Special thanks to all the contributors and users who have supported this project.

Contact
=======

For any questions or suggestions, please contact us at meshlite.developers@gmail.com.

Support
=======

If you find this project useful, please consider giving it a star on GitHub. Your support is appreciated!