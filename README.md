# TriLite: A C++23 Trimesh Library

TriLite is a modern C++23 library designed for efficient and easy manipulation and analysis of triangular meshes. It's built to be fast, flexible, and straightforward, making it suitable for both academic research and industrial applications in graphics, computational geometry, and 3D modeling. Each release is tested through a python binding on the Tringi10K dataset to ensure industrial robustness (see python/tests).

## Features

- **Modern C++23 Codebase**: Utilizes the latest language features for optimal performance and easier maintenance.
- **Header-Only Library**: Easy to integrate with other projects without the need for complicated build systems.
- **Efficient Memory Management**: Optimized for low memory overhead while handling large mesh datasets.
- **Extensive Mesh Formats**: Supports a variety of mesh formats to be easily incorporated in an existing project.
- **Robust API Documentation**: Comprehensive documentation making it easier for new users to get started.

## Getting Started

### Prerequisites

- C++23 compiler that supports std::generator (e.g., GCC 14+)
- If using python bindings (Makefile) : GCC 14+ and Eigen located in /usr/include/eigen3

### Installing

You can integrate TriLite into your project by including it as a submodule or directly adding the header files to your project.

```bash
git clone https://github.com/MeshLite/TriLite.git
```

You can test the python binding of this library by running make.
The virtual env where the trilite module is installed can be activated with :
```bash
source python/bindings/venv/bin/activate
```

### Usage

Here is a simple example on how to read a mesh, perform basic operation and write the result:

```cpp
#include "TriLite/Core/Trimesh.hpp"
#include "TriLite/Modules/IO.hpp"
#include "TriLite/Modules/MeshProcessing.hpp"
int main() {
  TL::Trimesh mesh = TL::ReadMeshFile("bunny.obj");
  TL::DecimateMesh(mesh, mesh.NumFaces() / 2);  // Simplify mesh by 50%
  TL::WriteMeshFile(mesh, "out.stl");
}
```

Or directly with python binding :
```python
    import trilite as TL
    mesh = TL.ReadMeshFile("bunny.obj)
    TL.DecimateMesh(mesh, mesh.NumFaces() // 2)
    TL.WriteMeshFile(mesh, "out.stl")
```

## Code Formatting

The C++ code in this library follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) recommandations.
The python code in this library follows the [PEP8](https://peps.python.org/pep-0008/) recommandations.

## Documentation

For more detailed information and API descriptions, visit the [documentation page](https://MeshLite.github.io/TriLite).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback, please reach out to us at [meshlite.developers@gmail.com](mailto:meshlite.developers@gmail.com).

## Support

If you find this library useful and would like to support its development, consider starring it on GitHub.