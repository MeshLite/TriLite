// Copyright (c) 2024 TriLite
// This file is part of the TriLite project, a C++23 library for triangular mesh
// processing. Distributed under the MIT License. The full license text can be
// found at: https://github.com/MeshLite/TriLite/blob/main/LICENSE
// This notice must remain intact in all copies or substantial portions of the
// file.

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../TriLite/Core/Trimesh.hpp"
#include "../../TriLite/Modules/IO.hpp"
#include "../../TriLite/Modules/MeshProcessing.hpp"

namespace py = pybind11;
using namespace TL;

PYBIND11_MODULE(trilite, m) {
  py::class_<Trimesh>(m, "Trimesh")
      .def(py::init<>())
      .def(py::init<const std::vector<std::variant<V, Eigen::Vector3d>>&>())
      .def(py::init<const Trimesh&>())
      .def("NumHalfedges", &Trimesh::NumHalfedges)
      .def("NumVertices", &Trimesh::NumVertices)
      .def("NumFaces", &Trimesh::NumFaces)
      .def("Halfedges",
           [](const Trimesh& trimesh) {
             return std::ranges::to<std::vector>(trimesh.Halfedges());
           })
      .def("Vertices",
           [](const Trimesh& trimesh) {
             return std::ranges::to<std::vector>(trimesh.Vertices());
           })
      .def("Faces",
           [](const Trimesh& trimesh) {
             return std::ranges::to<std::vector>(trimesh.Faces());
           })
      .def("Positions",
           [](const Trimesh& trimesh) {
             return std::ranges::to<std::vector>(trimesh.Positions());
           })
      .def("HNext", &Trimesh::HNext)
      .def("HPrev", &Trimesh::HPrev)
      .def("HStart", &Trimesh::HStart)
      .def("HEnd", &Trimesh::HEnd)
      .def("HFace", &Trimesh::HFace)
      .def("HOpposite", &Trimesh::HOpposite)
      .def("HNextAroundStart", &Trimesh::HNextAroundStart)
      .def("HPrevAroundStart", &Trimesh::HPrevAroundStart)
      .def("HNextAroundEnd", &Trimesh::HNextAroundEnd)
      .def("HPrevAroundEnd", &Trimesh::HPrevAroundEnd)
      .def("HGeometry", &Trimesh::HGeometry)
      .def("HConnectionsAroundStart",
           [](const Trimesh& trimesh, H h) {
             return std::ranges::to<std::vector>(
                 trimesh.HConnectionsAroundStart(h));
           })
      .def("HHalfedgesAroundHole",
           [](const Trimesh& trimesh, H h) {
             return std::ranges::to<std::vector>(
                 trimesh.HHalfedgesAroundHole(h));
           })
      .def("HLength", &Trimesh::HLength)
      .def("VStarting", &Trimesh::VStarting)
      .def("VEnding", &Trimesh::VEnding)
      .def("VStartings",
           [](const Trimesh& trimesh, V v) {
             return std::ranges::to<std::vector>(trimesh.VStartings(v));
           })
      .def("VEndings",
           [](const Trimesh& trimesh, V v) {
             return std::ranges::to<std::vector>(trimesh.VEndings(v));
           })
      .def("VFaces",
           [](const Trimesh& trimesh, V v) {
             return std::ranges::to<std::vector>(trimesh.VFaces(v));
           })
      .def("VPosition", py::overload_cast<V>(&Trimesh::VPosition, py::const_))
      .def("VNormal", &Trimesh::VNormal)
      .def("VValence", &Trimesh::VValence)
      .def("VIsManifold", &Trimesh::VIsManifold)
      .def("FHalfedge", &Trimesh::FHalfedge)
      .def("FHalfedges",
           [](const Trimesh& trimesh, F f) {
             return std::ranges::to<std::vector>(trimesh.FHalfedges(f));
           })
      .def("FVertices",
           [](const Trimesh& trimesh, F f) {
             return std::ranges::to<std::vector>(trimesh.FVertices(f));
           })
      .def("FPositions",
           [](const Trimesh& trimesh, F f) {
             return std::ranges::to<std::vector>(trimesh.FPositions(f));
           })
      .def("FNormal", &Trimesh::FNormal)
      .def("FArea", &Trimesh::FArea)
      .def("EdgeHalfedges",
           [](const Trimesh& trimesh, H h) {
             return std::ranges::to<std::vector>(trimesh.EdgeHalfedges(h));
           })
      .def("EdgeFaces",
           [](const Trimesh& trimesh, H h) {
             return std::ranges::to<std::vector>(trimesh.EdgeFaces(h));
           })
      .def("EdgeIsManifold", &Trimesh::EdgeIsManifold)
      .def("BoundaryHalfedges",
           [](const Trimesh& trimesh) {
             return std::ranges::to<std::vector>(trimesh.BoundaryHalfedges());
           })
      .def("AddFace", &Trimesh::AddFace)
      .def("RemoveFace", &Trimesh::RemoveFace)
      .def("CollapseEdge", &Trimesh::CollapseEdge)
      .def("DisconnectFace", &Trimesh::DisconnectFace)
      .def("DisconnectFacesUntilManifoldEdges",
           &Trimesh::DisconnectFacesUntilManifoldEdges)
      .def("DisconnectFacesUntilManifoldVertices",
           &Trimesh::DisconnectFacesUntilManifoldVertices)
      .def("DisconnectFacesUntilManifold",
           &Trimesh::DisconnectFacesUntilManifold)
      .def("__copy__", [](const Trimesh& self) { return Trimesh(self); })
      .def("__deepcopy__",
           [](const Trimesh& self, py::dict) { return Trimesh(self); });

  m.def("ReadMeshFile", &TL::ReadMeshFile,
        "Read a mesh file and return a Trimesh object", py::arg("filepath"));
  m.def("WriteMeshFile", &TL::WriteMeshFile,
        "Write a Trimesh object to a mesh file", py::arg("mesh"),
        py::arg("filepath"), py::arg("binary_mode") = true);
  m.def("DecimateMesh", &TL::DecimateMesh,
        "A function to simplify a triangular mesh", py::arg("mesh"),
        py::arg("target_face_count"));
  m.def("FillMeshHoles", &TL::FillMeshHoles,
        "A function to remove holes from triangular mesh", py::arg("mesh"),
        py::arg("target_hole_count") = 0);
  m.def("TaubinSmoothing", &TL::TaubinSmoothing,
        "Apply Laplacian smoothing to a mesh", py::arg("mesh"),
        py::arg("iterations") = 1, py::arg("lambda") = 0.5,
        py::arg("mu") = -0.53);
}