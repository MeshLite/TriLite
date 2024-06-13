// MIT License
//
// Copyright (c) 2024 TriLite https://github.com/MeshLite/TriLite
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../TriLite/Core/Trimesh.hpp"
#include "../../TriLite/Modules/Distance.hpp"
#include "../../TriLite/Modules/IO.hpp"
#include "../../TriLite/Modules/Processing.hpp"

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
             return ToVector<H>(trimesh.Halfedges());
           })
      .def("Vertices",
           [](const Trimesh& trimesh) {
             return ToVector<V>(trimesh.Vertices());
           })
      .def("Faces",
           [](const Trimesh& trimesh) { return ToVector<F>(trimesh.Faces()); })
      .def("Positions",
           [](const Trimesh& trimesh) {
             return ToVector<Vector3d>(trimesh.Positions());
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
      .def("HCentroid", &Trimesh::HCentroid)
      .def("HConnectionsAroundStart",
           [](const Trimesh& trimesh, H h) {
             return ToVector<H>(trimesh.HConnectionsAroundStart(h));
           })
      .def("HHalfedgesAroundHole",
           [](const Trimesh& trimesh, H h) {
             return ToVector<H>(trimesh.HHalfedgesAroundHole(h));
           })
      .def("HLength", &Trimesh::HLength)
      .def("VStarting", &Trimesh::VStarting)
      .def("VEnding", &Trimesh::VEnding)
      .def("VStartings", [](const Trimesh& trimesh,
                            V v) { return ToVector<H>(trimesh.VStartings(v)); })
      .def("VEndings", [](const Trimesh& trimesh,
                          V v) { return ToVector<H>(trimesh.VEndings(v)); })
      .def("VFaces", [](const Trimesh& trimesh,
                        V v) { return ToVector<F>(trimesh.VFaces(v)); })
      .def("VPosition", py::overload_cast<V>(&Trimesh::VPosition, py::const_))
      .def("VNormal", &Trimesh::VNormal)
      .def("VValence", &Trimesh::VValence)
      .def("VIsManifold", &Trimesh::VIsManifold)
      .def("VIsBoundary", &Trimesh::VIsBoundary)
      .def("FHalfedge", &Trimesh::FHalfedge)
      .def("FHalfedges", [](const Trimesh& trimesh,
                            F f) { return ToVector<H>(trimesh.FHalfedges(f)); })
      .def("FNeighbors", [](const Trimesh& trimesh,
                            F f) { return ToVector<H>(trimesh.FNeighbors(f)); })
      .def("FVertices", [](const Trimesh& trimesh,
                           F f) { return ToVector<V>(trimesh.FVertices(f)); })
      .def("FPositions",
           [](const Trimesh& trimesh, F f) {
             return ToVector<Vector3d>(trimesh.FPositions(f));
           })
      .def("FNormal", &Trimesh::FNormal)
      .def("FBoundingBox", &Trimesh::FBoundingBox)
      .def("FCentroid", &Trimesh::FCentroid)
      .def("FArea", &Trimesh::FArea)
      .def("EdgeHalfedges",
           [](const Trimesh& trimesh, H h) {
             return ToVector<H>(trimesh.EdgeHalfedges(h));
           })
      .def("EdgeFaces", [](const Trimesh& trimesh,
                           H h) { return ToVector<F>(trimesh.EdgeFaces(h)); })
      .def("EdgeIsManifold", &Trimesh::EdgeIsManifold)
      .def("BoundaryHalfedges",
           [](const Trimesh& trimesh) {
             return ToVector<H>(trimesh.BoundaryHalfedges());
           })
      .def("MedianEdgeLength", &Trimesh::MedianEdgeLength)
      .def("BoundingBox", &Trimesh::BoundingBox)
      .def("Centroid", &Trimesh::Centroid)
      .def("AddFace", &Trimesh::AddFace)
      .def("RemoveFace", &Trimesh::RemoveFace)
      .def("CollapseEdge", &Trimesh::CollapseEdge)
      .def("FlipHalfedgeWithOpposite", &Trimesh::FlipHalfedgeWithOpposite)
      .def("SplitEdge", &Trimesh::SplitEdge)
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

  py::class_<IO>(m, "IO")
      .def_static("ReadMeshFile", &IO::ReadMeshFile,
                  "Read a mesh file and return a Trimesh object",
                  py::arg("filepath"))
      .def_static("WriteMeshFile", &IO::WriteMeshFile,
                  "Write a Trimesh object to a mesh file", py::arg("mesh"),
                  py::arg("filepath"), py::arg("binary_mode") = true);

  py::class_<Processing>(m, "Processing")
      .def_static("Decimate", &Processing::Decimate,
                  "External decimation function of a triangular mesh by edge "
                  "length ordering.",
                  py::arg(" mesh "), py::arg("target_face_count"))
      .def_static(
          "Simplify", &Processing::Simplify,
          "Simplifies a given triangular mesh by collapsing edges based on a "
          "quadric error metric. Only valid collapses are performed.",
          py::arg(" mesh "), py::arg("max_collapse_cost_ratio") = 0.05,
          py::arg("preserve_boundaries") = false)
      .def_static("FillHoles", &Processing::FillHoles,
                  "A function to remove holes from triangular mesh",
                  py::arg("mesh"), py::arg("target_hole_count") = 0)
      .def_static("TaubinSmoothing", &Processing::TaubinSmoothing,
                  "Apply Laplacian smoothing to a mesh", py::arg("mesh"),
                  py::arg("iterations") = 1, py::arg("lambda") = 0.5,
                  py::arg("mu") = -0.53)
      .def_static(
          "RemoveSelfIntersections", &Processing::RemoveSelfIntersections,
          "A function to remove self-intersections in a triangular mesh",
          py::arg("mesh"))
      .def_static("PrintabilityHeuristics", &Processing::PrintabilityHeuristics,
                  "Prepares a triangular mesh for 3D printing by iteratively "
                  "closing it while applying several cleaning and repair steps "
                  "(only the largest connected component is preserved)",
                  py::arg("mesh"), py::arg("niters") = 10);

  py::class_<Distance>(m, "Distance")
      .def_static(
          "AsymmetricHausdorff", &Distance::AsymmetricHausdorff,
          "Compute the asymmetric Hausdorff distance between two meshes",
          py::arg("mesh"), py::arg("target_mesh"), py::arg("precision"))
      .def_static("Hausdorff", &Distance::Hausdorff,
                  "Compute the Hausdorff distance between two meshes",
                  py::arg("mesh1"), py::arg("mesh2"), py::arg("precision"))
      .def_static(
          "AsymmetricMeanEuclidean", &Distance::AsymmetricMeanEuclidean,
          "Compute the asymmetric mean Euclidean distance between two meshes",
          py::arg("mesh"), py::arg("target_mesh"), py::arg("precision"))
      .def_static("MeanEuclidean", &Distance::MeanEuclidean,
                  "Compute the mean Euclidean distance between two meshes",
                  py::arg("mesh1"), py::arg("mesh2"), py::arg("precision"));

  py::class_<Distance::Tree>(m.attr("Distance"), "Tree")
      .def(py::init<const Trimesh&>())
      .def("Distance", &Distance::Tree::Distance,
           "Compute the unsigned Euclidean distance from a point to the "
           "nearest point on the mesh",
           py::arg("point"))
      .def("ClosestPoint", &Distance::Tree::ClosestPoint,
           "Find the closest point on the mesh to a given point",
           py::arg("point"));
}