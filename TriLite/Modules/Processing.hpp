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

#ifndef MESH_PROCESSING_HPP
#define MESH_PROCESSING_HPP

#include <set>

#include "../Core/Trimesh.hpp"
namespace TL {

class Processing {
 public:
  /**
   * @brief External decimation function to simplify a triangular mesh.
   * @param mesh Reference to the Trimesh object to be simplified.
   * @param target_face_count The target number of faces after decimation.
   */
  static void DecimateMesh(Trimesh& mesh, size_t target_face_count);

  /**
   * @brief Fills holes in a triangular mesh by detecting boundary edges and
   * adding triangles.
   * @param mesh The triangular mesh to be processed.
   */
  static void FillMeshHoles(Trimesh& mesh, size_t target_hole_count = 0);

  /**
   * @brief Applies volume-preserving Laplacian smoothing (Taubin smoothing) to
   * a triangular mesh.
   * @param mesh The triangular mesh to be smoothed.
   * @param iterations The number of smoothing iterations.
   * @param lambda The smoothing factor for the Laplacian step.
   * @param mu The inverse smoothing factor for the inverse Laplacian step.
   */
  static void TaubinSmoothing(Trimesh& mesh, int iterations = 1,
                              double lambda = 0.5, double mu = -0.53);
};

void Processing::DecimateMesh(Trimesh& mesh, size_t target_face_count) {
  std::vector<double> edge_lengths(mesh.NumHalfedges());
  auto compare = [&edge_lengths](H h, H g) {
    return edge_lengths[h] < edge_lengths[g] ||
           (edge_lengths[h] == edge_lengths[g] && h < g);
  };
  std::set<H, decltype(compare)> hset(compare);
  for (H h : mesh.Halfedges()) {
    edge_lengths[h] = mesh.HLength(h);
    assert(std::isfinite(edge_lengths[h]));
    hset.insert(h);
  }
  while (mesh.NumFaces() > target_face_count && !hset.empty()) {
    auto minh = *hset.begin();
    hset.erase(minh);
    if (minh >= mesh.NumHalfedges()) {
      continue;
    }
    std::array<V, 2> verts = {mesh.HStart(minh), mesh.HEnd(minh)};
    V last_vert_id = mesh.NumVertices() - 1;
    auto [rm_faces, rm_vertices] = mesh.CollapseEdge(minh);
    for (const V& rem_v : rm_vertices) {
      for (V& v : verts) {
        if (v == rem_v) {
          v = kInvalidId;
        } else if (v == last_vert_id) {
          v = rem_v;
        }
      }
      last_vert_id--;
    }
    for (F f : rm_faces) {
      for (H h : {3 * f, 3 * f + 1, 3 * f + 2}) {
        if (h < mesh.NumHalfedges()) {
          hset.erase(h);
        }
      }
    }
    for (V v : verts) {
      if (v != kInvalidId) {
        for (H h : mesh.VStartings(v)) {
          hset.erase(h);
        }
      }
    }
    for (F f : rm_faces) {
      for (H h : {3 * f, 3 * f + 1, 3 * f + 2}) {
        if (h < mesh.NumHalfedges()) {
          edge_lengths[h] = mesh.HLength(h);
          hset.insert(h);
        }
      }
    }
    for (V v : verts) {
      if (v != kInvalidId) {
        for (H h : mesh.VStartings(v)) {
          edge_lengths[h] = mesh.HLength(h);
          hset.insert(h);
        }
      }
    }
  }
}

void Processing::FillMeshHoles(Trimesh& mesh, size_t target_hole_count) {
  std::unordered_set<H> boundary_edges;
  std::vector<std::vector<H>> polygons;
  for (H st_h : mesh.BoundaryHalfedges()) {
    if (!boundary_edges.count(st_h)) {
      polygons.push_back(
          std::ranges::to<std::vector>(mesh.HHalfedgesAroundHole(st_h)));
      boundary_edges.insert(polygons.back().begin(), polygons.back().end());
    }
  }
  std::sort(polygons.begin(), polygons.end(),
            [](auto& A, auto& B) { return A.size() < B.size(); });
  for (const auto& [id, polygon] : std::views::enumerate(polygons)) {
    if (polygons.size() - id == target_hole_count) {
      break;
    }
    assert(polygon.size() >= 3);
    V v0 = mesh.HStart(polygon[0]);
    for (size_t i = 1; i < polygon.size() - 1; ++i) {
      V v1 = mesh.HStart(polygon[i]);
      V v2 = mesh.HStart(polygon[i + 1]);
      mesh.AddFace({v0, v2, v1});
    }
  }
}

void Processing::TaubinSmoothing(Trimesh& mesh, int iterations, double lambda,
                                 double mu) {
  std::vector<Eigen::Vector3d> new_positions(mesh.NumVertices());

  for (int iter = 0; iter < iterations; ++iter) {
    // Laplacian smoothing step (positive weight)
    for (V v : mesh.Vertices()) {
      Eigen::Vector3d laplacian(0, 0, 0);
      size_t valence = 0;
      for (H h : mesh.VStartings(v)) {
        laplacian += mesh.VPosition(mesh.HEnd(h));
        ++valence;
      }
      laplacian /= static_cast<double>(valence);
      new_positions[v] =
          mesh.VPosition(v) + lambda * (laplacian - mesh.VPosition(v));
    }
    for (V v : mesh.Vertices()) {
      mesh.VPosition(v) = new_positions[v];
    }

    // Laplacian smoothing step (negative weight, inverse smoothing)
    for (V v : mesh.Vertices()) {
      Eigen::Vector3d laplacian(0, 0, 0);
      size_t valence = 0;
      for (H h : mesh.VStartings(v)) {
        laplacian += mesh.VPosition(mesh.HEnd(h));
        ++valence;
      }
      laplacian /= static_cast<double>(valence);
      new_positions[v] =
          mesh.VPosition(v) + mu * (laplacian - mesh.VPosition(v));
    }
    for (V v : mesh.Vertices()) {
      mesh.VPosition(v) = new_positions[v];
    }
  }
}
}  // namespace TL

#endif  // MESH_PROCESSING_HPP