// MIT License
//
// Copyright (c) 2024 TriLite https://github.com/MeshLite/TriLite
// Copyright 2020 Tomas Akenine-Möller
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

  /**
   * @brief Removes self-intersections from a triangular mesh.
   * @param mesh The triangular mesh to be processed.
   */
  static void RemoveSelfIntersections(Trimesh& mesh);

 private:
  struct BVHNode {
    Eigen::AlignedBox3d bbox_;
    std::optional<std::pair<TL::F, std::array<Vector3d, 3>>> triangle_;
    BVHNode* left_;
    BVHNode* right_;

    BVHNode() : left_(nullptr), right_(nullptr) {}
    ~BVHNode() {
      delete left_;
      delete right_;
    }
  };

  static BVHNode* ConstructBVH(
      std::vector<std::pair<TL::F, std::array<Vector3d, 3>>>& triangles,
      size_t start, size_t end, int depth = 0) {
    BVHNode* node = new BVHNode();
    assert(end > start);
    if (end - start == 1) {
      node->triangle_ = triangles[start];
      for (const auto& vertex : triangles[start].second) {
        node->bbox_.extend(vertex);
      }
    } else {
      for (size_t i = start; i < end; ++i) {
        for (const auto& vertex : triangles[i].second) {
          node->bbox_.extend(vertex);
        }
      }
      Vector3d extents = node->bbox_.sizes();
      int axis = 0;
      if (extents[1] > extents[0]) axis = 1;
      if (extents[2] > extents[axis]) axis = 2;
      std::sort(triangles.begin() + start, triangles.begin() + end,
                [axis](const std::pair<TL::F, std::array<Vector3d, 3>>& a,
                       const std::pair<TL::F, std::array<Vector3d, 3>>& b) {
                  return a.second[0][axis] < b.second[0][axis];
                });

      size_t mid = start + (end - start) / 2;
      node->left_ = ConstructBVH(triangles, start, mid, depth + 1);
      node->right_ = ConstructBVH(triangles, mid, end, depth + 1);
    }
    return node;
  }

  static bool DoesIntersect(
      BVHNode* node, const std::pair<TL::F, std::array<Vector3d, 3>>& tri) {
    if (!node || !node->bbox_.intersects(Eigen::AlignedBox3d(tri.second[0])
                                             .extend(tri.second[1])
                                             .extend(tri.second[2]))) {
      return false;
    }
    if (node->triangle_) {
      const auto& node_tri = node->triangle_.value();
      if (tri.first != node_tri.first &&
          Intersection::TriIntersectTri2(tri.second, node_tri.second)) {
        return true;
      }
    }
    return DoesIntersect(node->left_, tri) || DoesIntersect(node->right_, tri);
  }

  // Main function to find all auto-intersections in a mesh
  static std::vector<TL::F> FindSelfIntersections(TL::Trimesh& mesh,
                                                  double shrink_factor = 1e-8) {
    if (mesh.NumFaces() == 0) {
      return {};
    }
    std::vector<std::pair<TL::F, std::array<Vector3d, 3>>> triangles;
    for (TL::F f : mesh.Faces()) {
      if (mesh.FArea(f) > 1e-10) {
        Vector3d centroid = mesh.FCentroid(f);
        std::array<Vector3d, 3> tri;
        for (auto [i, position] : std::views::enumerate(mesh.FPositions(f))) {
          tri[i] = position + (centroid - position) * shrink_factor;
        }
        triangles.emplace_back(f, tri);
      }
    }
    std::vector<TL::F> self_intersection_faces;
    BVHNode* bvh_root = ConstructBVH(triangles, 0, triangles.size());
    for (const auto& tri : triangles) {
      if (DoesIntersect(bvh_root, tri)) {
        self_intersection_faces.push_back(tri.first);
      }
    }
    delete bvh_root;
    return self_intersection_faces;
  }
  class Intersection {
   public:
    static bool TriIntersectTri2(std::array<Vector3d, 3> tri1,
                                 std::array<Vector3d, 3> tri2) {
      for (int i : {0, 1, 2}) {
        if (SegmentIntersectsTriangle({tri1[i], tri1[(i + 1) % 3]}, tri2)) {
          return true;
        }
        // if (SegmentIntersectsTriangle({tri2[i], tri2[(i + 1) % 3]}, tri1)) {
        //   return true;
        // }
      }
      return false;
    }

   private:
    // Function to check if a point is inside a triangle using barycentric
    // coordinates

    static bool SegmentIntersectsTriangle(const std::array<Vector3d, 2>& seg,
                                          const std::array<Vector3d, 3>& tri) {
      Vector3d N = (tri[1] - tri[0]).cross(tri[2] - tri[0]).normalized();
      Vector3d ray_vector = (seg[1] - seg[0]).normalized();
      Vector3d edge1 = tri[1] - tri[0];
      Vector3d edge2 = tri[2] - tri[0];
      Vector3d ray_cross_e2 = ray_vector.cross(edge2);
      double det = edge1.dot(ray_cross_e2);
      if (abs(det) < 1e-11) {
        // Edge is coplanar with the triangle
        if (std::abs((tri[0] - seg[0]).dot(N)) < 1e-11) {
          return CoplanarSegmentIntersectsTriangle(seg, tri);
        }
      } else {
        // Edge is not coplanar with the triangle
        return NoCoplanarSegmentIntersectsTriangle(seg, tri);
      }
      return false;
    }
    static bool NoCoplanarSegmentIntersectsTriangle(
        const std::array<Vector3d, 2>& seg,
        const std::array<Vector3d, 3>& tri) {
      Vector3d ray_origin = seg[0];
      Vector3d ray_vector = (seg[1] - seg[0]).normalized();
      Vector3d edge1 = tri[1] - tri[0];
      Vector3d edge2 = tri[2] - tri[0];
      Vector3d ray_cross_e2 = ray_vector.cross(edge2);
      double det = edge1.dot(ray_cross_e2);
      assert(abs(det) > 1e-12);
      double inv_det = 1.0 / det;
      Vector3d s = ray_origin - tri[0];
      double u = inv_det * s.dot(ray_cross_e2);
      if (u < 0.0 || u > 1.0) {
        return false;
      }
      Vector3d s_cross_e1 = s.cross(edge1);
      double v = inv_det * ray_vector.dot(s_cross_e1);
      if (v < 0.0 || u + v > 1.0) {
        return false;
      }
      double t = inv_det * edge2.dot(s_cross_e1);
      if (t >= 0.0 && t <= (seg[1] - seg[0]).norm()) {
        return true;
      } else {
        return false;
      }
    }
    static bool EdgeEdgeTest(const Vector3d& v0, const Vector3d& u0,
                             const Vector3d& u1, short i0, short i1, double ax,
                             double ay) {
      double bx = u0[i0] - u1[i0];
      double by = u0[i1] - u1[i1];
      double cx = v0[i0] - u0[i0];
      double cy = v0[i1] - u0[i1];
      double f = ay * bx - ax * by;
      double d = by * cx - bx * cy;
      if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
        double e = ax * cy - ay * cx;
        if (f > 0) {
          if (e >= 0 && e <= f) return true;
        } else {
          if (e <= 0 && e >= f) return true;
        }
      }
      return false;
    }

    static bool EdgeAgainstTriEdges(const Vector3d& v0, const Vector3d& v1,
                                    const Vector3d& u0, const Vector3d& u1,
                                    const Vector3d& u2, short i0, short i1) {
      double ax = v1[i0] - v0[i0];
      double ay = v1[i1] - v0[i1];
      if (EdgeEdgeTest(v0, u0, u1, i0, i1, ax, ay)) {
        return true;
      }
      if (EdgeEdgeTest(v0, u1, u2, i0, i1, ax, ay)) {
        return true;
      }
      if (EdgeEdgeTest(v0, u2, u0, i0, i1, ax, ay)) {
        return true;
      }
      return false;
    }

    static bool CoplanarSegmentIntersectsTriangle(
        const std::array<Vector3d, 2>& seg,
        const std::array<Vector3d, 3>& tri) {
      short i0, i1;
      Vector3d edge = seg[1] - seg[0];
      Vector3d n = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
      Vector3d a = n.cwiseAbs();
      if (a[0] > a[1]) {
        if (a[0] > a[2]) {
          i0 = 1;
          i1 = 2;
        } else {
          i0 = 0;
          i1 = 1;
        }
      } else {
        if (a[1] > a[2]) {
          i0 = 0;
          i1 = 2;
        } else {
          i0 = 0;
          i1 = 1;
        }
      }
      if (EdgeEdgeTest(seg[0], tri[0], tri[1], i0, i1, edge[i0], edge[i1])) {
        return true;
      }
      if (EdgeEdgeTest(seg[0], tri[1], tri[2], i0, i1, edge[i0], edge[i1])) {
        return true;
      }
      if (EdgeEdgeTest(seg[0], tri[2], tri[0], i0, i1, edge[i0], edge[i1])) {
        return true;
      }
      if (PointInTri(seg[0], tri[0], tri[1], tri[2], i0, i1)) {
        return true;
      }
      return false;
    }

    static bool PointInTri(const Vector3d& p, const Vector3d& u0,
                           const Vector3d& u1, const Vector3d& u2, short i0,
                           short i1) {
      double a, b, c, d0, d1, d2;
      a = u1[i1] - u0[i1];
      b = -(u1[i0] - u0[i0]);
      c = -a * u0[i0] - b * u0[i1];
      d0 = a * p[i0] + b * p[i1] + c;

      a = u2[i1] - u1[i1];
      b = -(u2[i0] - u1[i0]);
      c = -a * u1[i0] - b * u1[i1];
      d1 = a * p[i0] + b * p[i1] + c;

      a = u0[i1] - u2[i1];
      b = -(u0[i0] - u2[i0]);
      c = -a * u2[i0] - b * u2[i1];
      d2 = a * p[i0] + b * p[i1] + c;

      if (d0 * d1 > 0.0 && d0 * d2 > 0.0) {
        return true;
      }
      return false;
    }
  };
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
  std::vector<Vector3d> new_positions(mesh.NumVertices());

  for (int iter = 0; iter < iterations; ++iter) {
    // Laplacian smoothing step (positive weight)
    for (V v : mesh.Vertices()) {
      Vector3d laplacian(0, 0, 0);
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
      Vector3d laplacian(0, 0, 0);
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

void Processing::RemoveSelfIntersections(Trimesh& mesh) {
  std::vector<TL::F> rm_faces;
  for (TL::F f : FindSelfIntersections(mesh)) {
    rm_faces.push_back(f);
  }
  mesh.RemoveFaces(rm_faces);
}

}  // namespace TL

#endif  // MESH_PROCESSING_HPP