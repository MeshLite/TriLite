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

  /**
   * @brief Prepares a triangular mesh for 3D printing by iteratively closing it
   * while applying several cleaning and repair steps (only the largest
   * connected component is preserved).
   * @param mesh Reference to the Trimesh object to be prepared.
   * @param niters The number of iterations to perform (default is 10).
   */
  static void PrintabilityHeuristics(Trimesh& mesh, int niters = 10);

 private:
  static void ClampEdgeLengths(Trimesh& mesh, double min_length,
                               double max_length, int niters = 5) {
    for (int i = 0; i < niters; i++) {
      bool stop = true;
      for (H h : mesh.Halfedges()) {
        if (mesh.HLength(h) > max_length || mesh.HLength(h) < min_length) {
          stop = false;
          break;
        }
      }
      if (stop) {
        break;
      }
      F nf = mesh.NumFaces();
      std::vector<V> edge_to_v(mesh.NumHalfedges(), kInvalidId);
      std::vector<F> rm_faces;
      std::vector<std::vector<H>> edge_halfedges(mesh.NumHalfedges());
      for (H h : mesh.Halfedges()) {
        edge_halfedges[h] = std::ranges::to<std::vector>(mesh.EdgeHalfedges(h));
      }
      for (F f = 0; f < nf; f++) {
        int ct = 0;
        H he{kInvalidId};
        for (H h : mesh.FHalfedges(f)) {
          if (edge_to_v[h] != kInvalidId ||
              mesh.HLength(edge_halfedges[h].front()) > max_length) {
            ++ct;
            he = h;
          }
        }
        if (ct > 1) {
          rm_faces.push_back(f);
          std::vector<V> verts = {mesh.HStart(3 * f), mesh.HStart(3 * f + 1),
                                  mesh.HStart(3 * f + 2)};
          std::array<std::variant<V, Vector3d>, 3> centroids = {
              mesh.HCentroid(3 * f), mesh.HCentroid(3 * f + 1),
              mesh.HCentroid(3 * f + 2)};
          V num_v = mesh.NumVertices();
          for (int j : {0, 1, 2}) {
            if (edge_to_v[3 * f + j] != kInvalidId) {
              assert(edge_to_v[3 * f + j] < mesh.NumVertices());
              centroids[j] = edge_to_v[3 * f + j];
            } else {
              bool flag = false;
              for (H g : edge_halfedges[3 * f + j]) {
                flag |= (g == 3 * f + j);
                edge_to_v[g] = num_v;
              }
              assert(flag);
              ++num_v;
            }
          }
          mesh.AddFace(centroids);
          mesh.AddFace({verts[0], edge_to_v[3 * f], edge_to_v[3 * f + 2]});
          mesh.AddFace({verts[1], edge_to_v[3 * f + 1], edge_to_v[3 * f]});
          mesh.AddFace({verts[2], edge_to_v[3 * f + 2], edge_to_v[3 * f + 1]});
        } else if (ct == 1) {
          rm_faces.push_back(f);
          std::variant<V, Vector3d> centroid{mesh.HCentroid(he)};
          if (edge_to_v[he] != kInvalidId) {
            centroid = edge_to_v[he];
          } else {
            for (H g : edge_halfedges[he]) {  // mesh.EdgeHalfedges(he)) {
              edge_to_v[g] = mesh.NumVertices();
            }
          }
          mesh.AddFace(
              {mesh.HStart(mesh.HPrev(he)), mesh.HStart(he), centroid});
          centroid = edge_to_v[he];
          mesh.AddFace({mesh.HEnd(he), mesh.HStart(mesh.HPrev(he)), centroid});
        }
      }
      mesh.RemoveFaces(rm_faces);
      CollapseSmallEdges(mesh, min_length);
    }
  }
  static void CollapseSmallEdges(Trimesh& mesh, double min_length) {
    for (F f = 0; f < mesh.NumFaces();) {
      bool flag = true;
      for (H h : mesh.FHalfedges(f)) {
        if (mesh.HLength(h) < min_length) {
          mesh.CollapseEdge(h);
          flag = false;
          break;
        }
      }
      f += flag;
    }
  }
  static std::pair<TL::Vector3d, double> CalculateCircumsphere(
      const TL::Vector3d& a, const TL::Vector3d& b, const TL::Vector3d& c) {
    TL::Vector3d ac = c - a;
    TL::Vector3d ab = b - a;
    TL::Vector3d abxac = ab.cross(ac);

    TL::Vector3d to_circumsphere_center = (abxac.cross(ab) * ac.squaredNorm() +
                                           ac.cross(abxac) * ab.squaredNorm()) /
                                          (2.0 * abxac.squaredNorm());
    double circumsphere_radius = to_circumsphere_center.norm();
    TL::Vector3d ccs = a + to_circumsphere_center;

    return std::make_pair(ccs, circumsphere_radius);
  }
  static bool IsDelaunay(const TL::Trimesh& mesh, TL::H h) {
    TL::H hopp = mesh.HOpposite(h);
    if (hopp == TL::kInvalidId) {
      return true;
    }
    TL::Vector3d a = mesh.VPosition(mesh.HStart(h));
    TL::Vector3d b = mesh.VPosition(mesh.HEnd(h));
    TL::Vector3d c = mesh.VPosition(mesh.HEnd(mesh.HNext(h)));
    TL::Vector3d d = mesh.VPosition(mesh.HEnd(mesh.HNext(hopp)));
    std::pair<TL::Vector3d, double> sphere = CalculateCircumsphere(a, b, c);
    return (sphere.first - d).norm() >= (sphere.second);
  }
  static void MakeDelaunay(TL::Trimesh& mesh) {
    std::queue<TL::H> edge_queue;
    for (TL::H h : mesh.Halfedges()) {
      if (mesh.HOpposite(h) != TL::kInvalidId &&
          mesh.HOpposite(h) / 3 != mesh.HOpposite(mesh.HNext(h)) / 3) {
        edge_queue.push(h);
      }
    }
    size_t count = 0;
    size_t limit = 3 * mesh.NumFaces();
    while (!edge_queue.empty()) {
      if (count++ == limit) {
        break;
      }
      TL::H h = edge_queue.front();
      edge_queue.pop();
      TL::H hopp = mesh.HOpposite(h);
      if (hopp != TL::kInvalidId && !IsDelaunay(mesh, h)) {
        TL::V v1 = mesh.HStart(mesh.HPrev(h));
        TL::V v2 = mesh.HStart(mesh.HPrev(hopp));
        bool valid = true;
        for (TL::H he : mesh.VStartings(v1)) {
          valid &= (mesh.HEnd(he) != v2);
        }
        for (TL::H he : mesh.VStartings(v2)) {
          valid &= (mesh.HEnd(he) != v1);
        }
        if (valid) {
          mesh.FlipHalfedgeWithOpposite(h);
          edge_queue.push(mesh.HNext(h));
          edge_queue.push(mesh.HPrev(h));
          edge_queue.push(mesh.HNext(hopp));
          edge_queue.push(mesh.HPrev(hopp));
        }
      }
    }
  }
  static void RetainLargestComponent(TL::Trimesh& mesh) {
    std::vector<bool> visited(mesh.NumFaces(), false);
    std::vector<std::vector<F>> components;
    for (F f = 0; f < mesh.NumFaces(); ++f) {
      if (!visited[f]) {
        std::vector<F> component;
        std::vector<F> stack;
        stack.push_back(f);
        while (!stack.empty()) {
          F current = stack.back();
          stack.pop_back();
          if (visited[current]) continue;
          visited[current] = true;
          component.push_back(current);
          for (F neighbor : mesh.FNeighbors(current)) {
            if (!visited[neighbor]) {
              stack.push_back(neighbor);
            }
          }
        }
        components.push_back(component);
      }
    }
    size_t max_id = 0;
    for (size_t i = 1; i < components.size(); i++) {
      if (components[i].size() > components[max_id].size()) {
        max_id = i;
      }
    }
    std::vector<F> faces_to_delete;
    for (size_t i = 0; i < components.size(); i++) {
      if (i != max_id) {
        for (F f : components[i]) {
          faces_to_delete.push_back(f);
        }
      }
    }
    mesh.RemoveFaces(faces_to_delete);
  }

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
    std::vector<double> sum_angle(polygon.size());
    auto set_sum_angle = [&mesh, &polygon, &sum_angle](size_t i) {
      H h = polygon[i];
      sum_angle[i] = 0;
      for (H he : mesh.HConnectionsAroundStart(h)) {
        sum_angle[i] += std::acos(
            std::clamp(mesh.HGeometry(he).normalized().dot(
                           -mesh.HGeometry(mesh.HPrev(he)).normalized()),
                       -1.0, 1.0));
      }
    };
    auto cmp = [&mesh, &sum_angle](std::pair<size_t, double> a,
                                   std::pair<size_t, double> b) {
      return a.second < b.second;
    };
    std::priority_queue<std::pair<size_t, double>,
                        std::vector<std::pair<size_t, double>>, decltype(cmp)>
        q(cmp);
    for (size_t i = 0; i < polygon.size(); i++) {
      set_sum_angle(i);
      q.push(std::make_pair(i, sum_angle[i]));
    }

    std::set<std::pair<size_t, H>> s;
    for (size_t i = 0; i < polygon.size(); i++) {
      s.insert(std::make_pair(i, polygon[i]));
    }
    while (!q.empty() && s.size() >= 3) {
      auto [i1, angle] = q.top();
      q.pop();
      if (angle != sum_angle[i1]) {
        continue;
      }
      H h1 = polygon[i1];
      auto it = s.find(std::make_pair(i1, h1));
      if (it == s.end()) {
        continue;
      }
      auto [i0, h0] = *((it == s.begin()) ? std::prev(s.end()) : std::prev(it));
      auto [i2, h2] = *((std::next(it) == s.end()) ? s.begin() : std::next(it));

      mesh.AddFace({mesh.HStart(h0), mesh.HStart(h2), mesh.HStart(h1)});
      sum_angle[i0] += std::acos(std::clamp(
          mesh.HGeometry(mesh.NumHalfedges() - 3)
              .normalized()
              .dot(-mesh.HGeometry(mesh.NumHalfedges() - 1).normalized()),
          -1.0, 1.0));
      q.push(std::make_pair(i0, sum_angle[i0]));
      sum_angle[i2] += std::acos(std::clamp(
          mesh.HGeometry(mesh.NumHalfedges() - 2)
              .normalized()
              .dot(-mesh.HGeometry(mesh.NumHalfedges() - 3).normalized()),
          -1.0, 1.0));
      q.push(std::make_pair(i2, sum_angle[i2]));
      s.erase(it);
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

void Processing::PrintabilityHeuristics(Trimesh& mesh, int niters) {
  if (mesh.NumFaces() == 0) {
    return;
  }
  RetainLargestComponent(mesh);
  std::vector<double> lengths(mesh.NumHalfedges());
  for (H h : mesh.Halfedges()) {
    lengths[h] = mesh.HLength(h);
  }
  std::sort(lengths.begin(), lengths.end());
  double max_length = 2.0 * lengths[lengths.size() / 2];
  double min_length = 0.5 * lengths[lengths.size() / 2];
  for (int i = 0; i < niters; i++) {
    ClampEdgeLengths(mesh, min_length, max_length);
    TL::Processing::RemoveSelfIntersections(mesh);
    mesh.DisconnectFacesUntilManifold();
    RetainLargestComponent(mesh);
    TL::Processing::FillMeshHoles(mesh, 0);
    mesh.DisconnectFacesUntilManifold();
    RetainLargestComponent(mesh);
    MakeDelaunay(mesh);
    TL::Processing::TaubinSmoothing(mesh);
  }
}
}  // namespace TL

#endif  // MESH_PROCESSING_HPP