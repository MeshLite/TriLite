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

#ifndef TRILITE_PROCESSING_HPP
#define TRILITE_PROCESSING_HPP

#include <set>

#include "../Core/Trimesh.hpp"

namespace TL {

class Processing {
 public:
  /**
   * @brief External decimation function of a triangular mesh by edge length
   * ordering.
   * @param mesh Reference to the Trimesh object to be simplified.
   * @param target_face_count The target number of faces after decimation.
   */
  static void Decimate(Trimesh& mesh, size_t target_face_count);

  /**
   * @brief Simplifies a given triangular mesh by collapsing edges based on a
   * quadric error metric. Only valid collapses are performed by default.
   * @param mesh The triangular mesh to be simplified.
   * @param simplification_ratio The ratio between the number of faces after vs
   * before the mesh simplification (expected between 0 and 1).
   * @param preserve_boundaries Whether the collapse are prevented on the
   * boundaries.
   * @param prevent_invalid_collapse Whether a collapse is prevented if it leads
   * to an invalid triangle mesh.
   */
  static void Simplify(Trimesh& mesh, double simplification_ratio = 0.1,
                       bool preserve_boundaries = false,
                       bool prevent_invalid_collapse = true);

  /**
   * @brief Fills holes in a triangular mesh by detecting boundary edges and
   * adding triangles.
   * @param mesh The triangular mesh to be processed.
   */
  static void FillHoles(Trimesh& mesh, size_t target_hole_count = 0);

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
  static void MakeWatertight(Trimesh& mesh, int niters = 10);

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
        edge_halfedges[h] = ToVector<H>(mesh.EdgeHalfedges(h));
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
          mesh.CollapseEdge(h, mesh.HCentroid(h));
          flag = false;
          break;
        }
      }
      f += flag;
    }
  }
  static std::pair<Vector3d, double> CalculateCircumsphere(const Vector3d& a,
                                                           const Vector3d& b,
                                                           const Vector3d& c) {
    Vector3d ac = c - a;
    Vector3d ab = b - a;
    Vector3d abxac = ab.cross(ac);

    Vector3d to_circumsphere_center = (abxac.cross(ab) * ac.squaredNorm() +
                                       ac.cross(abxac) * ab.squaredNorm()) /
                                      (2.0 * abxac.squaredNorm());
    double circumsphere_radius = to_circumsphere_center.norm();
    Vector3d ccs = a + to_circumsphere_center;

    return std::make_pair(ccs, circumsphere_radius);
  }
  static bool IsDelaunay(const Trimesh& mesh, H h) {
    H hopp = mesh.HOpposite(h);
    if (hopp == kInvalidId) {
      return true;
    }
    Vector3d a = mesh.VPosition(mesh.HStart(h));
    Vector3d b = mesh.VPosition(mesh.HEnd(h));
    Vector3d c = mesh.VPosition(mesh.HEnd(mesh.HNext(h)));
    Vector3d d = mesh.VPosition(mesh.HEnd(mesh.HNext(hopp)));
    std::pair<Vector3d, double> sphere = CalculateCircumsphere(a, b, c);
    return (sphere.first - d).norm() >= (sphere.second);
  }
  static void MakeDelaunay(Trimesh& mesh) {
    std::queue<H> edge_queue;
    for (H h : mesh.Halfedges()) {
      if (mesh.HOpposite(h) != kInvalidId &&
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
      H h = edge_queue.front();
      edge_queue.pop();
      H hopp = mesh.HOpposite(h);
      if (hopp != kInvalidId && !IsDelaunay(mesh, h)) {
        V v1 = mesh.HStart(mesh.HPrev(h));
        V v2 = mesh.HStart(mesh.HPrev(hopp));
        bool valid = true;
        for (H he : mesh.VStartings(v1)) {
          valid &= (mesh.HEnd(he) != v2);
        }
        for (H he : mesh.VStartings(v2)) {
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
  static void RetainLargestComponent(Trimesh& mesh) {
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
    std::optional<std::pair<F, std::array<Vector3d, 3>>> triangle_;
    BVHNode* left_;
    BVHNode* right_;

    BVHNode() : left_(nullptr), right_(nullptr) {}
    ~BVHNode() {
      delete left_;
      delete right_;
    }
  };

  static BVHNode* ConstructBVH(
      std::vector<std::pair<F, std::array<Vector3d, 3>>>& triangles,
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
                [axis](const std::pair<F, std::array<Vector3d, 3>>& a,
                       const std::pair<F, std::array<Vector3d, 3>>& b) {
                  return a.second[0][axis] < b.second[0][axis];
                });

      size_t mid = start + (end - start) / 2;
      node->left_ = ConstructBVH(triangles, start, mid, depth + 1);
      node->right_ = ConstructBVH(triangles, mid, end, depth + 1);
    }
    return node;
  }

  static bool DoesIntersect(BVHNode* node,
                            const std::pair<F, std::array<Vector3d, 3>>& tri) {
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
  static std::vector<F> FindSelfIntersections(Trimesh& mesh,
                                              double shrink_factor = 1e-8) {
    if (mesh.NumFaces() == 0) {
      return {};
    }
    std::vector<std::pair<F, std::array<Vector3d, 3>>> triangles;
    for (F f : mesh.Faces()) {
      if (mesh.FArea(f) > 1e-10) {
        Vector3d centroid = mesh.FCentroid(f);
        std::array<Vector3d, 3> tri;
        int i = 0;
        for (const Vector3d& position : mesh.FPositions(f)) {
          tri[i++] = position + (centroid - position) * shrink_factor;
        }
        triangles.emplace_back(f, tri);
      }
    }
    std::vector<F> self_intersection_faces;
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

void Processing::Decimate(Trimesh& mesh, size_t target_face_count) {
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
    auto [rm_faces, rm_vertices] =
        mesh.CollapseEdge(minh, mesh.HCentroid(minh));
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

void Processing::Simplify(TL::Trimesh& mesh, double simplification_ratio,
                          bool preserve_boundaries,
                          bool prevent_invalid_collapse) {
  using Eigen::Matrix4d;
  using Eigen::Vector4d;
  // double max_collapse_cost =
  //     max_collapse_cost_ratio * std::pow(mesh.MedianEdgeLength(), 2);
  auto compute_vertex_quadric = [](const TL::Trimesh& mesh, TL::V v) {
    Matrix4d Q = Matrix4d::Zero();
    for (TL::F f : mesh.VFaces(v)) {
      Vector3d normal = mesh.FNormal(f);
      Vector3d point_on_plane = mesh.FCentroid(f);
      double d = -normal.dot(point_on_plane);  // Plane offset from origin
      Vector4d plane;
      plane << normal(0), normal(1), normal(2), d;
      Q += plane * plane.transpose();
    }
    return Q;
  };
  auto collapse_cost = [&](TL::H h) {
    TL::V v1 = mesh.HStart(h);
    TL::V v2 = mesh.HEnd(h);
    if (std::max(mesh.VValence(v1), mesh.VValence(v2)) > 20) {
      return std::make_pair(0.0, mesh.HCentroid(h));
    }
    if (preserve_boundaries && (mesh.VIsBoundary(v1) || mesh.VIsBoundary(v2))) {
      return std::make_pair(std::numeric_limits<double>::max(),
                            Vector3d{0.0, 0.0, 0.0});
    }
    Matrix4d Q =
        compute_vertex_quadric(mesh, v1) + compute_vertex_quadric(mesh, v2);
    Matrix4d Q_bar = Q;
    Q_bar(3, 0) = Q_bar(3, 1) = Q_bar(3, 2) = 0;
    Q_bar(3, 3) = 1;

    Vector4d v_bar;
    if (Q_bar.determinant() > 1e-5) {
      v_bar = Q_bar.inverse() * Vector4d(0, 0, 0, 1);
    } else {
      Vector3d midpoint = (mesh.VPosition(v1) + mesh.VPosition(v2)) / 2.0;
      v_bar << midpoint, 1.0;
    }
    if (mesh.VIsBoundary(v1) || mesh.VIsBoundary(v2)) {
      v_bar << (mesh.VIsBoundary(v1) ? mesh.VIsBoundary(v2) ? mesh.HCentroid(h)
                                                            : mesh.VPosition(v1)
                                     : mesh.VPosition(v2)),
          1.0;
    }
    double cost = v_bar.transpose() * Q * v_bar;
    return std::make_pair(cost, Vector3d{v_bar[0], v_bar[1], v_bar[2]});
  };
  auto is_valid_collapse = [&](TL::H h, const Vector3d& p) -> bool {
    if (!prevent_invalid_collapse) {
      return true;
    }
    V v1 = mesh.HStart(h);
    V v2 = mesh.HEnd(h);
    for (int i = 0; i < 2; ++i) {
      std::swap(v1, v2);
      for (F f : mesh.VFaces(v1)) {
        std::vector<Vector3d> triangle;
        for (V v : mesh.FVertices(f)) {
          if (v == v2) {
            break;
          }
          if (v == v1) {
            triangle.push_back(p);
          } else {
            triangle.push_back(mesh.VPosition(v));
          }
        }
        if (triangle.size() == 3) {
          Vector3d n =
              (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
          if (n.dot(mesh.FNormal(f)) <= 0.0 ||
              n.norm() < 1e-2 * mesh.FArea(f)) {
            return false;
          }
        }
      }
    }
    return true;
  };
  H id = (H)(simplification_ratio * mesh.NumHalfedges());
  if (id < mesh.NumHalfedges()) {
    std::vector<double> costs;
    for (H h : mesh.Halfedges()) {
      costs.push_back(collapse_cost(h).first);
    }
    std::nth_element(costs.begin(), costs.begin() + id, costs.end(),
                     std::greater<double>());
    double max_collapse_cost = costs[id];
    for (F st_f = 0; st_f < mesh.NumFaces(); st_f++) {
      std::vector<F> vf{st_f};
      while (!vf.empty() && mesh.NumFaces()) {
        F f = vf.back();
        vf.pop_back();
        if (f <= std::min(st_f, mesh.NumFaces() - 1)) {
          for (H h : mesh.FHalfedges(f)) {
            if (mesh.NumHalfedges() <= id) {
              return;
            }
            auto [cost, midpoint] = collapse_cost(h);
            if (cost <= max_collapse_cost && is_valid_collapse(h, midpoint)) {
              auto [deleted_faces, deleted_verts] =
                  mesh.CollapseEdge(h, midpoint);
              vf.insert(vf.end(), deleted_faces.begin(), deleted_faces.end());
              break;
            }
          }
        }
      }
    }
  }
}

void Processing::FillHoles(Trimesh& mesh, size_t target_hole_count) {
  std::unordered_set<H> boundary_edges;
  std::vector<std::vector<H>> polygons;
  for (H st_h : mesh.BoundaryHalfedges()) {
    if (!boundary_edges.count(st_h)) {
      polygons.push_back(ToVector<H>(mesh.HHalfedgesAroundHole(st_h)));
      boundary_edges.insert(polygons.back().begin(), polygons.back().end());
    }
  }
  std::sort(polygons.begin(), polygons.end(),
            [](auto& A, auto& B) { return A.size() < B.size(); });
  size_t id = 0;
  for (const auto& polygon : polygons) {
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
    id++;
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
  std::vector<F> rm_faces;
  for (F f : FindSelfIntersections(mesh)) {
    rm_faces.push_back(f);
  }
  mesh.RemoveFaces(rm_faces);
}

void Processing::MakeWatertight(Trimesh& mesh, int niters) {
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
  double noise_lambda{1e-2 * lengths[lengths.size() / 2]};
  for (int i = 0; i < niters; i++) {
    ClampEdgeLengths(mesh, min_length, max_length);
    for (Vector3d& positions : mesh.Positions()) {
      positions += noise_lambda * Vector3d::Random();
    }
    Processing::RemoveSelfIntersections(mesh);
    mesh.DisconnectFacesUntilManifold();
    RetainLargestComponent(mesh);
    Processing::FillHoles(mesh, 0);
    mesh.DisconnectFacesUntilManifold();
    RetainLargestComponent(mesh);
    MakeDelaunay(mesh);
    Processing::TaubinSmoothing(mesh);
  }
}
}  // namespace TL

#endif  // TRILITE_PROCESSING_HPP