// The MIT License (MIT)
//
// Copyright (c) 2021 José Antonio Fernández Fernández
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

#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "../Core/Trimesh.hpp"

namespace TL {

class Distance {
 public:
  /**
   * @brief Computes the asymmetric Hausdorff distance between two triangle
   * meshes.
   *
   * This function computes the asymmetric Hausdorff distance between two
   * triangle meshes, `mesh` and `target_mesh`. The distance is calculated by
   * finding the maximum distance from any point on `mesh` to its closest point
   * on `target_mesh`.
   *
   * @param mesh The first triangle mesh from which points are sampled.
   * @param target_mesh The second triangle mesh to which distances are
   * calculated.
   * @param precision Desired precision for the computed asymmetric Hausdorff
   * distance, influencing the sampling density on the mesh.
   * @return double The computed asymmetric Hausdorff distance.
   */
  static double AsymmetricHausdorff(const Trimesh& mesh,
                                    const Trimesh& target_mesh,
                                    double precision);

  /**
   * @brief Computes the Hausdorff distance between two triangle meshes.
   *
   * This function computes the Hausdorff distance between two triangle meshes,
   * `mesh1` and `mesh2`. The distance is calculated as the maximum of the
   * asymmetric Hausdorff distances computed in both directions: from `mesh1`
   * to `mesh2` and from `mesh2` to `mesh1`.
   *
   * @param mesh1 The first triangle mesh.
   * @param mesh2 The second triangle mesh.
   * @param precision Desired precision for the computed Hausdorff distance,
   * influencing the sampling density on the mesh.
   * @return double The computed Hausdorff distance.
   */
  static double Hausdorff(const Trimesh& mesh1, const Trimesh& mesh2,
                          double precision);
  /**
   * @brief Computes the asymmetric mean Euclidean distance between two triangle
   * meshes.
   *
   * This function computes the asymmetric mean Euclidean distance between two
   * triangle meshes, `mesh1` and `mesh2`. The distance is calculated by
   * sampling points from `mesh` and finding the closest points on
   * `target_mesh`.
   *
   * @param mesh The first triangle mesh from which points are sampled.
   * @param target_mesh The second triangle mesh to which distances are
   * calculated.
   * @param precision Desired precision for the computed asymmetric mean
   * Euclidean distance, influencing the sampling density on the mesh.
   * @return double The computed asymmetric mean Euclidean distance.
   */
  static double AsymmetricMeanEuclidean(const Trimesh& mesh,
                                        const Trimesh& target_mesh,
                                        double precision);

  /**
   * @brief Computes the mean Euclidean distance between two triangle meshes.
   *
   * This function computes the mean Euclidean distance between two triangle
   * meshes, `mesh1` and `mesh2`. The distance is calculated by averaging the
   * asymmetric mean Euclidean distances computed in both directions: from
   * `mesh1` to `mesh2` and from `mesh2` to `mesh1`.
   *
   * @param mesh1 The first triangle mesh.
   * @param mesh2 The second triangle mesh.
   * @param precision Desired precision for the computed mean Euclidean
   * distance, influencing the sampling density on the mesh.
   * @return double The computed mean Euclidean distance.
   */
  static double MeanEuclidean(const Trimesh& mesh1, const Trimesh& mesh2,
                              double precision);

  class Tree {
   public:
    /**
     * @brief Constructs a Tree object with a given triangle mesh.
     *
     * Constructor that initializes the Tree object with a given triangle
     * mesh.
     *
     * @param mesh The triangle mesh used to construct the tree.
     */
    Tree(const Trimesh& mesh);

    /**
     * @brief Computes the unsigned Euclidean distance from a given point to
     * the nearest point on the mesh.
     *
     * This method computes the unsigned Euclidean distance from a given point
     * to the nearest point on the mesh.
     *
     * @param point The point from which the distance to the mesh is
     * calculated.
     * @return double The unsigned distance from the given point to the
     * closest point on the mesh.
     */
    double Distance(const Vector3d& point);

    /**
     * @brief Finds the closest point on the mesh to a given point.
     *
     * This method finds the closest point on the mesh to a given point.
     *
     * @param point The point from which the closest point on the mesh is
     * found.
     * @return Vector3d The closest point on the mesh to the given point.
     */
    Vector3d ClosestPoint(const Vector3d& point);

   private:
    // Point-Triangle distance declarations
    enum class NearestEntity { V0, V1, V2, E01, E12, E02, F };

    /**
     * Computes the squared distance, the nearest entity (vertex, edge or
     * face) and the nearest point from a point to a triangle.
     */
    static void PointTriangleSqUnsigned(double& distance_sq,
                                        NearestEntity& nearestEntity,
                                        Vector3d& barycentric,
                                        Vector3d& nearestPoint,
                                        const Vector3d& p, const Vector3d& a,
                                        const Vector3d& b, const Vector3d& c) {
      // This function is a modified version of the one found in the Real-Time
      // Collision Detection book by Ericson.
      Vector3d ab = b - a;
      Vector3d ac = c - a;
      Vector3d bc = c - b;

      // Compute parametric position s for projection P’ of P on AB
      double snom = (p - a).dot(ab), sdenom = (p - b).dot(a - b);
      // Compute parametric position t for projection P’ of P on AC
      double tnom = (p - a).dot(ac), tdenom = (p - c).dot(a - c);
      if (snom <= 0.0 && tnom <= 0.0) {
        nearestEntity = NearestEntity::V0;
        barycentric = {1.0, 0.0, 0.0};
        nearestPoint = a;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }

      // Compute parametric position u for projection P’ of P on BC
      double unom = (p - b).dot(bc), udenom = (p - c).dot(b - c);
      if (sdenom <= 0.0 && unom <= 0.0) {
        nearestEntity = NearestEntity::V1;
        barycentric = {0.0, 1.0, 0.0};
        nearestPoint = b;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }
      if (tdenom <= 0.0 && udenom <= 0.0) {
        nearestEntity = NearestEntity::V2;
        barycentric = {0.0, 0.0, 1.0};
        nearestPoint = c;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }

      // Normal for the triangle
      Vector3d n = ab.cross(ac);

      // Check if P is outside AB
      double vc = n.dot((a - p).cross(b - p));
      if (vc <= 0.0 && snom >= 0.0 && sdenom >= 0.0) {
        double arc = snom / (snom + sdenom);
        nearestEntity = NearestEntity::E01;
        barycentric = {1.0 - arc, arc, 0.0};
        nearestPoint = barycentric[0] * a + barycentric[1] * b;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }

      // Check if P is outside BC
      double va = n.dot((b - p).cross(c - p));
      if (va <= 0.0 && unom >= 0.0 && udenom >= 0.0) {
        double arc = unom / (unom + udenom);
        nearestEntity = NearestEntity::E12;
        barycentric = {0.0, 1.0 - arc, arc};
        nearestPoint = barycentric[1] * b + barycentric[2] * c;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }

      // Check if P is outside AC
      double vb = n.dot((c - p).cross(a - p));
      if (vb <= 0.0 && tnom >= 0.0 && tdenom >= 0.0) {
        double arc = tnom / (tnom + tdenom);
        nearestEntity = NearestEntity::E02;
        barycentric = {1.0 - arc, 0.0, arc};
        nearestPoint = barycentric[0] * a + barycentric[2] * c;
        distance_sq = (p - nearestPoint).squaredNorm();
        return;
      }

      // P must project inside the triangle; compute using barycentric
      // coordinates
      double u = va / (va + vb + vc);
      double v = vb / (va + vb + vc);
      double w = 1.0 - u - v;  // = vc / (va + vb + vc)
      nearestEntity = NearestEntity::F;
      barycentric = {u, v, w};
      nearestPoint = u * a + v * b + w * c;
      distance_sq = (p - nearestPoint).squaredNorm();
      return;
    }
    // -----------------------------------

    // Struct that contains the result of a distance query
    struct Result {
      double distance_ = std::numeric_limits<double>::max();
      Vector3d nearestPoint_;
      NearestEntity nearestEntity_;
      int triangleId_ = -1;
      Vector3d barycentric_;
    };
    struct BoundingSphere {
      Vector3d center_{0., 0., 0.};
      double radius_;
    };

    struct Node {
      BoundingSphere bvLeft_;
      BoundingSphere bvRight_;
      int left_ = -1;  // If left == -1, right is the triangle_id
      int right_ = -1;
    };

    struct Triangle {
      std::array<Vector3d, 3> vertices_;
      int id_ = -1;
    };

    std::vector<Vector3d> vertices_;
    std::vector<std::array<int, 3>> triangles_;
    std::vector<Node> nodes_;
    std::vector<Vector3d> pseudonormalsTriangles_;
    std::vector<std::array<Vector3d, 3>> pseudonormalsEdges_;
    std::vector<Vector3d> pseudonormalsVertices_;
    BoundingSphere rootBv_;
    bool isConstructed_ = false;

    void Construct();
    void BuildTree(const int node_id, BoundingSphere& bounding_sphere,
                   std::vector<Triangle>& triangles, const int begin,
                   const int end);
    void Query(Result& result, const Node& node, const Vector3d& point) const;

    Tree() = default;

    void Construct(const std::vector<Vector3d>& vertices,
                   const std::vector<std::array<int, 3>>& triangles);

    // Result signed_distance(const Vector3d& point) const;

    Result UnsignedResult(const Vector3d& point) const;
  };
};  // namespace

double Distance::AsymmetricHausdorff(const Trimesh& mesh,
                                     const Trimesh& target_mesh,
                                     double precision) {
  Tree tree(target_mesh);
  double hausdorff_d = 0;
  std::queue<std::pair<std::array<Vector3d, 3>, double>> q;
  for (F f : mesh.Faces()) {
    q.push(std::make_pair(
        std::array<Vector3d, 3>{mesh.VPosition(mesh.HStart(3 * f)),
                                mesh.VPosition(mesh.HStart(3 * f + 1)),
                                mesh.VPosition(mesh.HStart(3 * f + 2))},
        std::numeric_limits<double>::max()));
  }
  while (!q.empty()) {
    auto [tri, cur_max] = q.front();
    q.pop();
    if (cur_max < hausdorff_d + precision) {
      continue;
    }
    Vector3d bary = (tri[0] + tri[1] + tri[2]) / 3.;
    double bary_dist = tree.Distance(bary);
    hausdorff_d = std::max(hausdorff_d, bary_dist);
    for (int i : {0, 1, 2}) {
      double max_feasible = bary_dist + (tri[i] - bary).norm();
      if (max_feasible > hausdorff_d + precision) {
        Vector3d a = (tri[0] + tri[1]) / 2.0;
        Vector3d b = (tri[1] + tri[2]) / 2.0;
        Vector3d c = (tri[2] + tri[0]) / 2.0;
        for (std::array<Vector3d, 3> sub_tri :
             std::array<std::array<Vector3d, 3>, 4>{
                 std::array<Vector3d, 3>{a, b, c},
                 std::array<Vector3d, 3>{tri[0], a, c},
                 std::array<Vector3d, 3>{tri[1], b, a},
                 std::array<Vector3d, 3>{tri[2], c, b}}) {
          q.push(std::make_pair(std::move(sub_tri), max_feasible));
        }
        break;
      }
    }
  }
  return hausdorff_d;
}
double Distance::Hausdorff(const Trimesh& mesh1, const Trimesh& mesh2,
                           double precision) {
  return std::max(AsymmetricHausdorff(mesh1, mesh2, precision),
                  AsymmetricHausdorff(mesh2, mesh1, precision));
}
double Distance::AsymmetricMeanEuclidean(const Trimesh& mesh,
                                         const Trimesh& target_mesh,
                                         double precision) {
  // 18 / (2\sqrt{3} + \ln(\sqrt{3} + 2));
  const double ratio_precision_to_length = 3.764855876524;
  double distance = 0;
  double sum_area = 0;
  Tree distance_tree(target_mesh);
  double sq_length = pow(precision * ratio_precision_to_length, 2);
  std::function<void(std::array<Vector3d, 3>&&, double)> process_triangle =
      [&](std::array<Vector3d, 3>&& tri, double area) {
        if (std::max({(tri[1] - tri[0]).squaredNorm(),
                      (tri[2] - tri[1]).squaredNorm(),
                      (tri[0] - tri[2]).squaredNorm()}) <= sq_length) {
          distance +=
              area * distance_tree.Distance((tri[0] + tri[1] + tri[2]) / 3.0);
          return;
        }
        area /= 4.0;
        Vector3d a = (tri[0] + tri[1]) / 2.0;
        Vector3d b = (tri[1] + tri[2]) / 2.0;
        Vector3d c = (tri[2] + tri[0]) / 2.0;
        for (std::array<Vector3d, 3> sub_tri :
             std::array<std::array<Vector3d, 3>, 4>{
                 std::array<Vector3d, 3>{a, b, c},
                 std::array<Vector3d, 3>{tri[0], a, c},
                 std::array<Vector3d, 3>{tri[1], b, a},
                 std::array<Vector3d, 3>{tri[2], c, b}}) {
          process_triangle(std::move(sub_tri), area);
        }
      };

  for (TL::F f : mesh.Faces()) {
    std::array<Vector3d, 3> tri{mesh.VPosition(mesh.HStart(3 * f)),
                                mesh.VPosition(mesh.HStart(3 * f + 1)),
                                mesh.VPosition(mesh.HStart(3 * f + 2))};
    double area = mesh.FArea(f);
    process_triangle(std::move(tri), area);
    sum_area += area;
  }
  return distance / sum_area;
}
double Distance::MeanEuclidean(const Trimesh& mesh1, const Trimesh& mesh2,
                               double precision) {
  double distance1 = Distance::AsymmetricMeanEuclidean(mesh1, mesh2, precision);
  double distance2 = Distance::AsymmetricMeanEuclidean(mesh2, mesh1, precision);
  return (distance1 + distance2) / 2.0;
}
Distance::Tree::Tree(const Trimesh& mesh) {
  std::vector<Vector3d> vertices;
  for (const Vector3d& p : mesh.Positions()) {
    vertices.push_back({p[0], p[1], p[2]});
  }
  std::vector<std::array<int, 3>> triangles(mesh.NumFaces());
  for (F f : mesh.Faces()) {
    for (auto [i, v] : std::views::enumerate(mesh.FVertices(f))) {
      triangles[f][i] = v;
    }
  }
  this->Construct(vertices, triangles);
}
double Distance::Tree::Distance(const Vector3d& point) {
  return UnsignedResult(point).distance_;
}
Vector3d Distance::Tree::ClosestPoint(const Vector3d& point) {
  Vector3d p = UnsignedResult(point).nearestPoint_;
  return {p[0], p[1], p[2]};
}
void Distance::Tree::Query(Result& result, const Node& node,
                           const Vector3d& point) const {
  // End of recursion
  if (node.left_ == -1) {
    const int triangle_id = node.right_;
    const std::array<int, 3>& triangle =
        this->triangles_[node.right_];  // If left == -1, right is the
                                        // triangle_id
    const Vector3d& v0 = this->vertices_[triangle[0]];
    const Vector3d& v1 = this->vertices_[triangle[1]];
    const Vector3d& v2 = this->vertices_[triangle[2]];

    double distance_sq;
    NearestEntity nearestEntity;
    Vector3d barycentric;
    Vector3d nearestPoint;
    PointTriangleSqUnsigned(distance_sq, nearestEntity, barycentric,
                            nearestPoint, point, v0, v1, v2);

    if (distance_sq < result.distance_ * result.distance_) {
      result.nearestEntity_ = nearestEntity;
      result.nearestPoint_ = nearestPoint;
      result.barycentric_ = barycentric;
      result.distance_ = std::sqrt(distance_sq);
      result.triangleId_ = triangle_id;
    }
  }

  // Recursion
  else {
    // Find which child bounding volume is closer
    const double d_left =
        (point - node.bvLeft_.center_).norm() - node.bvLeft_.radius_;
    const double d_right =
        (point - node.bvRight_.center_).norm() - node.bvRight_.radius_;

    if (d_left < d_right) {
      // Overlap test
      if (d_left < result.distance_) {
        this->Query(result, this->nodes_[node.left_], point);
      }

      if (d_right < result.distance_) {
        this->Query(result, this->nodes_[node.right_], point);
      }
    } else {
      if (d_right < result.distance_) {
        this->Query(result, this->nodes_[node.right_], point);
      }
      if (d_left < result.distance_) {
        this->Query(result, this->nodes_[node.left_], point);
      }
    }
  }
}
void Distance::Tree::Construct(
    const std::vector<Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles) {
  this->vertices_.resize(vertices.size());
  for (size_t i = 0; i < vertices.size(); i++) {
    this->vertices_[i][0] = (double)vertices[i][0];
    this->vertices_[i][1] = (double)vertices[i][1];
    this->vertices_[i][2] = (double)vertices[i][2];
  }
  this->triangles_.resize(triangles.size());
  for (size_t i = 0; i < triangles.size(); i++) {
    this->triangles_[i][0] = (int)triangles[i][0];
    this->triangles_[i][1] = (int)triangles[i][1];
    this->triangles_[i][2] = (int)triangles[i][2];
  }
  this->Construct();
}

// Result Tree::signed_distance(const Vector3d& point) const
// {
//   const Vector3d p(point[0], point[1], point[2]);
//   Result result = this->UnsignedResult(point);

//   const std::array<int, 3>& triangle =
//   this->triangles_[result.triangle_id]; Vector3d pseudonormal; switch
//   (result.nearestEntity) {
//     case NearestEntity::V0:
//       pseudonormal = this->pseudonormalsVertices[triangle[0]];
//       break;
//     case NearestEntity::V1:
//       pseudonormal = this->pseudonormalsVertices[triangle[1]];
//       break;
//     case NearestEntity::V2:
//       pseudonormal = this->pseudonormalsVertices[triangle[2]];
//       break;
//     case NearestEntity::E01:
//       pseudonormal = this->pseudonormalsEdges[result.triangle_id][0];
//       break;
//     case NearestEntity::E12:
//       pseudonormal = this->pseudonormalsEdges[result.triangle_id][1];
//       break;
//     case NearestEntity::E02:
//       pseudonormal = this->pseudonormalsEdges[result.triangle_id][2];
//       break;
//     case NearestEntity::F:
//       pseudonormal = this->pseudonormalsTriangles[result.triangle_id];
//       break;

//     default:
//       break;
//   }

//   const Vector3d nearestPoint(
//       result.barycentric[0] * this->vertices_[triangle[0]] +
//       result.barycentric[1] * this->vertices_[triangle[1]] +
//       result.barycentric[2] * this->vertices_[triangle[2]]);
//   const Vector3d u = p - nearestPoint;
//   result.distance_ *= (u.dot(pseudonormal) >= 0.0) ? 1.0 : -1.0;

//   return result;
// }

Distance::Tree::Result Distance::Tree::UnsignedResult(
    const Vector3d& point) const {
  if (!this->isConstructed_) {
    std::cout << "DistanceTriangleMesh error: not constructed." << std::endl;
    exit(-1);
  }

  const Vector3d p(point[0], point[1], point[2]);
  Result result;
  result.distance_ = std::numeric_limits<double>::max();
  this->Query(result, this->nodes_[0], p);
  return result;
}

void Distance::Tree::Construct() {
  if (this->triangles_.size() == 0) {
    std::cout << "DistanceTriangleMesh error: Empty triangle list."
              << std::endl;
    exit(-1);
  }

  // Build the tree containing the triangles
  std::vector<Triangle> triangles;

  triangles.resize(this->triangles_.size());
  for (int i = 0; i < (int)this->triangles_.size(); i++) {
    triangles[i].id_ = i;

    const std::array<int, 3>& triangle = this->triangles_[i];
    triangles[i].vertices_[0] = this->vertices_[triangle[0]];
    triangles[i].vertices_[1] = this->vertices_[triangle[1]];
    triangles[i].vertices_[2] = this->vertices_[triangle[2]];
  }

  this->nodes_.push_back(Node());
  this->BuildTree(0, this->rootBv_, triangles, 0, (int)triangles.size());

  // Compute pseudonormals
  //// Edge data structure
  std::unordered_map<uint64_t, Vector3d> edge_normals;
  std::unordered_map<uint64_t, int> edges_count;
  const uint64_t n_vertices = (uint64_t)this->vertices_.size();
  auto add_edge_normal = [&](const int i, const int j,
                             const Vector3d& triangle_normal) {
    const uint64_t key = std::min(i, j) * n_vertices + std::max(i, j);
    if (edge_normals.find(key) == edge_normals.end()) {
      edge_normals[key] = triangle_normal;
      edges_count[key] = 1;
    } else {
      edge_normals[key] += triangle_normal;
      edges_count[key] += 1;
    }
  };
  auto get_edge_normal = [&](const int i, const int j) {
    const uint64_t key = std::min(i, j) * n_vertices + std::max(i, j);
    return edge_normals.find(key)->second;
  };

  //// Compute
  this->pseudonormalsTriangles_.resize(this->triangles_.size());
  this->pseudonormalsEdges_.resize(this->triangles_.size());
  this->pseudonormalsVertices_.resize(this->vertices_.size(), {0, 0, 0});
  for (int i = 0; i < (int)this->triangles_.size(); i++) {
    // Triangle
    const std::array<int, 3>& triangle = this->triangles_[i];
    const Vector3d& a = this->vertices_[triangle[0]];
    const Vector3d& b = this->vertices_[triangle[1]];
    const Vector3d& c = this->vertices_[triangle[2]];

    const Vector3d triangle_normal = (b - a).cross(c - a).normalized();
    this->pseudonormalsTriangles_[i] = triangle_normal;

    // Vertex
    const double alpha_0 =
        std::acos((b - a).normalized().dot((c - a).normalized()));
    const double alpha_1 =
        std::acos((a - b).normalized().dot((c - b).normalized()));
    const double alpha_2 =
        std::acos((b - c).normalized().dot((a - c).normalized()));
    this->pseudonormalsVertices_[triangle[0]] += alpha_0 * triangle_normal;
    this->pseudonormalsVertices_[triangle[1]] += alpha_1 * triangle_normal;
    this->pseudonormalsVertices_[triangle[2]] += alpha_2 * triangle_normal;

    // Edge
    add_edge_normal(triangle[0], triangle[1], triangle_normal);
    add_edge_normal(triangle[1], triangle[2], triangle_normal);
    add_edge_normal(triangle[0], triangle[2], triangle_normal);
  }

  for (Vector3d& n : this->pseudonormalsVertices_) {
    n.normalize();
  }

  for (int tri_i = 0; tri_i < (int)this->triangles_.size(); tri_i++) {
    const std::array<int, 3>& triangle = this->triangles_[tri_i];
    this->pseudonormalsEdges_[tri_i][0] =
        get_edge_normal(triangle[0], triangle[1]).normalized();
    this->pseudonormalsEdges_[tri_i][1] =
        get_edge_normal(triangle[1], triangle[2]).normalized();
    this->pseudonormalsEdges_[tri_i][2] =
        get_edge_normal(triangle[0], triangle[2]).normalized();
  }

  this->isConstructed_ = true;
}

void Distance::Tree::BuildTree(const int node_id,
                               BoundingSphere& bounding_sphere,
                               std::vector<Triangle>& triangles,
                               const int begin, const int end) {
  const int n_triangles = end - begin;

  if (n_triangles == 0) {
    std::cout << "DistanceTriangleMesh::Construct error: Empty leave."
              << std::endl;
    exit(-1);
  } else if (n_triangles == 1) {
    // Build node leaf
    this->nodes_[node_id].left_ = -1;
    this->nodes_[node_id].right_ = triangles[begin].id_;

    //// Bounding sphere
    const Triangle& tri = triangles[begin];
    const Vector3d center =
        (tri.vertices_[0] + tri.vertices_[1] + tri.vertices_[2]) / 3.0;
    const double radius = std::max(std::max((tri.vertices_[0] - center).norm(),
                                            (tri.vertices_[1] - center).norm()),
                                   (tri.vertices_[2] - center).norm());
    bounding_sphere.center_ = center;
    bounding_sphere.radius_ = radius;
  } else {
    // Compute AxisAligned Bounding Box center and largest dimension of all
    // current triangles
    Vector3d top = {std::numeric_limits<double>::lowest(),
                    std::numeric_limits<double>::lowest(),
                    std::numeric_limits<double>::lowest()};
    Vector3d bottom = {std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max()};
    Vector3d center = {0, 0, 0};
    for (int tri_i = begin; tri_i < end; tri_i++) {
      for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
        const Vector3d& p = triangles[tri_i].vertices_[vertex_i];
        center += p;

        for (int coord_i = 0; coord_i < 3; coord_i++) {
          top[coord_i] = std::max(top[coord_i], p[coord_i]);
          bottom[coord_i] = std::min(bottom[coord_i], p[coord_i]);
        }
      }
    }
    center /= 3 * n_triangles;
    const Vector3d diagonal = top - bottom;
    const int split_dim =
        (int)(std::max_element(&diagonal[0], &diagonal[0] + 3) - &diagonal[0]);

    // Set node bounding sphere
    double radius_sq = 0.0;
    for (int tri_i = begin; tri_i < end; tri_i++) {
      for (int i = 0; i < 3; i++) {
        radius_sq = std::max(
            radius_sq, (center - triangles[tri_i].vertices_[i]).squaredNorm());
      }
    }
    bounding_sphere.center_ = center;
    bounding_sphere.radius_ = std::sqrt(radius_sq);

    // Sort the triangles according to their center along the split
    // dimension
    std::sort(triangles.begin() + begin, triangles.begin() + end,
              [split_dim](const Triangle& a, const Triangle& b) {
                return a.vertices_[0][split_dim] < b.vertices_[0][split_dim];
              });

    // Children
    const int mid = (int)(0.5 * (begin + end));

    this->nodes_[node_id].left_ = (int)this->nodes_.size();
    this->nodes_.push_back(Node());
    this->BuildTree(this->nodes_[node_id].left_, this->nodes_[node_id].bvLeft_,
                    triangles, begin, mid);

    this->nodes_[node_id].right_ = (int)this->nodes_.size();
    this->nodes_.push_back(Node());
    this->BuildTree(this->nodes_[node_id].right_,
                    this->nodes_[node_id].bvRight_, triangles, mid, end);
  }
}
}  // namespace TL