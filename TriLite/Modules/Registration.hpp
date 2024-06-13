// The MIT License (MIT)
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

#ifndef TRILITE_INTERPOLATION_HPP
#define TRILITE_INTERPOLATION_HPP

#include <Eigen/Eigen>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

#include "../Core/Trimesh.hpp"
#include "../Core/nanoflann.hpp"

namespace TL {

class Registration {
 public:
  /**
   * @brief Finds the closest points in the target mesh for each point in the
   * source mesh.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @return A vector of indices of the closest points in the target mesh, and a
   * vector containing the corresponding distances.
   */
  static std::pair<std::vector<int>, std::vector<double>> FindClosestPoints(
      const Trimesh& source, const Trimesh& target);

  /**
   * @brief Applies the Kabsch algorithm to find the best rotation matrix
   * between source and target points.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @param already_centered Whether the points are already centered.
   * @return The optimal rotation matrix.
   */
  static Eigen::Matrix3d FindBestRotation(const Trimesh& source,
                                          const Trimesh& target,
                                          bool already_centered = false);

  /**
   * @brief Finds the best rigid transformation (rotation and translation)
   * between source and target meshes.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @return A tuple containing the rotation matrix and translation vector.
   */
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d>
  FindBestRigidTransformation(const Trimesh& source, const Trimesh& target);

  /**
   * @brief Finds the best similarity transformation (rotation, translation, and
   * scaling) between source and target meshes.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @return A tuple containing the rotation matrix, translation vector, and
   * scaling factor.
   */
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d, double>
  FindBestSimilarityTransformation(const Trimesh& source,
                                   const Trimesh& target);

  /**
   * @brief Applies the Iterative Closest Point (ICP) algorithm to find the best
   * alignment between source and target meshes.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @param max_iterations The maximum number of iterations.
   * @param tolerance The convergence tolerance.
   * @return A tuple containing the rotation matrix and translation vector.
   */
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d> ICP(
      const Trimesh& source, const Trimesh& target, int max_iterations = 1000,
      double tolerance = 1e-6);

  /**
   * @brief Applies the Iterative Closest Point (ICP) algorithm to find the best
   * alignment between source and target meshes with different preprocessing on
   * the source mesh to select the best result.
   * @param source The source Trimesh.
   * @param target The target Trimesh.
   * @return A tuple containing the rotation matrix and translation vector.
   */
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d>
  RigidRegistrationHeuristics(const Trimesh& source, const Trimesh& target);

 private:
  struct PointCloudAdaptor {
    const Eigen::MatrixXd& pointCloud;

    PointCloudAdaptor(const Eigen::MatrixXd& pointCloud)
        : pointCloud(pointCloud) {}

    inline size_t kdtree_get_point_count() const { return pointCloud.cols(); }

    inline double kdtree_distance(const double* p1, const size_t idx_p2,
                                  size_t /*size*/) const {
      return (pointCloud.col(idx_p2) -
              Eigen::Map<const Eigen::VectorXd>(p1, pointCloud.rows()))
          .squaredNorm();
    }

    inline double kdtree_get_pt(const size_t idx, int dim) const {
      return pointCloud(dim, idx);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const {
      return false;
    }
  };

  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor>,
      PointCloudAdaptor, -1 /* dim */
      >
      KDTree;

  static void FindClosestPoints(const Eigen::MatrixXd& source,
                                const Eigen::MatrixXd& target,
                                std::vector<int>& indices,
                                std::vector<double>& distances) {
    PointCloudAdaptor pcAdaptor(target);
    KDTree kdTree(target.rows(), pcAdaptor,
                  nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    indices.resize(source.cols());
    distances.resize(source.cols());

    for (Eigen::Index i = 0; i < source.cols(); ++i) {
      std::vector<unsigned int> ret_index(1);
      std::vector<double> out_dist_sqr(1);

      kdTree.knnSearch(source.col(i).data(), 1, &ret_index[0],
                       &out_dist_sqr[0]);

      indices[i] = ret_index[0];
      distances[i] = out_dist_sqr[0];
    }
  }
  // Kabsch algorithm to find the best rotation matrix
  static Eigen::Matrix3d FindBestRotation(const Eigen::MatrixXd& source,
                                          const Eigen::MatrixXd& target,
                                          bool already_centered = false) {
    if (source.cols() != target.cols()) {
      throw std::invalid_argument(
          "Source and target must have the same number of points");
    }
    Eigen::MatrixXd covariance_matrix;
    if (!already_centered) {
      // Compute centroids
      Eigen::Vector3d centroid_source = source.rowwise().mean();
      Eigen::Vector3d centroid_target = target.rowwise().mean();

      // Center the points and compute the covariance matrix
      covariance_matrix = (target.colwise() - centroid_target) *
                          (source.colwise() - centroid_source).transpose();
    } else {
      covariance_matrix = target * source.transpose();
    }

    // Compute the SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        covariance_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d rotation_matrix = svd.matrixU() * svd.matrixV().transpose();

    // Ensure a proper rotation (no reflection)
    if (rotation_matrix.determinant() < 0) {
      Eigen::Matrix3d V = svd.matrixV();
      V.col(2) *= -1;
      rotation_matrix = svd.matrixU() * V.transpose();
    }

    return rotation_matrix;
  }

  // Find the best rotation and translation using the Kabsch algorithm
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d>
  FindBestRigidTransformation(const Eigen::MatrixXd& source,
                              const Eigen::MatrixXd& target) {
    if (source.cols() != target.cols()) {
      throw std::invalid_argument(
          "Source and target must have the same number of points");
    }
    // Compute centroids
    Eigen::Vector3d centroid_source = source.rowwise().mean();
    Eigen::Vector3d centroid_target = target.rowwise().mean();

    Eigen::Matrix3d rotation_matrix =
        FindBestRotation(source.colwise() - centroid_source,
                         target.colwise() - centroid_target, true);

    // Compute the translation
    Eigen::Vector3d translation_vector =
        centroid_target - rotation_matrix * centroid_source;

    return {rotation_matrix, translation_vector};
  }

  // Find the best rotation, translation, and scaling
  static std::tuple<Eigen::Matrix3d, Eigen::Vector3d, double>
  FindBestSimilarityTransformation(const Eigen::MatrixXd& source,
                                   const Eigen::MatrixXd& target) {
    if (source.cols() != target.cols()) {
      throw std::invalid_argument(
          "Source and target must have the same number of points");
    }

    // Compute centroids
    Eigen::Vector3d centroid_source = source.rowwise().mean();
    Eigen::Vector3d centroid_target = target.rowwise().mean();

    // Center the points
    Eigen::MatrixXd centered_source = source.colwise() - centroid_source;
    Eigen::MatrixXd centered_target = target.colwise() - centroid_target;

    // Compute the scale
    double norm_source = centered_source.norm();
    double norm_target = centered_target.norm();
    double scale = norm_target / norm_source;

    // Scale the source points
    centered_source *= scale;

    // Compute the rotation using the Kabsch algorithm
    Eigen::Matrix3d rotation_matrix =
        FindBestRotation(centered_source, centered_target, true);

    // Compute the translation
    Eigen::Vector3d translation_vector =
        centroid_target - scale * (rotation_matrix * centroid_source);

    return {rotation_matrix, translation_vector, scale};
  }
  static std::pair<Eigen::Matrix3d, Vector3d> ICP(const Eigen::MatrixXd& source,
                                                  const Eigen::MatrixXd& target,
                                                  int max_iterations = 1000,
                                                  double tolerance = 1e-6) {
    Eigen::MatrixXd current_source = source;
    Eigen::Matrix4d total_transformation = Eigen::Matrix4d::Identity();
    std::vector<int> indices;
    std::vector<double> distances;

    for (int iter = 0; iter < max_iterations; ++iter) {
      FindClosestPoints(current_source, target, indices, distances);

      Eigen::MatrixXd tmp_target(target.rows(), source.cols());

      for (Eigen::Index i = 0; i < source.cols(); ++i) {
        tmp_target.col(i) = target.col(indices[i]);
      }
      auto [rotation, translation] =
          FindBestRigidTransformation(current_source, tmp_target);

      current_source = (rotation * current_source).colwise() + translation;
      Eigen::Matrix4d transformation = Eigen::Matrix4d::Identity();
      transformation.block<3, 3>(0, 0) = rotation;
      transformation.block<3, 1>(0, 3) = translation;
      total_transformation = transformation * total_transformation;

      double diff = (transformation - Eigen::Matrix4d::Identity()).norm();
      if (diff < tolerance) {
        break;
      }
    }

    return std::make_pair(total_transformation.block<3, 3>(0, 0),
                          total_transformation.block<3, 1>(0, 3));
  }
  static std::pair<Eigen::Matrix3d, Vector3d> RigidRegistrationHeuristics(
      const Eigen::MatrixXd& source, const Eigen::MatrixXd& target) {
    Vector3d source_centroid = source.rowwise().mean();
    Vector3d target_centroid = target.rowwise().mean();

    Eigen::MatrixXd source_centered = source.colwise() - source_centroid;
    Eigen::MatrixXd target_centered = target.colwise() - target_centroid;

    double min_distance = std::numeric_limits<double>::max();
    Eigen::Matrix3d best_rotation;
    Vector3d best_translation;

    std::vector<Eigen::Matrix3d> rotations;
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          if ((bool)i + (bool)j + (bool)k == 1) {
            Eigen::Matrix3d rot;
            rot.col(0) = Vector3d(i, j, k).normalized();
            for (int ii = -1; ii <= 1; ii++) {
              for (int jj = -1; jj <= 1; jj++) {
                for (int kk = -1; kk <= 1; kk++) {
                  if ((bool)ii + (bool)jj + (bool)kk == 1 &&
                      (bool)(ii && i) + (bool)(jj && j) + (bool)(kk && k) ==
                          0) {
                    rot.col(1) = Vector3d(ii, jj, kk).normalized();
                    rot.col(2) = rot.col(0).cross(rot.col(1));
                    rotations.push_back(rot);
                  }
                }
              }
            }
          }
        }
      }
    }
    Eigen::Matrix4d final_transformation = Eigen::Matrix4d::Identity();
    for (const Eigen::Matrix3d& R : rotations) {
      Eigen::MatrixXd rotated_source = R * source_centered;
      auto [rotation, translation] = ICP(rotated_source, target_centered);
      Eigen::MatrixXd transformed_source =
          (rotation * rotated_source).colwise() + translation;

      std::vector<int> indices;
      std::vector<double> distances1;
      FindClosestPoints(transformed_source, target_centered, indices,
                        distances1);
      std::vector<double> distances2;
      FindClosestPoints(target_centered, transformed_source, indices,
                        distances2);
      double distance =
          std::accumulate(distances1.begin(), distances1.end(), 0.0) +
          std::accumulate(distances2.begin(), distances2.end(), 0.0);
      if (distance < min_distance) {
        min_distance = distance;
        best_rotation = R * rotation;
        best_translation = translation;
        final_transformation = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
        T1.block<3, 1>(0, 3) = -source_centroid;
        Eigen::Matrix4d R1 = Eigen::Matrix4d::Identity();
        R1.block<3, 3>(0, 0) = R;
        Eigen::Matrix4d A = Eigen::Matrix4d::Identity();
        A.block<3, 3>(0, 0) = rotation;
        A.block<3, 1>(0, 3) = translation;
        Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
        T2.block<3, 1>(0, 3) = target_centroid;
        final_transformation = T2 * A * R1 * T1;
      }
    }

    return std::make_pair(final_transformation.block<3, 3>(0, 0),
                          final_transformation.block<3, 1>(0, 3));
  }
  static Eigen::MatrixXd ConvertTrimeshToMatrix(const Trimesh& mesh) {
    Eigen::MatrixXd matrix(3, mesh.NumVertices());
    for (auto v : mesh.Vertices()) {
      matrix.col(v) = mesh.VPosition(v);
    }
    return matrix;
  }
};
std::pair<std::vector<int>, std::vector<double>>
Registration::FindClosestPoints(const Trimesh& source, const Trimesh& target) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  std::vector<int> indices;
  std::vector<double> distances;
  FindClosestPoints(source_matrix, target_matrix, indices, distances);
  return std::make_pair(indices, distances);
}

Eigen::Matrix3d Registration::FindBestRotation(const Trimesh& source,
                                               const Trimesh& target,
                                               bool already_centered) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  return FindBestRotation(source_matrix, target_matrix, already_centered);
}

std::tuple<Eigen::Matrix3d, Eigen::Vector3d>
Registration::FindBestRigidTransformation(const Trimesh& source,
                                          const Trimesh& target) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  return FindBestRigidTransformation(source_matrix, target_matrix);
}

std::tuple<Eigen::Matrix3d, Eigen::Vector3d, double>
Registration::FindBestSimilarityTransformation(const Trimesh& source,
                                               const Trimesh& target) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  return FindBestSimilarityTransformation(source_matrix, target_matrix);
}

std::tuple<Eigen::Matrix3d, Eigen::Vector3d> Registration::ICP(
    const Trimesh& source, const Trimesh& target, int max_iterations,
    double tolerance) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  return ICP(source_matrix, target_matrix, max_iterations, tolerance);
}

std::tuple<Eigen::Matrix3d, Eigen::Vector3d>
Registration::RigidRegistrationHeuristics(const Trimesh& source,
                                          const Trimesh& target) {
  Eigen::MatrixXd source_matrix = ConvertTrimeshToMatrix(source);
  Eigen::MatrixXd target_matrix = ConvertTrimeshToMatrix(target);
  return RigidRegistrationHeuristics(source_matrix, target_matrix);
}
}  // namespace TL

#endif  // TRILITE_INTERPOLATION_HPP