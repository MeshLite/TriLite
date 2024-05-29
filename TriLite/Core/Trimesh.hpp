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

#ifndef TRILITE_MESH_H
#define TRILITE_MESH_H

#include <Eigen/Eigen>
#include <generator>
#include <memory>
#include <queue>
#include <unordered_set>

namespace TL {

enum MeshComponent { kHalfedge, kVertex, kFace };
using Index = unsigned int;
using H = Index;
using V = Index;
using F = Index;
constexpr Index kInvalidId = std::numeric_limits<Index>::max();
using Eigen::Vector3d;
/**
 * @class Trimesh
 * @brief Represents a triangle mesh using a half-edge data structure.
 */
class Trimesh {
 public:
  /**
   * @brief Constructs an empty Trimesh.
   */
  Trimesh();

  /**
   * @brief Constructs a Trimesh from a vector of vertices and their positions.
   * @param tri_points Vector of vertices and positions defining the triangle
   * mesh.
   */
  Trimesh(const std::vector<std::variant<V, Vector3d>>& tri_points);

  /**
   * @brief Copy constructor for Trimesh.
   * @param other The Trimesh object to copy.
   */
  Trimesh(const Trimesh& other);

  /**
   * @brief Assignment operator for Trimesh.
   * @param other The Trimesh object to assign.
   * @return A reference to the assigned Trimesh object.
   */
  Trimesh& operator=(const Trimesh& other);

  /**
   * @brief Returns the number of halfedges in the mesh.
   * @return The number of halfedges.
   */
  H NumHalfedges() const;

  /**
   * @brief Returns the number of vertices in the mesh.
   * @return The number of vertices.
   */
  V NumVertices() const;

  /**
   * @brief Returns the number of faces in the mesh.
   * @return The number of faces.
   */
  F NumFaces() const;

  /**
   * @brief Returns a view of the halfedges in the mesh.
   * @return A range of halfedge indices.
   */
  std::ranges::iota_view<H, H> Halfedges() const;

  /**
   * @brief Returns a view of the vertices in the mesh.
   * @return A range of vertex indices.
   */
  std::ranges::iota_view<V, V> Vertices() const;

  /**
   * @brief Returns a view of the faces in the mesh.
   * @return A range of face indices.
   */
  std::ranges::iota_view<F, F> Faces() const;

  /**
   * @brief Returns a view of the vertex positions.
   * @return A reference view of the vertex positions.
   */
  std::ranges::ref_view<std::vector<Vector3d>> Positions();

  /**
   * @brief Returns a constant view of the vertex positions.
   * @return A constant reference view of the vertex positions.
   */
  std::ranges::ref_view<const std::vector<Vector3d>> Positions() const;

  /**
   * @brief Returns the next halfedge in the face.
   * @param h The current halfedge.
   * @return The next halfedge in the face.
   */
  H HNext(H h) const;

  /**
   * @brief Returns the previous halfedge in the face.
   * @param h The current halfedge.
   * @return The previous halfedge in the face.
   */
  H HPrev(H h) const;

  /**
   * @brief Returns the starting vertex of the halfedge.
   * @param h The halfedge index.
   * @return The starting vertex index.
   */
  V HStart(H h) const;

  /**
   * @brief Returns the ending vertex of the halfedge.
   * @param h The halfedge index.
   * @return The ending vertex index.
   */
  V HEnd(H h) const;

  /**
   * @brief Returns the face associated with the halfedge.
   * @param h The halfedge index.
   * @return The face index.
   */
  F HFace(H h) const;

  /**
   * @brief Returns the opposite halfedge of the given halfedge.
   * @param h The halfedge index.
   * @return The opposite halfedge index or kInvalidId if it does not exist.
   */
  H HOpposite(H h) const;

  /**
   * @brief Returns the next halfedge around the starting vertex.
   * @param h The current halfedge.
   * @return The next halfedge around the starting vertex.
   */
  H HNextAroundStart(H h) const;

  /**
   * @brief Returns the previous halfedge around the starting vertex.
   * @param h The current halfedge.
   * @return The previous halfedge around the starting vertex.
   */
  H HPrevAroundStart(H h) const;

  /**
   * @brief Returns the next halfedge around the ending vertex.
   * @param h The current halfedge.
   * @return The next halfedge around the ending vertex.
   */
  H HNextAroundEnd(H h) const;

  /**
   * @brief Returns the previous halfedge around the ending vertex.
   * @param h The current halfedge.
   * @return The previous halfedge around the ending vertex.
   */
  H HPrevAroundEnd(H h) const;

  /**
   * @brief Returns the geometric vector of the halfedge.
   * @param h The halfedge index.
   * @return The geometric vector.
   */
  Vector3d HGeometry(H h) const;

  /**
   * @brief Returns a generator for the halfedges connected around the starting
   * vertex.
   * @param st_h The starting halfedge.
   * @return A generator for the connected halfedges.
   */
  std::generator<H> HConnectionsAroundStart(H st_h) const;

  /**
   * @brief Returns a generator for the halfedges around a hole.
   * @param st_h The starting halfedge of the hole.
   * @return A generator for the halfedges around the hole.
   */
  std::generator<H> HHalfedgesAroundHole(H st_h) const;

  /**
   * @brief Returns the length of the halfedge.
   * @param h The halfedge index.
   * @return The length of the halfedge.
   */
  double HLength(H h) const;

  /**
   * @brief Returns the starting halfedge of the vertex.
   * @param v The vertex index.
   * @return The starting halfedge index.
   */
  H VStarting(V v) const;

  /**
   * @brief Returns the ending halfedge of the vertex.
   * @param v The vertex index.
   * @return The ending halfedge index.
   */
  H VEnding(V v) const;

  /**
   * @brief Returns a generator for the starting halfedges of the vertex.
   * @param v The vertex index.
   * @return A generator for the starting halfedges.
   */
  std::generator<H> VStartings(V v) const;

  /**
   * @brief Returns a generator for the ending halfedges of the vertex.
   * @param v The vertex index.
   * @return A generator for the ending halfedges.
   */
  std::generator<H> VEndings(V v) const;

  /**
   * @brief Returns a view of the faces connected to the vertex.
   * @param v The vertex index.
   * @return A view of the face indices.
   */
  auto VFaces(V v) const;

  /**
   * @brief Returns the position of the vertex.
   * @param v The vertex index.
   * @return The position of the vertex.
   */
  Vector3d& VPosition(V v);

  /**
   * @brief Returns the constant position of the vertex.
   * @param v The vertex index.
   * @return The constant position of the vertex.
   */
  const Vector3d& VPosition(V v) const;

  /**
   * @brief Returns the normal vector of the vertex.
   * @param v The vertex index.
   * @return The normal vector.
   */
  Vector3d VNormal(V v) const;

  /**
   * @brief Returns the valence (degree) of the vertex.
   * @param v The vertex index.
   * @return The valence of the vertex.
   */
  size_t VValence(V v) const;

  /**
   * @brief Checks if the vertex is manifold.
   * @param v The vertex index.
   * @return True if the vertex is manifold, false otherwise.
   */
  bool VIsManifold(V v) const;

  /**
   * @brief Returns the halfedge associated with the face.
   * @param f The face index.
   * @return The halfedge index.
   */
  H FHalfedge(F f) const;

  /**
   * @brief Returns a view of the halfedges of the face.
   * @param f The face index.
   * @return A range of halfedge indices.
   */
  std::ranges::iota_view<H, H> FHalfedges(F f) const;

  /**
   * @brief Returns a generator for neighboring faces of a given face.
   * @param f The face index for which to find neighboring faces.
   * @return std::generator<F> A generator for the indices of neighboring faces.
   */
  std::generator<F> FNeighbors(F f) const;

  /**
   * @brief Returns a view of the vertices of the face.
   * @param f The face index.
   * @return A view of the vertex indices.
   */
  auto FVertices(F f) const;

  /**
   * @brief Returns a view of the vertex positions of the face.
   * @param f The face index.
   * @return A view of the vertex positions.
   */
  auto FPositions(F f) const;

  /**
   * @brief Returns the normal vector of the face.
   * @param f The face index.
   * @return The normal vector.
   */
  Vector3d FNormal(F f) const;

  /**
   * @brief Returns the area of the face.
   * @param f The face index.
   * @return The area of the face.
   */
  double FArea(F f) const;

  /**
   * @brief Computes the axis-aligned bounding box (AABB) of a specific face.
   * @param f The face index for which to compute the bounding box.
   * @return A pair of Vector3d, where the first element is the minimum corner
   * of the bounding box and the second element is the maximum corner of the
   * bounding box.
   */
  inline std::pair<Vector3d, Vector3d> FBoundingBox(F f) const;

  /**
   * @brief Computes the centroid (mean) of a specific face.
   * @param f The face index for which to compute the centroid.
   * @return A Vector3d representing the centroid of the specified face.
   */
  inline Vector3d FCentroid(F f) const;

  /**
   * @brief Returns a generator for the halfedges forming an edge.
   * @param h The halfedge index.
   * @return A generator for the halfedges forming the edge.
   */
  std::generator<H> EdgeHalfedges(H h) const;

  /**
   * @brief Returns a view of the faces connected by an edge.
   * @param h The halfedge index.
   * @return A view of the face indices.
   */
  auto EdgeFaces(H h) const;

  /**
   * @brief Checks if the edge is manifold.
   * @param h The halfedge index.
   * @return True if the edge is manifold, false otherwise.
   */
  bool EdgeIsManifold(H h) const;

  /**
   * @brief Returns a view of the boundary halfedges.
   * @return A view of the boundary halfedge indices.
   */
  auto BoundaryHalfedges() const;

  /**
   * @brief Computes the axis-aligned bounding box (AABB) of the entire mesh.
   * @return A pair of Vector3d, where the first element is the minimum corner
   * of the bounding box and the second element is the maximum corner of the
   * bounding box.
   * @throws std::runtime_error if the mesh has no vertices.
   */
  std::pair<Vector3d, Vector3d> BoundingBox() const;

  /**
   * @brief Computes the centroid (mean) of all vertices in the mesh.
   * @return A Vector3d representing the centroid of all vertices in the mesh.
   * @throws std::runtime_error if the mesh has no vertices.
   */
  Vector3d Centroid() const;

  /**
   * @brief Adds a face to the mesh.
   * @param triangle Array of three vertices or positions defining the triangle.
   * @return The face index.
   */
  F AddFace(const std::array<std::variant<V, Vector3d>, 3>& triangle);

  /**
   * @brief Removes a face from the mesh.
   * @param f The face index.
   * @return A vector of vertices removed along with the face.
   */
  std::vector<V> RemoveFace(F f);

  /**
   * @brief Removes multiple face indices from the mesh.
   * @param f_set A vector of face indices to remove (duplicates are ignored).
   * @return A vector of vertices indices removed along with the faces in the
   * order of deletion (vertex indices must be updated accordingly on the fly).
   */
  std::vector<V> RemoveFaces(std::vector<F> f_set);

  /**
   * @brief Collapses an edge and merges its vertices.
   * @param h The halfedge index.
   * @return A pair of vectors containing the removed faces and vertices.
   */
  std::pair<std::vector<F>, std::vector<V>> CollapseEdge(H h);

  /**
   * @brief Disconnects a face from the mesh, converting it to a boundary face.
   * @param f The face index.
   */
  void DisconnectFace(F f);

  /**
   * @brief Disconnects faces until all edges are manifold.
   */
  void DisconnectFacesUntilManifoldEdges();

  /**
   * @brief Disconnects faces until all vertices are manifold.
   */
  void DisconnectFacesUntilManifoldVertices();

  /**
   * @brief Disconnects faces until both edges and vertices are manifold.
   */
  void DisconnectFacesUntilManifold();

  /**
   * @brief Creates an attribute for a specified mesh component.
   * @tparam C The mesh component (halfedge, vertex, or face).
   * @tparam T The attribute type.
   * @param key The key identifying the attribute.
   * @return A reference view of the created attribute.
   * @throws std::invalid_argument if the attribute already exists.
   */
  template <MeshComponent C, typename T>
  std::ranges::ref_view<std::vector<T>> CreateAttribute(const std::string& key);

  /**
   * @brief Retrieves an attribute for a specified mesh component.
   * @tparam C The mesh component (halfedge, vertex, or face).
   * @tparam T The attribute type.
   * @param key The key identifying the attribute.
   * @return A reference view of the attribute.
   * @throws std::invalid_argument if the attribute does not exist or type does
   * not match.
   */
  template <MeshComponent C, typename T>
  std::ranges::ref_view<std::vector<T>> GetAttribute(const std::string& key);

  /**
   * @brief Erases an attribute for a specified mesh component.
   * @tparam C The mesh component (halfedge, vertex, or face).
   * @param key The key identifying the attribute.
   * @throws std::invalid_argument if the attribute does not exist.
   */
  template <MeshComponent C>
  void EraseAttribute(const std::string& key);

 private:
  class BaseContainerWrapper {
   public:
    virtual ~BaseContainerWrapper() = default;
    virtual std::unique_ptr<BaseContainerWrapper> clone() const = 0;
    virtual void ReplaceErase(int i) = 0;
    virtual void IncrementSize() = 0;
  };
  template <typename Container>
  class ContainerWrapper : public BaseContainerWrapper {
   public:
    Container container_;
    explicit ContainerWrapper(const Container& container)
        : container_(container) {}
    std::unique_ptr<BaseContainerWrapper> clone() const override {
      return std::make_unique<ContainerWrapper<Container>>(*this);
    }
    void ReplaceErase(int i) override {
      container_[i] = std::move(container_.back());
      container_.pop_back();
    }
    void IncrementSize() override { container_.emplace_back(); }
  };
  std::vector<V> hStart_;
  std::vector<H> hCoStart_;
  std::vector<H> vStart_;
  std::vector<Vector3d> position_;
  std::array<std::map<std::string, std::unique_ptr<BaseContainerWrapper>>, 3>
      attributes_;
};

Trimesh::Trimesh() {}
Trimesh::Trimesh(const std::vector<std::variant<V, Vector3d>>& tri_points)
    : Trimesh() {
  for (V v = 0; v < (V)tri_points.size(); v += 3) {
    AddFace({tri_points[v], tri_points[v + 1], tri_points[v + 2]});
  }
}
Trimesh::Trimesh(const Trimesh& other)
    : hStart_(other.hStart_),
      hCoStart_(other.hCoStart_),
      vStart_(other.vStart_),
      position_(other.position_) {
  for (int i = 0; i < 3; ++i) {
    for (const auto& [key, wrapper] : other.attributes_[i]) {
      attributes_[i][key] = wrapper->clone();
    }
  }
}
Trimesh& Trimesh::operator=(const Trimesh& other) {
  if (this != &other) {
    hStart_ = other.hStart_;
    hCoStart_ = other.hCoStart_;
    vStart_ = other.vStart_;
    position_ = other.position_;
    for (auto& attr_map : attributes_) {
      attr_map.clear();
    }
    for (int i = 0; i < 3; ++i) {
      for (const auto& [key, wrapper] : other.attributes_[i]) {
        attributes_[i][key] = wrapper->clone();
      }
    }
  }
  return *this;
}
inline H Trimesh::NumHalfedges() const { return hStart_.size(); }
inline V Trimesh::NumVertices() const { return position_.size(); }
inline F Trimesh::NumFaces() const { return NumHalfedges() / 3; }
inline std::ranges::iota_view<H, H> Trimesh::Halfedges() const {
  return std::views::iota(H{0}, H{NumHalfedges()});
}
inline std::ranges::iota_view<V, V> Trimesh::Vertices() const {
  return std::views::iota(V{0}, V{NumVertices()});
}
inline std::ranges::iota_view<F, F> Trimesh::Faces() const {
  return std::views::iota(F{0}, F{NumFaces()});
}
inline std::ranges::ref_view<std::vector<Vector3d>> Trimesh::Positions() {
  return std::views::all(position_);
}
inline std::ranges::ref_view<std::vector<Vector3d> const> Trimesh::Positions()
    const {
  return std::views::all(position_);
}
inline H Trimesh::HNext(H h) const {
  assert(h < NumHalfedges());
  return h - h % 3 + (h + 1) % 3;
}
inline H Trimesh::HPrev(H h) const {
  assert(h < NumHalfedges());
  return h - h % 3 + (h + 2) % 3;
}
inline V Trimesh::HStart(H h) const { return hStart_.at(h); }
inline V Trimesh::HEnd(H h) const { return HStart(HNext(h)); }
inline F Trimesh::HFace(H h) const {
  assert(h < NumHalfedges());
  return h / 3;
}
H Trimesh::HOpposite(H h) const {
  for (H nh : VStartings(HStart(HNext(h)))) {
    if (HEnd(nh) == HStart(h) && HStart(nh) == HEnd(h)) {
      return nh;
    }
  }
  return kInvalidId;
}
inline H Trimesh::HNextAroundStart(H h) const { return HOpposite(HPrev(h)); }
inline H Trimesh::HPrevAroundStart(H h) const {
  return HOpposite(h) == kInvalidId ? kInvalidId : HNext(HOpposite(h));
}
inline H Trimesh::HNextAroundEnd(H h) const {
  return HOpposite(h) == kInvalidId ? kInvalidId : HPrev(HOpposite(h));
}
inline H Trimesh::HPrevAroundEnd(H h) const { return HOpposite(HNext(h)); }
inline Vector3d Trimesh::HGeometry(H h) const {
  return VPosition(HEnd(h)) - VPosition(HStart(h));
}
std::generator<H> Trimesh::HConnectionsAroundStart(H st_h) const {
  H h = st_h;
  while (HOpposite(h) != kInvalidId) {
    assert(EdgeIsManifold(h));
    h = HPrevAroundStart(h);
    if (h == st_h) {
      break;
    }
  }
  st_h = h;
  do {
    assert(EdgeIsManifold(h));
    co_yield h;
    h = HNextAroundStart(h);
  } while (h != st_h && h != kInvalidId);
}
std::generator<H> Trimesh::HHalfedgesAroundHole(H st_h) const {
  assert(HOpposite(st_h) == kInvalidId);
  H h = st_h;
  do {
    co_yield h;
    h = HNext(h);
    while (HOpposite(h) != kInvalidId) {
      assert(EdgeIsManifold(h));
      h = HPrevAroundStart(h);
    }
  } while (h != st_h);
}
inline double Trimesh::HLength(H h) const { return HGeometry(h).norm(); }
inline H Trimesh::VStarting(V v) const { return vStart_.at(v); }
inline H Trimesh::VEnding(V v) const { return HPrev(VStarting(v)); }
std::generator<H> Trimesh::VStartings(V v) const {
  for (H h = VStarting(v); h != kInvalidId; h = hCoStart_[h]) {
    co_yield h;
  }
}
std::generator<H> Trimesh::VEndings(V v) const {
  for (H h : VStartings(v)) {
    co_yield HPrev(h);
  }
}
auto Trimesh::VFaces(TL::V v) const {
  return std::views::transform(VStartings(v), [this](H h) { return HFace(h); });
}
inline Vector3d& Trimesh::VPosition(V v) { return position_.at(v); }
inline const Vector3d& Trimesh::VPosition(V v) const { return position_.at(v); }
Vector3d Trimesh::VNormal(V v) const {
  Vector3d normal{Vector3d::Zero()};
  for (F f : VFaces(v)) {
    normal += FNormal(f);
  }
  return normal.normalized();
}
size_t Trimesh::VValence(V v) const {
  size_t ans = 0;
  for ([[maybe_unused]] H h : VStartings(v)) ans++;
  return ans;
}
bool Trimesh::VIsManifold(V v) const {
  return (std::ranges::to<std::vector>(HConnectionsAroundStart(VStarting(v))))
             .size() == (std::ranges::to<std::vector>(VStartings(v))).size();
}
inline H Trimesh::FHalfedge(F f) const {
  assert(f < NumFaces());
  return 3 * f;
}
inline std::ranges::iota_view<H, H> Trimesh::FHalfedges(F f) const {
  return std::views::iota(FHalfedge(f), H{3 * f + 3});
}
std::generator<F> Trimesh::FNeighbors(F f) const {
  for (H h : FHalfedges(f)) {
    H opp = HOpposite(h);
    if (opp != kInvalidId) {
      co_yield HFace(opp);
    }
  }
}
inline auto Trimesh::FVertices(F f) const {
  return std::views::transform(FHalfedges(f),
                               [this](H h) { return HStart(h); });
}
inline auto Trimesh::FPositions(F f) const {
  return std::views::transform(FVertices(f),
                               [this](V v) { return VPosition(v); });
}
inline Vector3d Trimesh::FNormal(F f) const {
  return HGeometry(3 * f).cross(-HGeometry(3 * f + 2)).normalized();
}
inline double Trimesh::FArea(F f) const {
  return 0.5 * HGeometry(3 * f).cross(HGeometry(3 * f + 2)).norm();
}
inline std::pair<Vector3d, Vector3d> Trimesh::FBoundingBox(F f) const {
  return {VPosition(HStart(3 * f))
              .cwiseMin(VPosition(HStart(3 * f + 1)))
              .cwiseMin(VPosition(HStart(3 * f + 2))),
          VPosition(HStart(3 * f))
              .cwiseMax(VPosition(HStart(3 * f + 1)))
              .cwiseMax(VPosition(HStart(3 * f + 2)))};
}
inline Vector3d Trimesh::FCentroid(F f) const {
  return (VPosition(HStart(3 * f)) + VPosition(HStart(3 * f + 1)) +
          VPosition(HStart(3 * f + 2))) /
         3.0;
}
std::generator<H> Trimesh::EdgeHalfedges(H h) const {
  for (H he : VStartings(HStart(h))) {
    if (HEnd(he) == HEnd(h)) {
      co_yield (he);
    }
  }
  for (H he : VStartings(HEnd(h))) {
    if (HEnd(he) == HStart(h)) {
      co_yield (he);
    }
  }
}
auto Trimesh::EdgeFaces(H h) const {
  return std::views::transform(EdgeHalfedges(h),
                               [this](H h) { return HFace(h); });
}
bool Trimesh::EdgeIsManifold(H h) const {
  for (H g : EdgeHalfedges(h)) {
    if (g != h && g != HOpposite(h)) {
      return false;
    }
  }
  return true;
}
inline auto Trimesh::BoundaryHalfedges() const {
  return std::views::filter(Halfedges(),
                            [this](H h) { return HOpposite(h) == kInvalidId; });
}
std::pair<Vector3d, Vector3d> Trimesh::BoundingBox() const {
  if (position_.empty()) {
    throw std::runtime_error("Mesh has no vertices.");
  }
  Vector3d min = position_.front();
  Vector3d max = min;
  for (const auto& pos : position_) {
    min = min.cwiseMin(pos);
    max = max.cwiseMax(pos);
  }
  return {min, max};
}
Vector3d Trimesh::Centroid() const {
  if (position_.empty()) {
    throw std::runtime_error("Mesh has no vertices.");
  }
  Vector3d centroid = Vector3d::Zero();
  for (const auto& pos : position_) {
    centroid += pos;
  }
  centroid /= static_cast<double>(position_.size());
  return centroid;
}
F Trimesh::AddFace(const std::array<std::variant<V, Vector3d>, 3>& triangle) {
  for (int i = 0; i < 3; i++) {
    if (std::holds_alternative<V>(triangle[i])) {
      for (int j = 0; j < i; j++) {
        if (triangle[i] == triangle[j]) {
          throw std::invalid_argument(
              "AddFace called with a degenerated triangle parameter");
        }
      }
    }
  }
  for (const auto& [key, f_attr] : attributes_[kFace]) {
    f_attr->IncrementSize();
  }
  H h_offset = hStart_.size();
  for (const auto& elem : triangle) {
    H h = hStart_.size();
    V v = position_.size();
    if (std::holds_alternative<V>(elem)) {
      v = std::get<V>(elem);
    } else {
      for (const auto& [key, v_attr] : attributes_[kVertex]) {
        v_attr->IncrementSize();
      }
      position_.push_back(std::get<Vector3d>(elem));
      vStart_.push_back(kInvalidId);
    }
    for (const auto& [key, h_attr] : attributes_[kHalfedge]) {
      h_attr->IncrementSize();
    }
    hStart_.push_back(v);
    hCoStart_.push_back(vStart_[v]);
    vStart_[v] = h;
  }
  return HFace(h_offset);
}
std::vector<V> Trimesh::RemoveFace(F f) {
  std::vector<V> removed_vertices;
  for (const auto& [key, f_attr] : attributes_[kFace]) {
    f_attr->ReplaceErase(f);
  }
  for (H h : std::views::reverse(FHalfedges(f))) {
    H g = hStart_.size() - 1;
    V v = hStart_[h];
    if (vStart_[v] == h) {
      vStart_[v] = hCoStart_[h];
    } else {
      for (H he : VStartings(v)) {
        if (hCoStart_[he] == h) {
          hCoStart_[he] = hCoStart_[h];
          break;
        }
      }
    }
    if (h != g) {
      V u = hStart_[g];
      hStart_[h] = u;
      hCoStart_[h] = hCoStart_[g];
      H* p = nullptr;
      if (vStart_[u] == g) {
        p = &vStart_[u];
      } else {
        for (H he : VStartings(u)) {
          if (hCoStart_[he] == g) {
            p = &hCoStart_[he];
            break;
          }
        }
      }
      *p = h;  // replace g with h and ensure the structure remains sorted
               // from greater to lower
      for (H nh = hCoStart_[h]; hCoStart_[h] > h && hCoStart_[h] != kInvalidId;
           nh = hCoStart_[h]) {
        *p = nh;
        hCoStart_[h] = hCoStart_[nh];
        hCoStart_[nh] = h;
        p = &hCoStart_[nh];
      }
    }
    hStart_.pop_back();
    hCoStart_.pop_back();
    for (const auto& [key, h_attr] : attributes_[kHalfedge]) {
      h_attr->ReplaceErase(h);
    }
    if (vStart_[v] == kInvalidId) {
      removed_vertices.push_back(v);
      for (H he : VStartings(vStart_.size() - 1)) {
        hStart_[he] = v;
      }
      vStart_[v] = vStart_.back(), vStart_.pop_back();
      position_[v] = std::move(position_.back()), position_.pop_back();
      for (const auto& [key, v_attr] : attributes_[kVertex]) {
        v_attr->ReplaceErase(v);
      }
    }
  }
  return removed_vertices;
}
std::vector<V> Trimesh::RemoveFaces(std::vector<F> faces) {
  std::ranges::sort(faces);
  const auto [first, last] = std::ranges::unique(faces);
  faces.erase(first, last);
  std::vector<V> removed_vertices;
  for (std::size_t i = 0; i < faces.size();) {
    F f = faces[i];
    std::vector<V> removed = RemoveFace(f);
    removed_vertices.insert(removed_vertices.end(), removed.begin(),
                            removed.end());
    if (faces.back() == NumFaces()) {
      faces.pop_back();
    } else {
      i++;
    }
  }
  return removed_vertices;
}

std::pair<std::vector<F>, std::vector<V>> Trimesh::CollapseEdge(H h) {
  std::vector<V> removed_vertices;
  std::vector<F> removed_faces = EdgeFaces(h) | std::ranges::to<std::vector>();
  std::array<V, 2> verts = {HStart(h), HEnd(h)};
  Vector3d midpoint = (VPosition(verts[0]) + VPosition(verts[1])) / 2.0;
  V last_vert_id = NumVertices() - 1;
  for (size_t i = 0; i < removed_faces.size(); i++) {
    std::vector<V> rm_verts = RemoveFace(removed_faces[i]);
    removed_vertices.insert(removed_vertices.end(), rm_verts.begin(),
                            rm_verts.end());
    for (size_t j = i + 1; j < removed_faces.size(); j++) {
      if (removed_faces[j] == NumFaces()) {
        removed_faces[j] = removed_faces[i];
      }
    }
  }
  for (const V& rem_v : removed_vertices) {
    for (V& v : verts) {
      if (v == rem_v) {
        v = kInvalidId;
      } else if (v == last_vert_id) {
        v = rem_v;
      }
    }
    last_vert_id--;
  }
  if (verts[0] == kInvalidId) {
    std::swap(verts[0], verts[1]);
  }
  if (verts[0] == kInvalidId) {
    return std::make_pair(std::move(removed_faces),
                          std::move(removed_vertices));
  }
  position_[verts[0]] = midpoint;
  std::vector<H> v_starts;
  for (V vert : verts) {
    if (vert != kInvalidId) {
      for (H he : VStartings(vert)) {
        v_starts.push_back(he);
      }
    }
  }
  std::sort(v_starts.begin(), v_starts.end());
  vStart_[verts[0]] = kInvalidId;
  for (H he : v_starts) {
    hStart_[he] = verts[0];
    hCoStart_[he] = vStart_[verts[0]];
    vStart_[verts[0]] = he;
  }
  if (verts[1] != kInvalidId) {
    removed_vertices.push_back(verts[1]);
    if (verts[1] != vStart_.size() - 1) {
      for (H he : VStartings(vStart_.size() - 1)) {
        hStart_[he] = verts[1];
      }
      vStart_[verts[1]] = vStart_.back();
      position_[verts[1]] = std::move(position_.back());
    }
    vStart_.pop_back();
    position_.pop_back();
    for (const auto& [key, v_attr] : attributes_[kVertex]) {
      v_attr->ReplaceErase(verts[1]);
    }
  }
  return std::make_pair(std::move(removed_faces), std::move(removed_vertices));
}
void Trimesh::DisconnectFace(F f) {
  AddFace({VPosition(HStart(3 * f)), VPosition(HStart(3 * f + 1)),
           VPosition(HStart(3 * f + 2))});
  RemoveFace(f);
}
void Trimesh::DisconnectFacesUntilManifoldEdges() {
  for (F f : std::views::reverse(Faces())) {
    for (H h : FHalfedges(f)) {
      if (!EdgeIsManifold(h)) {
        DisconnectFace(f);
        break;
      }
    }
  }
}
void Trimesh::DisconnectFacesUntilManifoldVertices() {
  std::queue<V> q;
  for (V v : Vertices()) {
    if (!VIsManifold(v)) {
      q.push(v);
    }
  }
  while (!q.empty()) {
    V v = q.front();
    q.pop();
    if (VIsManifold(v)) {
      continue;
    }
    std::unordered_set<V> s;
    for (H h : HConnectionsAroundStart(VStarting(v))) {
      s.insert(h);
    }
    for (H h : VStartings(v)) {
      if (!s.count(h)) {
        F f{HFace(h)};
        for (V u : FVertices(f)) {
          q.push(u);
        }
        DisconnectFace(f);
      }
    }
  }
}
void Trimesh::DisconnectFacesUntilManifold() {
  DisconnectFacesUntilManifoldEdges();
  DisconnectFacesUntilManifoldVertices();
}
template <MeshComponent C, typename T>
std::ranges::ref_view<std::vector<T>> Trimesh::CreateAttribute(
    const std::string& key) {
  if (attributes_[C].find(key) != attributes_[C].end()) {
    throw std::invalid_argument("Attribute already exists.");
  }
  size_t N = C == kHalfedge ? NumHalfedges()
             : C == kVertex ? NumVertices()
                            : NumFaces();
  attributes_[C][key] =
      std::make_unique<ContainerWrapper<std::vector<T>>>(std::vector<T>(N));
  return GetAttribute<C, T>(key);
}
template <MeshComponent C, typename T>
std::ranges::ref_view<std::vector<T>> Trimesh::GetAttribute(
    const std::string& key) {
  auto it = attributes_[C].find(key);
  if (it == attributes_[C].end()) {
    throw std::invalid_argument("Attribute does not exist.");
  }
  auto containerPtr =
      dynamic_cast<ContainerWrapper<std::vector<T>>*>(it->second.get());
  if (!containerPtr) {
    throw std::invalid_argument("Attribute exists but type does not match.");
  }
  return std::views::all(containerPtr->container_);
}
template <MeshComponent C>
void Trimesh::EraseAttribute(const std::string& key) {
  auto it = attributes_[C].find(key);
  if (it == attributes_[C].end()) {
    throw std::invalid_argument("Attribute does not exist.");
  }
  attributes_[C].erase(it);
}

}  // namespace TL

#endif  // TRILITE_GENERATOR_H
