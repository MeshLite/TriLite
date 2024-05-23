
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

#ifndef TRILITE_SAMPLING_H
#define TRILITE_SAMPLING_H

#include "Trimesh.hpp"

namespace TL {

class Sampling {
 public:
  static std::generator<std::pair<Vector3d, double>> BarycentersAndAreas(
      const TL::Trimesh& mesh, double max_edge_length) {
    double sq_length = max_edge_length * max_edge_length;
    std::function<std::generator<std::pair<Vector3d, double>>(
        std::array<Vector3d, 3>, double)>
        process_triangle =
            [&](std::array<Vector3d, 3>&& tri,
                double st_area) -> std::generator<std::pair<Vector3d, double>> {
      if (std::max({(tri[1] - tri[0]).squaredNorm(),
                    (tri[2] - tri[1]).squaredNorm(),
                    (tri[0] - tri[2]).squaredNorm()}) <= sq_length) {
        co_yield std::make_pair((tri[0] + tri[1] + tri[2]) / 3.0, st_area);
        co_return;
      }
      st_area /= 4.0;
      Vector3d a = (tri[0] + tri[1]) / 2.0;
      Vector3d b = (tri[1] + tri[2]) / 2.0;
      Vector3d c = (tri[2] + tri[0]) / 2.0;
      for (std::array<Vector3d, 3> sub_tri :
           std::array<std::array<Vector3d, 3>, 4>{
               std::array<Vector3d, 3>{a, b, c},
               std::array<Vector3d, 3>{tri[0], a, c},
               std::array<Vector3d, 3>{tri[1], b, a},
               std::array<Vector3d, 3>{tri[2], c, b}}) {
        for (std::pair<Vector3d, double> point_area :
             process_triangle(sub_tri, st_area)) {
          co_yield point_area;
        }
      }
    };

    for (TL::F f : mesh.Faces()) {
      std::array<Vector3d, 3> tri;
      for (auto [i, position] : std::views::enumerate(mesh.FPositions(f))) {
        tri[i] = position;
      }
      for (std::pair<Vector3d, double> point_area :
           process_triangle(std::move(tri), mesh.FArea(f))) {
        co_yield point_area;
      }
    }
  }
};
}  // namespace TL
#endif