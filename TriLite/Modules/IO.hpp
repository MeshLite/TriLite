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

#ifndef TRILITE_IO_H
#define TRILITE_IO_H

#include <fstream>
#include <iomanip>

#include "../Core/Trimesh.hpp"
namespace TL {

class IO {
 public:
  /**
   * @brief Determines the file format and calls the corresponding mesh reader
   * function.
   * @param filepath The path to the mesh file.
   * @return A Trimesh object representing the mesh.
   * @throws std::runtime_error if the file cannot be read or has an unsupported
   * format.
   */
  static Trimesh ReadMeshFile(const std::string& filepath);

  /**
   * @brief Writes a triangle mesh to a file, determining the appropriate format
   * from the file extension.
   * @param mesh The Trimesh object representing the mesh.
   * @param filepath The path to the output file.
   * @param binary_mode Whether the file must be written in binary mode (if
   * compatible with file extension).
   * @throws std::runtime_error if the file format is unsupported or if the file
   * cannot be opened for writing.
   */
  static void WriteMeshFile(const Trimesh& mesh, const std::string& filepath,
                            bool binary_mode = true);

 private:
  // Internal structures and functions
  static const int kAsciiDigitsPrecision_ = 10;
  struct Vector3Hash {
    std::size_t operator()(const Vector3d& key) const {
      const std::hash<double> hasher;
      size_t result = 0;
      for (int i = 0; i < 3; ++i) {
        result ^= hasher(key[i]) + 0x9e3779b9 + (result << 6) + (result >> 2);
      }
      return result;
    }
  };
  struct Float3Hash {
    std::size_t operator()(const std::array<float, 3>& key) const {
      const std::hash<float> hasher;
      size_t result = 0;
      for (int i = 0; i < 3; ++i) {
        result ^= hasher(key[i]) + 0x9e3779b9 + (result << 6) + (result >> 2);
      }
      return result;
    }
  };
  static std::vector<std::variant<V, Vector3d>> ConvertedToTrimeshPoints(
      const std::vector<Vector3d>& points) {
    std::vector<std::variant<V, Vector3d>> triangle_points;
    std::unordered_map<Vector3d, V, Vector3Hash> position_to_vertex;
    int n_duplicates = 0;
    for (size_t i = 0; i < points.size(); i += 3) {
      for (int j = 0; j < 3; j++) {
        auto [it, inserted] = position_to_vertex.emplace(
            points[i + j],
            static_cast<V>(position_to_vertex.size() + n_duplicates));
        for (int k = 0; k < j; k++) {
          if (points[i + j] == points[i + k]) {
            inserted = true;
            n_duplicates++;
            break;
          }
        }
        if (inserted) {
          triangle_points.push_back(points[i + j]);
        } else {
          triangle_points.push_back(it->second);
        }
      }
    }
    return triangle_points;
  }
  static std::vector<unsigned char> ReadBinaryFile(
      const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<unsigned char> buffer(size);
    file.read(reinterpret_cast<char*>(buffer.data()), size);
    return buffer;
  }
  static Trimesh ReadAsciiSTL(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line, token;
    std::vector<Vector3d> positions;
    Vector3d normal, vertex;
    while (std::getline(file, line)) {
      std::istringstream linestream(line);
      linestream >> token;
      if (token == "facet") {
        linestream >> token;
        linestream >> normal[0] >> normal[1] >> normal[2];
      } else if (token == "vertex") {
        linestream >> vertex[0] >> vertex[1] >> vertex[2];
        positions.push_back(vertex);
      } else if (token == "endfacet") {
        // Remove vertices if they don't count to three (error in file)
        positions.resize(positions.size() - (positions.size() % 3));
      }
    }
    positions.resize(positions.size() - (positions.size() % 3));
    return Trimesh(ConvertedToTrimeshPoints(positions));
  }
  static Trimesh ReadBinarySTL(const std::string& filepath) {
    const std::vector<unsigned char> binaryData = ReadBinaryFile(filepath);
    int nTriangles = *reinterpret_cast<const int*>(binaryData.data() + 80);
    std::vector<Vector3d> points(3 * nTriangles);
    float x;
    for (int i = 0; i < nTriangles; ++i) {
      auto ptr = binaryData.data() + 80 + 4 + 50 * i;
      for (int j : std::views::iota(1, 4)) {
        for (int k : std::views::iota(0, 3)) {
          memcpy(&x, ptr + j * 12 + k * 4, 4);
          points[3 * i + j - 1][k] = x;
        }
      }
    }
    return Trimesh(ConvertedToTrimeshPoints(points));
  }
  static Trimesh ReadSTL(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    char header[80];
    int nTriangles;
    file.read(header, 80);
    file.read(reinterpret_cast<char*>(&nTriangles), sizeof(nTriangles));
    file.seekg(0, std::ios::end);
    std::streampos fileSize = file.tellg();
    if (nTriangles * 50 + 84 == fileSize || strncmp(header, "solid", 5)) {
      return ReadBinarySTL(filepath);
    } else {
      return ReadAsciiSTL(filepath);
    }
  }
  static Trimesh ReadOBJ(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;
    std::vector<Vector3d> points;
    std::vector<std::variant<V, Vector3d>> trimesh_points;
    std::vector<V> id_to_vid;
    Vector3d p;
    std::string s_index;
    V count_vids = 0;
    while (std::getline(file, line)) {
      std::istringstream lineStream(line);
      std::string lineType;
      lineStream >> lineType;
      if (lineType == "v") {
        lineStream >> p.x() >> p.y() >> p.z();
        points.push_back(std::move(p));
        id_to_vid.push_back({kInvalidId});
      } else if (lineType == "f") {
        std::vector<int> indices;
        while (lineStream >> s_index) {
          s_index = s_index.substr(0, s_index.find('/'));
          indices.push_back(std::stoi(s_index) - 1);
        }
        for (size_t i = 1; i < indices.size() - 1; ++i) {
          int idx[] = {indices[0], indices[i], indices[i + 1]};
          if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2]) {
            continue;
          }
          for (int j : idx) {
            if (id_to_vid[j] != kInvalidId) {
              trimesh_points.push_back(id_to_vid[j]);
            } else {
              trimesh_points.push_back(points[j]);
              id_to_vid[j] = {count_vids++};
            }
          }
        }
      }
    }

    return Trimesh(trimesh_points);
  }
  static Trimesh ReadOFF(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;
    std::getline(file, line);
    if (line != "OFF") {
      throw std::runtime_error("Invalid OFF file format");
    }
    int numVertices, numFaces, numEdges;
    file >> numVertices >> numFaces >> numEdges;
    std::vector<Vector3d> points(numVertices);
    std::vector<std::variant<V, Vector3d>> trimesh_points;
    std::vector<V> id_to_vid(numVertices, {kInvalidId});
    V count_vids = 0;
    for (Vector3d& p : points) {
      file >> p.x() >> p.y() >> p.z();
    }
    int numVerticesPerFace;
    for (int f = 0; f < numFaces; ++f) {
      file >> numVerticesPerFace;
      std::vector<int> indices(numVerticesPerFace);
      for (int& index : indices) {
        file >> index;
      }
      for (size_t i = 1; i < indices.size() - 1; ++i) {
        int idx[] = {indices[0], indices[i], indices[i + 1]};
        if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2]) {
          continue;
        }
        for (int j : idx) {
          if (id_to_vid[j] != kInvalidId) {
            trimesh_points.push_back(id_to_vid[j]);
          } else {
            trimesh_points.push_back(points[j]);
            id_to_vid[j] = {count_vids++};
          }
        }
      }
    }
    return Trimesh(trimesh_points);
  }
  static Trimesh ReadBinaryPLY(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    std::string line;
    bool header_ended = false;
    int point_count = 0;
    int face_count = 0;
    bool is_big_endian = false;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;
      if (line.find("binary_big_endian") != std::string::npos) {
        is_big_endian = true;
      } else if (token == "end_header") {
        header_ended = true;
        break;
      } else if (token == "element") {
        iss >> token;
        if (token == "vertex") {
          iss >> point_count;
        } else if (token == "face") {
          iss >> face_count;
        }
      }
    }
    if (!header_ended) {
      throw std::runtime_error("PLY end_header is missing.");
    }
    auto convert_endian = [](float& value) {
      char* value_ptr = reinterpret_cast<char*>(&value);
      std::reverse(value_ptr, value_ptr + sizeof(float));
    };
    std::vector<Vector3d> points(point_count);
    for (int i = 0; i < point_count; ++i) {
      std::array<float, 3> point;
      for (float& xyz : point) {
        file.read(reinterpret_cast<char*>(&xyz), sizeof(float));
        if (is_big_endian) {
          convert_endian(xyz);
        }
      }
      points[i] = Vector3d{point[0], point[1], point[2]};
    }
    V count_vids = 0;
    std::vector<V> id_to_vid(point_count, kInvalidId);
    std::vector<std::variant<V, Vector3d>> trimesh_points;
    for (int f = 0; f < face_count; ++f) {
      unsigned char vertex_count;
      file.read(reinterpret_cast<char*>(&vertex_count), sizeof(vertex_count));
      std::vector<int> indices(vertex_count);
      for (int j = 0; j < vertex_count; ++j) {
        int vertex_index;
        file.read(reinterpret_cast<char*>(&vertex_index), sizeof(vertex_index));
        indices[j] = vertex_index;
      }
      for (size_t i = 1; i < indices.size() - 1; ++i) {
        int idx[] = {indices[0], indices[i], indices[i + 1]};
        if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2]) {
          continue;
        }
        for (int j : idx) {
          if (id_to_vid[j] != kInvalidId) {
            trimesh_points.push_back(id_to_vid[j]);
          } else {
            trimesh_points.push_back(points[j]);
            id_to_vid[j] = count_vids++;
          }
        }
      }
    }
    return Trimesh(trimesh_points);
  }
  static Trimesh ReadAsciiPLY(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;
    bool header_ended = false;
    int point_count = 0;
    int face_count = 0;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;
      if (token == "end_header") {
        header_ended = true;
        break;
      } else if (token == "element") {
        iss >> token;
        if (token == "vertex") {
          iss >> point_count;
        } else if (token == "face") {
          iss >> face_count;
        }
      }
    }
    if (!header_ended) {
      throw std::runtime_error("PLY end_header is missing.");
    }
    std::vector<Vector3d> points(point_count);
    for (Vector3d& p : points) {
      file >> p.x() >> p.y() >> p.z();
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::vector<std::variant<V, Vector3d>> trimesh_points;
    std::vector<V> id_to_vid(point_count, kInvalidId);
    V count_vids = 0;
    for (int i = 0; i < face_count; ++i) {
      int vertex_count;
      file >> vertex_count;
      std::vector<int> indices(vertex_count);
      for (int& index : indices) {
        file >> index;
      }
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      for (size_t j = 1; j + 1 < indices.size(); ++j) {
        int idx[] = {indices[0], indices[j], indices[j + 1]};
        if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2]) {
          continue;
        }
        for (int index : idx) {
          if (id_to_vid[index] != kInvalidId) {
            trimesh_points.push_back(id_to_vid[index]);
          } else {
            trimesh_points.push_back(points[index]);
            id_to_vid[index] = count_vids++;
          }
        }
      }
    }
    return Trimesh(trimesh_points);
  }
  static Trimesh ReadPLY(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    std::string line;
    getline(file, line);
    getline(file, line);
    file.close();
    if (line.find("ascii") != std::string::npos) {
      return ReadAsciiPLY(filepath);
    } else if (line.find("binary_little_endian") != std::string::npos) {
      return ReadBinaryPLY(filepath);
    } else if (line.find("binary_big_endian") != std::string::npos) {
      return ReadBinaryPLY(filepath);
    } else {
      throw std::runtime_error("Unsupported PLY format in file: " + filepath);
    }
  }
  static void WriteOBJ(const std::string& filepath, const Trimesh& mesh) {
    std::ofstream file(filepath);
    file << std::setprecision(kAsciiDigitsPrecision_);
    for (const auto& position : mesh.Positions()) {
      file << "v " << position.transpose() << '\n';
    }
    for (F f : mesh.Faces()) {
      file << "f";
      for (V v : mesh.FVertices(f)) {
        file << " " << v + 1;
      }
      file << '\n';
    }
    file.close();
  }
  static void WriteOFF(const std::string& filepath, const Trimesh& mesh) {
    std::ofstream file(filepath);
    file << std::setprecision(kAsciiDigitsPrecision_);
    file << "OFF\n" << mesh.NumVertices() << " " << mesh.NumFaces() << " 0\n";
    for (const auto& position : mesh.Positions()) {
      file << position.transpose() << '\n';
    }
    for (auto face : mesh.Faces()) {
      file << "3";
      for (V v : mesh.FVertices(face)) {
        file << " " << v;
      }
      file << '\n';
    }
    file.close();
  }
  static void WriteSTL(const std::string& filepath, const Trimesh& mesh,
                       bool binary_mode) {
    if (binary_mode) {
      std::unordered_map<std::array<float, 3>, V, Float3Hash> position_to_id;
      position_to_id.reserve(mesh.NumVertices());
      std::ofstream file(filepath, std::ios::binary);
      char header[80] = {};
      file.write(header, sizeof(header));
      uint32_t num_triangles = static_cast<uint32_t>(mesh.NumFaces());
      file.write(reinterpret_cast<const char*>(&num_triangles),
                 sizeof(num_triangles));
      auto write_float3 = [&file, &position_to_id](const Vector3d& vec,
                                                   V v = kInvalidId) {
        std::array<float, 3> data{static_cast<float>(vec[0]),
                                  static_cast<float>(vec[1]),
                                  static_cast<float>(vec[2])};
        if (v != kInvalidId) {
          // Displace the point until it is unique
          for (size_t i = 0;
               position_to_id.find(data) != position_to_id.end() &&
               position_to_id[data] != v;
               ++i) {
            if (!std::isfinite(data[i % 3])) {
              data[i % 3] = -std::numeric_limits<float>::infinity();
            }
            data[i % 3] = std::nextafter(
                data[i % 3], std::numeric_limits<float>::infinity());
          }
          position_to_id[data] = v;
        }
        file.write(reinterpret_cast<const char*>(&data), sizeof(data));
      };
      for (F f : mesh.Faces()) {
        write_float3(mesh.FNormal(f));
        for (V v : mesh.FVertices(f)) {
          write_float3(mesh.VPosition(v), v);
        }
        uint16_t attribute_byte_count = 0;
        file.write(reinterpret_cast<const char*>(&attribute_byte_count),
                   sizeof(attribute_byte_count));
      }
      file.close();
    } else {
      double increment = pow(10, 1 - kAsciiDigitsPrecision_);
      std::vector<Vector3d> positions{
          std::ranges::to<std::vector>(mesh.Positions())};
      std::vector<V> ids{std::ranges::to<std::vector>(mesh.Vertices())};
      std::sort(ids.begin(), ids.end(), [&positions](int u, int v) {
        return positions[u][0] < positions[v][0];
      });

      for (size_t i = 1; i < ids.size(); i++) {
        V u = ids[i - 1], v = ids[i];
        double ref =
            std::max(std::nextafter(positions[u][0],
                                    std::numeric_limits<double>::infinity()),
                     positions[u][0] + abs(positions[u][0]) * increment);
        positions[v][0] = std::max(positions[v][0], ref);
      }
      std::ofstream file(filepath);
      file << std::setprecision(kAsciiDigitsPrecision_) << "solid mesh\n";
      for (F f : mesh.Faces()) {
        file << "  facet normal " << mesh.FNormal(f).transpose() << "\n";
        file << "    outer loop\n";
        for (V v : mesh.FVertices(f)) {
          file << "    vertex " << positions[v].transpose() << "\n";
        }
        file << "    endloop\n";
        file << "  endfacet\n";
      }
      file << "endsolid mesh\n";
      file.close();
    }
  }
  static void WritePLY(const std::string& filepath, const Trimesh& mesh,
                       bool binary_mode) {
    std::ofstream file = binary_mode ? std::ofstream(filepath, std::ios::binary)
                                     : std::ofstream(filepath);
    if (!binary_mode) {
      file << std::setprecision(kAsciiDigitsPrecision_);
    }
    file << "ply\n";
    file << (binary_mode ? "format binary_little_endian 1.0\n"
                         : "format ascii 1.0\n");
    file << "element vertex " << mesh.NumVertices() << "\n";
    file << "property float x\n"
         << "property float y\n"
         << "property float z\n";
    file << "element face " << mesh.NumFaces() << "\n";
    file << "property list uchar int vertex_index\n";
    file << "end_header\n";
    for (const Vector3d& point : mesh.Positions()) {
      if (binary_mode) {
        for (float coord : point) {
          file.write(reinterpret_cast<const char*>(&coord), sizeof(float));
        }
      } else {
        file << point.transpose() << "\n";
      }
    }
    for (F f : mesh.Faces()) {
      unsigned char num_vertices = 3;
      if (binary_mode) {
        file.write(reinterpret_cast<const char*>(&num_vertices),
                   sizeof(num_vertices));
        for (V v : mesh.FVertices(f)) {
          file.write(reinterpret_cast<const char*>(&v), sizeof(v));
        }
      } else {
        file << 3 << " " << mesh.HStart(3 * f) << " " << mesh.HStart(3 * f + 1)
             << " " << mesh.HStart(3 * f + 2) << "\n";
      }
    }
    file.close();
  }
};

Trimesh IO::ReadMeshFile(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file for reading: " + filepath);
  }
  file.close();
  size_t dotPos = filepath.rfind('.');
  if (dotPos == std::string::npos) {
    throw std::runtime_error("No file extension found");
  }
  std::string extension = filepath.substr(dotPos + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (extension == "stl") {
    return ReadSTL(filepath);
  } else if (extension == "obj") {
    return ReadOBJ(filepath);
  } else if (extension == "off") {
    return ReadOFF(filepath);
  } else if (extension == "ply") {
    return ReadPLY(filepath);
  } else {
    throw std::runtime_error("Unsupported file format: " + extension);
  }
}
void IO::WriteMeshFile(const Trimesh& mesh, const std::string& filepath,
                       bool binary_mode) {
  size_t dotPos = filepath.rfind('.');
  if (dotPos == std::string::npos) {
    throw std::runtime_error("No file extension found");
  }
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file for writing: " + filepath);
  }
  file.close();
  std::string extension = filepath.substr(dotPos + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (extension == "stl") {
    WriteSTL(filepath, mesh, binary_mode);
  } else if (extension == "obj") {
    WriteOBJ(filepath, mesh);
  } else if (extension == "off") {
    WriteOFF(filepath, mesh);
  } else if (extension == "ply") {
    WritePLY(filepath, mesh, binary_mode);
  } else {
    throw std::runtime_error("Unsupported file format: " + extension);
  }
}

}  // namespace TL

#endif  // TRILITE_IO_H
