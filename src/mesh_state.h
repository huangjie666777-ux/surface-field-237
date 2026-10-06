#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "surface_field237/mesh_lod.h"

namespace surface_field237::internal {

using EdgeKey = std::pair<std::uint32_t, std::uint32_t>;
using FaceKey = std::array<std::uint32_t, 3>;

inline EdgeKey edge_key(std::uint32_t a, std::uint32_t b) {
  return a < b ? EdgeKey{a, b} : EdgeKey{b, a};
}

inline FaceKey face_key(const Tri& t) {
  FaceKey k{t[0], t[1], t[2]};
  if (k[0] > k[1]) std::swap(k[0], k[1]);
  if (k[1] > k[2]) std::swap(k[1], k[2]);
  if (k[0] > k[1]) std::swap(k[0], k[1]);
  return k;
}

struct EdgeRecord {
  std::uint32_t faces[2] = {UINT32_MAX, UINT32_MAX};
  int count = 0;
};

struct MeshState {
  std::vector<Vec3> pos;
  std::vector<Eigen::Matrix4d> quad;
  std::vector<bool> locked;
  std::vector<bool> v_alive;

  struct Face {
    Tri v;
    bool alive = true;
  };
  std::vector<Face> faces;

  std::vector<std::set<std::uint32_t>> v_faces;
  std::map<EdgeKey, EdgeRecord> edges;
  std::set<FaceKey> face_keys;
  std::size_t alive_faces = 0;

  void add_face(std::uint32_t fid);
  void remove_face(std::uint32_t fid);
  std::set<std::uint32_t> neighbors(std::uint32_t x) const;
  std::vector<std::uint32_t> opposite_vertices(std::uint32_t a,
                                               std::uint32_t b) const;
};

double tri_double_area(const Vec3& a, const Vec3& b, const Vec3& c);
bool tri_degenerate(const Vec3& a, const Vec3& b, const Vec3& c);

}  // namespace surface_field237::internal

