#include "mesh_state.h"

#include <cmath>

namespace surface_field237::internal {

namespace {
constexpr double kDegenerateEps = 1e-12;
}

double tri_double_area(const Vec3& a, const Vec3& b, const Vec3& c) {
  return (b - a).cross(c - a).norm();
}

bool tri_degenerate(const Vec3& a, const Vec3& b, const Vec3& c) {
  const double l0 = (b - a).norm();
  const double l1 = (c - b).norm();
  const double l2 = (a - c).norm();
  const double longest = std::max({l0, l1, l2});
  // 相对容差：按最长边平方缩放，避免误拒边长极小（如 1e-7）的有效三角形。
  const double tol = kDegenerateEps * longest * longest;
  return tri_double_area(a, b, c) <= tol;
}

void MeshState::add_face(std::uint32_t fid) {
  Face& f = faces[fid];
  f.alive = true;
  face_keys.insert(face_key(f.v));
  for (int i = 0; i < 3; ++i) {
    const std::uint32_t a = f.v[i];
    const std::uint32_t b = f.v[(i + 1) % 3];
    v_faces[a].insert(fid);
    EdgeRecord& rec = edges[edge_key(a, b)];
    rec.faces[rec.count++] = fid;
  }
  ++alive_faces;
}

void MeshState::remove_face(std::uint32_t fid) {
  Face& f = faces[fid];
  face_keys.erase(face_key(f.v));
  for (int i = 0; i < 3; ++i) {
    const std::uint32_t a = f.v[i];
    const std::uint32_t b = f.v[(i + 1) % 3];
    v_faces[a].erase(fid);
    auto it = edges.find(edge_key(a, b));
    EdgeRecord& rec = it->second;
    int slot = (rec.faces[0] == fid) ? 0 : 1;
    rec.faces[slot] = rec.faces[rec.count - 1];
    --rec.count;
    if (rec.count == 0) edges.erase(it);
  }
  f.alive = false;
  --alive_faces;
}

std::set<std::uint32_t> MeshState::neighbors(std::uint32_t x) const {
  std::set<std::uint32_t> out;
  for (std::uint32_t fid : v_faces[x]) {
    const Face& f = faces[fid];
    for (int i = 0; i < 3; ++i)
      if (f.v[i] != x) out.insert(f.v[i]);
  }
  return out;
}

std::vector<std::uint32_t> MeshState::opposite_vertices(std::uint32_t a,
                                                        std::uint32_t b) const {
  std::vector<std::uint32_t> out;
  auto it = edges.find(edge_key(a, b));
  if (it == edges.end()) return out;
  const EdgeRecord& rec = it->second;
  for (int i = 0; i < rec.count; ++i) {
    const Face& f = faces[rec.faces[i]];
    for (int k = 0; k < 3; ++k)
      if (f.v[k] != a && f.v[k] != b) out.push_back(f.v[k]);
  }
  return out;
}

}  // namespace surface_field237::internal
