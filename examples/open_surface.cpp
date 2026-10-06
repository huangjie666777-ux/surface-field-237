#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "surface_field237/mesh_lod.h"

using surface_field237::Tri;
using surface_field237::Vec3;

namespace {

// 开口曲面：NxN 顶点的波浪网格（带自由边界）。
void make_wave_grid(int n, std::vector<Vec3>& verts, std::vector<Tri>& tris) {
  for (int y = 0; y < n; ++y)
    for (int x = 0; x < n; ++x)
      verts.emplace_back(static_cast<double>(x), static_cast<double>(y),
                         0.3 * std::sin(0.5 * x) * std::cos(0.5 * y));
  auto id = [n](int x, int y) { return static_cast<std::uint32_t>(y * n + x); };
  for (int y = 0; y + 1 < n; ++y)
    for (int x = 0; x + 1 < n; ++x) {
      tris.push_back(Tri{id(x, y), id(x + 1, y), id(x + 1, y + 1)});
      tris.push_back(Tri{id(x, y), id(x + 1, y + 1), id(x, y + 1)});
    }
}

std::set<std::uint32_t> boundary_vertices(const std::vector<Tri>& tris) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_count;
  for (const Tri& t : tris)
    for (int k = 0; k < 3; ++k) {
      std::uint32_t a = t[k], b = t[(k + 1) % 3];
      if (a > b) std::swap(a, b);
      ++edge_count[{a, b}];
    }
  std::set<std::uint32_t> out;
  for (const auto& [e, c] : edge_count)
    if (c == 1) {
      out.insert(e.first);
      out.insert(e.second);
    }
  return out;
}

}  // namespace

int main() {
  std::vector<Vec3> verts;
  std::vector<Tri> tris;
  make_wave_grid(12, verts, tris);

  const std::size_t target = 60;
  const auto before_boundary = boundary_vertices(tris);

  const surface_field237::SimplifyResult r = surface_field237::simplify(verts, tris, target);

  std::cout << "input : " << verts.size() << " vertices, " << tris.size()
            << " triangles\n";
  std::cout << "output: " << r.vertices.size() << " vertices, "
            << r.actual_face_count << " triangles (target " << target << ")\n";
  std::cout << "collapses: " << r.collapse_count << "\n";
  std::cout << "stop reason: " << surface_field237::to_string(r.stop_reason) << "\n";

  // 边界保持验证：输出边界顶点坐标须为输入边界顶点坐标的子集。
  const auto after_boundary = boundary_vertices(r.triangles);
  std::set<std::tuple<double, double, double>> input_coords;
  for (std::uint32_t idx : before_boundary)
    input_coords.emplace(verts[idx].x(), verts[idx].y(), verts[idx].z());
  bool boundary_preserved = true;
  for (std::uint32_t idx : after_boundary) {
    const Vec3& p = r.vertices[idx];
    if (!input_coords.count({p.x(), p.y(), p.z()})) boundary_preserved = false;
  }
  std::cout << "boundary vertices: " << before_boundary.size() << " -> "
            << after_boundary.size()
            << ", coords preserved: " << (boundary_preserved ? "yes" : "NO")
            << "\n";
  return boundary_preserved ? 0 : 1;
}

