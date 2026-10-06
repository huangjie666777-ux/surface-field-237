#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "surface_field237/surface_field.h"

using surface_field237::Tri;
using surface_field237::Vec3;

namespace {

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

}  // namespace

int main() {
  std::vector<Vec3> verts;
  std::vector<Tri> tris;
  make_wave_grid(12, verts, tris);

  // 非恒定标量场：模拟原扫描表面的测量值。
  std::vector<double> field(verts.size());
  for (std::size_t i = 0; i < verts.size(); ++i)
    field[i] = verts[i].x() + 2.0 * verts[i].y() +
               0.5 * std::sin(verts[i].z());

  // 直接批量投射：包含曲面附近点、上方点与超限点。
  surface_field237::SurfaceField query(verts, tris, field);
  const std::vector<Vec3> points{
      {3.2, 4.7, 0.0},     // 曲面附近
      {8.9, 1.1, 2.0},     // 上方约 2 个单位
      {100.0, 100.0, 0.0}  // 远超上限
  };
  const auto hits = query.project(points, 1.0);
  std::cout << "direct projection (max distance 1.0):\n";
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (!h.projected) {
      std::cout << "  point " << i << ": not projected (beyond limit)\n";
      continue;
    }
    std::cout << "  point " << i << ": face " << h.face_id << " distance "
              << h.distance << " scalar " << h.scalar << "\n";
  }

  // 减面并携带投射结果。
  const auto r =
      surface_field237::simplify_with_field(verts, tris, field, 60, 0.5);
  std::size_t projected = 0;
  for (const auto& h : r.projections)
    if (h.projected) ++projected;
  std::cout << "simplify_with_field: " << r.mesh.vertices.size()
            << " vertices, " << r.mesh.actual_face_count << " triangles, "
            << projected << " projected, " << r.projections.size() - projected
            << " unprojected (kept, flagged)\n";
  std::cout << "stop reason: "
            << surface_field237::to_string(r.mesh.stop_reason) << "\n";
  return 0;
}
