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

  // 非恒定标量场：温度状分布。
  std::vector<double> field(verts.size());
  for (std::size_t i = 0; i < verts.size(); ++i)
    field[i] = 20.0 + 5.0 * std::sin(0.3 * verts[i].x()) *
                          std::cos(0.7 * verts[i].y());

  // 减面并将输出顶点投射回原表面，附带插值标量。
  const double max_dist = 2.0;
  const auto r = surface_field237::simplify_with_field(verts, tris, field, 60,
                                                       max_dist);
  std::size_t projected = 0;
  for (const auto& p : r.projections)
    if (p.projected) ++projected;
  std::cout << "simplified: " << r.mesh.vertices.size() << " vertices, "
            << r.mesh.actual_face_count << " triangles ("
            << surface_field237::to_string(r.mesh.stop_reason) << ")\n";
  std::cout << "projected " << projected << " / " << r.projections.size()
            << " output vertices onto original surface\n";
  for (std::size_t i = 0; i < 3 && i < r.projections.size(); ++i) {
    const auto& p = r.projections[i];
    std::cout << "  v" << i << ": face " << p.face_id << " dist "
              << p.distance << " scalar " << p.scalar << "\n";
  }

  // 独立批量查询：含一个超出最大投射距离的点。
  const surface_field237::SurfaceFieldQuery query(verts, tris, field);
  const std::vector<Vec3> points{{3.0, 4.0, 0.0},      // 曲面附近
                                 {100.0, 100.0, 0.0}};  // 远超 max_dist
  const auto res = query.project(points, 1.0);
  for (std::size_t i = 0; i < res.size(); ++i) {
    std::cout << "query " << i << ": "
              << (res[i].projected ? "projected" : "NOT projected");
    if (res[i].projected)
      std::cout << " face " << res[i].face_id << " dist " << res[i].distance
                << " scalar " << res[i].scalar;
    std::cout << "\n";
  }
  return (res[0].projected && !res[1].projected) ? 0 : 1;
}
