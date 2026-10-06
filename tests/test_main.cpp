#include <cmath>
#include <cstdint>
#include <map>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "surface_field237/mesh_lod.h"

using surface_field237::StopReason;
using surface_field237::Tri;
using surface_field237::Vec3;

namespace {

int g_failures = 0;

void expect(bool cond, const std::string& name) {
  if (cond) {
    std::cout << "[PASS] " << name << "\n";
  } else {
    std::cout << "[FAIL] " << name << "\n";
    ++g_failures;
  }
}

bool throws_invalid(const std::vector<Vec3>& v, const std::vector<Tri>& t,
                    std::size_t target) {
  try {
    surface_field237::simplify(v, t, target);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

std::vector<Vec3> tetra_verts() {
  return {{0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 0.5, 1}};
}
std::vector<Tri> tetra_tris() {
  return {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
}

void make_grid(int n, std::vector<Vec3>& verts, std::vector<Tri>& tris) {
  for (int y = 0; y < n; ++y)
    for (int x = 0; x < n; ++x)
      verts.emplace_back(x, y, 0.2 * std::sin(0.4 * x) * std::cos(0.4 * y));
  auto id = [n](int x, int y) { return static_cast<std::uint32_t>(y * n + x); };
  for (int y = 0; y + 1 < n; ++y)
    for (int x = 0; x + 1 < n; ++x) {
      tris.push_back(Tri{id(x, y), id(x + 1, y), id(x + 1, y + 1)});
      tris.push_back(Tri{id(x, y), id(x + 1, y + 1), id(x, y + 1)});
    }
}

}  // namespace

int main() {
  // 非有限坐标
  {
    auto v = tetra_verts();
    v[1].x() = std::nan("");
    expect(throws_invalid(v, tetra_tris(), 2), "reject non-finite coordinate");
  }
  // 越界索引
  {
    auto t = tetra_tris();
    t[0][0] = 99;
    expect(throws_invalid(tetra_verts(), t, 2), "reject out-of-range index");
  }
  // 重复面
  {
    auto v = tetra_verts();
    std::vector<Tri> t{{0, 1, 2}, {2, 1, 0}};
    expect(throws_invalid(v, t, 1), "reject duplicate face");
  }
  // 零面积
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    expect(throws_invalid(v, t, 1), "reject zero-area face");
  }
  // 绕向不一致（共享边同向）
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}, {1, 2, 3}};  // 边 (1,2) 同向
    expect(throws_invalid(v, t, 1), "reject inconsistent winding");
  }
  // 非流形边（3 面共边）
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, -1, 0}};
    std::vector<Tri> t{{0, 1, 2}, {1, 0, 3}, {0, 1, 4}};
    expect(throws_invalid(v, t, 1), "reject non-manifold edge");
  }
  // 非流形顶点（两个四面体共一个顶点）
  {
    std::vector<Vec3> v = tetra_verts();
    v.push_back({0.5, 0.5, 2});
    v.push_back({0.5, 1.5, 1});
    v.push_back({1.5, 0.5, 1});
    std::vector<Tri> t = tetra_tris();
    t.push_back({3, 4, 5});
    t.push_back({3, 5, 6});
    t.push_back({3, 6, 4});
    t.push_back({4, 6, 5});
    expect(throws_invalid(v, t, 2), "reject non-manifold vertex");
  }
  // 目标为 0
  {
    bool threw = false;
    try {
      surface_field237::simplify(tetra_verts(), tetra_tris(), 0);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "reject zero target");
  }
  // 目标不小于当前面数：原样返回
  {
    const auto v = tetra_verts();
    const auto t = tetra_tris();
    const auto r = surface_field237::simplify(v, t, 4);
    expect(r.stop_reason == StopReason::kAlreadyAtOrBelowTarget &&
               r.vertices == v && r.triangles == t && r.collapse_count == 0,
           "target >= current returns input unchanged");
  }
  // 封闭四面体：任何折叠都会产生重复面，应保持不变
  {
    const auto r = surface_field237::simplify(tetra_verts(), tetra_tris(), 2);
    expect(r.stop_reason == StopReason::kNoValidCandidate &&
               r.actual_face_count == 4 && r.collapse_count == 0,
           "tetrahedron has no legal collapse");
  }
  // 单个三角形：边界全锁定，无候选
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    const auto r = surface_field237::simplify(v, t, 1);
    expect(r.stop_reason == StopReason::kAlreadyAtOrBelowTarget ||
               (r.stop_reason == StopReason::kNoValidCandidate &&
                r.actual_face_count == 1),
           "single triangle boundary locked");
  }
  // 开口网格：减面成功、边界坐标保持、绕向保持
  {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    make_grid(10, v, t);
    const std::size_t before = t.size();
    const auto r = surface_field237::simplify(v, t, 50);
    expect(r.actual_face_count <= 50 && r.actual_face_count < before &&
               r.collapse_count > 0,
           "grid simplified to target");
    // 输出合法性：索引在界内、无重复面、无零面积
    bool valid = true;
    std::set<std::array<std::uint32_t, 3>> keys;
    for (const Tri& f : r.triangles) {
      for (auto x : f)
        if (x >= r.vertices.size()) valid = false;
      auto k = f;
      std::sort(k.begin(), k.end());
      if (!keys.insert(k).second) valid = false;
      const Vec3& a = r.vertices[f[0]];
      const Vec3& b = r.vertices[f[1]];
      const Vec3& c = r.vertices[f[2]];
      if ((b - a).cross(c - a).norm() <= 1e-12) valid = false;
    }
    expect(valid, "output faces valid (range/dup/area)");
    // 边界顶点坐标保持不变
    std::set<std::tuple<double, double, double>> in_coords, out_coords;
    auto boundary = [](const std::vector<Tri>& tris) {
      std::map<std::pair<std::uint32_t, std::uint32_t>, int> ec;
      for (const Tri& f : tris)
        for (int k = 0; k < 3; ++k) {
          std::uint32_t a = f[k], b = f[(k + 1) % 3];
          if (a > b) std::swap(a, b);
          ++ec[{a, b}];
        }
      std::set<std::uint32_t> out;
      for (const auto& [e, c] : ec)
        if (c == 1) {
          out.insert(e.first);
          out.insert(e.second);
        }
      return out;
    };
    for (auto idx : boundary(t))
      in_coords.emplace(v[idx].x(), v[idx].y(), v[idx].z());
    for (auto idx : boundary(r.triangles))
      out_coords.emplace(r.vertices[idx].x(), r.vertices[idx].y(),
                         r.vertices[idx].z());
    expect(in_coords == out_coords, "boundary coordinates preserved");
  }
  // 多分量：两个独立开口网格，互不连接
  {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    make_grid(6, v, t);
    const std::size_t off = v.size();
    std::vector<Vec3> v2;
    std::vector<Tri> t2;
    make_grid(6, v2, t2);
    for (auto& p : v2) p.x() += 100.0;
    v.insert(v.end(), v2.begin(), v2.end());
    for (const Tri& f : t2)
      t.push_back(Tri{static_cast<std::uint32_t>(f[0] + off),
                      static_cast<std::uint32_t>(f[1] + off),
                      static_cast<std::uint32_t>(f[2] + off)});
    const auto r = surface_field237::simplify(v, t, 40);
    bool cross = false;
    for (const Tri& f : r.triangles) {
      const bool a = r.vertices[f[0]].x() > 50;
      const bool b = r.vertices[f[1]].x() > 50;
      const bool c = r.vertices[f[2]].x() > 50;
      if (a != b || b != c) cross = true;
    }
    expect(!cross && r.actual_face_count <= 40,
           "components not connected, target reached");
  }

  std::cout << (g_failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return g_failures == 0 ? 0 : 1;
}
