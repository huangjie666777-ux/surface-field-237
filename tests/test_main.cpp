#include <cmath>
#include <cstdint>
#include <map>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "surface_field237/mesh_lod.h"
#include "surface_field237/surface_field.h"

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

  // 边长 1e-7 的有效三角形不得被面积下限误拒
  {
    const double s = 1e-7;
    std::vector<Vec3> v{{0, 0, 0}, {s, 0, 0}, {0, s, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    bool ok = true;
    try {
      const auto r = surface_field237::simplify(v, t, 1);
      ok = r.actual_face_count == 1;
    } catch (const std::invalid_argument&) {
      ok = false;
    }
    expect(ok, "tiny 1e-7 triangle accepted");
  }
  // 减面后不得保留未被引用的顶点
  {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    make_grid(10, v, t);
    const auto r = surface_field237::simplify(v, t, 30);
    std::vector<bool> used(r.vertices.size(), false);
    for (const Tri& f : r.triangles)
      for (auto x : f) used[x] = true;
    bool all_used = true;
    for (bool u : used)
      if (!u) all_used = false;
    expect(all_used, "no unreferenced vertices after simplify");
  }

  // ---- 表面标量场投射 ----
  // 单三角形：面内/边/顶点最近点与重心插值
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    std::vector<double> s{0.0, 10.0, 20.0};
    const surface_field237::SurfaceFieldQuery q(v, t, s);
    const std::vector<Vec3> pts{{0.25, 0.25, 1.0},   // 面内上方
                                {0.5, -0.5, 0.0},    // 边 AB 外侧
                                {2.0, 2.0, 0.0}};    // 顶点 B/C 之外
    const auto r = q.project(pts, 10.0);
    bool ok = r.size() == 3 && r[0].projected && r[1].projected &&
              r[2].projected;
    // 面内：最近点 (0.25,0.25,0)，重心 (0.5,0.25,0.25)，标量 7.5
    ok = ok && std::abs(r[0].closest_point.z()) < 1e-12 &&
         std::abs(r[0].distance - 1.0) < 1e-12 &&
         std::abs(r[0].barycentric[0] - 0.5) < 1e-12 &&
         std::abs(r[0].scalar - 7.5) < 1e-9 && r[0].face_id == 0;
    // 边：最近点 (0.5,0,0)，标量 5
    ok = ok && (r[1].closest_point - Vec3(0.5, 0, 0)).norm() < 1e-12 &&
         std::abs(r[1].scalar - 5.0) < 1e-9;
    // 顶点区域之外：最近点为边 BC 上某点或顶点，距离有限
    ok = ok && r[2].distance > 0.0;
    expect(ok, "exact closest point: face/edge/vertex regions");
  }
  // 等距取较小面 ID；边界距离等于上限可接受
  {
    // 两个共顶点三角形，查询点在公共边中垂线上等距。
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, -1, 0}};
    std::vector<Tri> t{{0, 1, 2}, {1, 0, 3}};
    std::vector<double> s{0, 1, 2, 3};
    const surface_field237::SurfaceFieldQuery q(v, t, s);
    const auto r = q.project({Vec3(0.5, 0.0, 1.0)}, 1.0);
    expect(r[0].projected && r[0].face_id == 0 &&
               std::abs(r[0].distance - 1.0) < 1e-12,
           "tie broken by smaller face id, boundary distance accepted");
  }
  // 超限标记未投射，标量非零冒充；非法批次整批拒绝且后续查询正常
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    std::vector<double> s{1, 2, 3};
    const surface_field237::SurfaceFieldQuery q(v, t, s);
    const auto r = q.project({Vec3(0.2, 0.2, 5.0)}, 1.0);
    bool ok = !r[0].projected && std::isnan(r[0].scalar);
    // 非法批次：非有限点 / 负距离
    bool threw = false;
    try {
      q.project({Vec3(0, 0, std::nan(""))}, 1.0);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    bool threw2 = false;
    try {
      q.project({Vec3(0, 0, 0)}, -1.0);
    } catch (const std::invalid_argument&) {
      threw2 = true;
    }
    // 后续查询不受影响
    const auto r2 = q.project({Vec3(0.2, 0.2, 0.5)}, 1.0);
    ok = ok && threw && threw2 && r2[0].projected &&
         std::abs(r2[0].distance - 0.5) < 1e-12;
    expect(ok, "out-of-range flagged, invalid batch rejected, query reusable");
  }
  // 快照独立：构造后修改输入不影响查询
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    std::vector<double> s{1, 2, 3};
    const surface_field237::SurfaceFieldQuery q(v, t, s);
    v[0] = Vec3(100, 100, 100);
    s[0] = 999.0;
    t[0] = Tri{2, 1, 0};
    const auto r = q.project({Vec3(0.0, 0.0, 1.0)}, 2.0);
    expect(r[0].projected && std::abs(r[0].scalar - 1.0) < 1e-12,
           "query object holds independent snapshot");
  }
  // 标量长度不匹配 / 非有限标量拒绝
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    bool threw = false;
    try {
      surface_field237::SurfaceFieldQuery q(v, t, {1.0, 2.0});
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    bool threw2 = false;
    try {
      surface_field237::SurfaceFieldQuery q(v, t, {1.0, 2.0, std::nan("")});
    } catch (const std::invalid_argument&) {
      threw2 = true;
    }
    expect(threw && threw2, "reject scalar length mismatch / non-finite");
  }
  // simplify_with_field：减面顶点投射回原表面，未投射顶点保留并标记
  {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    make_grid(10, v, t);
    std::vector<double> s(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
      s[i] = v[i].x() + 2.0 * v[i].y();
    const auto r =
        surface_field237::simplify_with_field(v, t, s, 50, 5.0);
    bool ok = r.projections.size() == r.mesh.vertices.size() &&
              r.mesh.actual_face_count <= 50;
    for (const auto& p : r.projections)
      ok = ok && p.projected && p.distance <= 5.0;
    // 输入未被修改
    ok = ok && v.size() == 100;
    // 极小的 max_distance：部分顶点可能未投射，须明确标记且仍保留
    const auto r2 = surface_field237::simplify_with_field(v, t, s, 50, 0.0);
    ok = ok && r2.projections.size() == r2.mesh.vertices.size();
    bool any_flag = false;
    for (const auto& p : r2.projections)
      if (!p.projected) any_flag = true;
    ok = ok && any_flag;
    expect(ok, "simplify_with_field projects output vertices, flags misses");
  }

  std::cout << (g_failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return g_failures == 0 ? 0 : 1;
}
