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

bool nearly(double a, double b, double tol = 1e-9) {
  return std::fabs(a - b) <= tol;
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
  // 边长 1e-7 的有效三角形：相对容差下不得被误拒
  {
    const double s = 1e-7;
    std::vector<Vec3> v{{0, 0, 0}, {s, 0, 0}, {0.5 * s, 0.8 * s, 0}};
    std::vector<Tri> t{{0, 1, 2}};
    bool accepted = true;
    try {
      surface_field237::simplify(v, t, 1);
    } catch (const std::invalid_argument&) {
      accepted = false;
    }
    expect(accepted, "tiny 1e-7 triangle accepted");
  }
  // 减面输出不保留未引用顶点
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
  // SurfaceField：构造校验
  {
    bool threw = false;
    try {
      surface_field237::SurfaceField f(tetra_verts(), tetra_tris(),
                                       {1.0, 2.0});  // 长度不匹配
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "field rejects scalar length mismatch");
    threw = false;
    try {
      surface_field237::SurfaceField f(tetra_verts(), tetra_tris(),
                                       {1.0, std::nan(""), 2.0, 3.0});
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "field rejects non-finite scalar");
  }
  // SurfaceField：面内/边/顶点投射、插值、快照独立、整批拒绝
  {
    std::vector<Vec3> v{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    std::vector<Tri> t{{0, 1, 2}, {2, 1, 3}};
    std::vector<double> s{0.0, 1.0, 2.0, 3.0};
    surface_field237::SurfaceField f(v, t, s);
    // 快照独立：修改输入不影响查询
    v[0] = {99, 99, 99};
    s[0] = 1e9;
    t[0] = {3, 2, 1};

    const auto hits = f.project(
        {{0.25, 0.25, 1.0},   // 面内（面 0），重心 (0.5,0.25,0.25)
         {0.5, 0.0, -0.5},    // 边 (0,1) 中点下方
         {0.0, 0.0, 2.0},     // 顶点 0 正上方
         {0.0, 0.0, 3.0}},    // 超出上限
        2.0);
    bool ok = hits.size() == 4;
    ok = ok && hits[0].projected && hits[0].face_id == 0 &&
         nearly(hits[0].distance, 1.0) &&
         nearly(hits[0].closest_point.z(), 0.0) &&
         nearly(hits[0].barycentric[0], 0.5) &&
         nearly(hits[0].barycentric[1], 0.25) &&
         nearly(hits[0].barycentric[2], 0.25) &&
         nearly(hits[0].scalar, 0.5 * 0.0 + 0.25 * 1.0 + 0.25 * 2.0);
    ok = ok && hits[1].projected && hits[1].face_id == 0 &&
         nearly(hits[1].closest_point.x(), 0.5) &&
         nearly(hits[1].closest_point.y(), 0.0) &&
         nearly(hits[1].scalar, 0.5);  // 边上两端标量均值
    ok = ok && hits[2].projected && nearly(hits[2].distance, 2.0) &&
         nearly(hits[2].scalar, 0.0);  // 距离等于上限可接收
    ok = ok && !hits[3].projected && std::isnan(hits[3].scalar) &&
         hits[3].face_id == UINT32_MAX;
    expect(ok, "projection interior/edge/vertex/limit/beyond");

    // 等距取面 ID 较小者：公共边中点
    const auto tie = f.project({{0.5, 0.5, 0.0}}, 1.0);
    expect(tie[0].projected && tie[0].face_id == 0,
           "tie broken by smaller face id");

    // 非法查询整批拒绝，且后续查询不受影响
    bool threw = false;
    try {
      f.project({{0.0, 0.0, std::nan("")}}, 1.0);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "reject batch with non-finite point");
    threw = false;
    try {
      f.project({{0.0, 0.0, 0.0}}, -1.0);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "reject negative max distance");
    const auto after = f.project({{0.25, 0.25, 0.0}}, 1.0);
    expect(after[0].projected && nearly(after[0].scalar, 0.75),
           "queries still valid after rejected batch");

    // 索引复用：多批查询结果一致
    const auto again = f.project({{0.25, 0.25, 1.0}}, 2.0);
    expect(again[0].projected && again[0].face_id == hits[0].face_id &&
               nearly(again[0].scalar, hits[0].scalar),
           "index reused across batches");
  }
  // simplify_with_field：投射与减面联动
  {
    std::vector<Vec3> v;
    std::vector<Tri> t;
    make_grid(10, v, t);
    std::vector<double> s(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
      s[i] = v[i].x() + 2.0 * v[i].y();  // 非恒定场
    const auto r = surface_field237::simplify_with_field(v, t, s, 40, 1.0);
    bool ok = r.projections.size() == r.mesh.vertices.size();
    std::size_t projected = 0;
    for (std::size_t i = 0; i < r.projections.size(); ++i) {
      const auto& h = r.projections[i];
      if (!h.projected) {
        ok = ok && std::isnan(h.scalar);  // 未投射保留并明确标记
        continue;
      }
      ++projected;
      ok = ok && h.distance <= 1.0 && h.face_id < t.size();
      ok = ok && nearly(h.barycentric[0] + h.barycentric[1] +
                            h.barycentric[2],
                        1.0, 1e-9);
    }
    expect(ok && projected > 0, "simplify_with_field projections consistent");
    // 输入未被修改
    expect(v.size() == 100 && t.size() == 162, "inputs unchanged");
  }

  std::cout << (g_failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return g_failures == 0 ? 0 : 1;
}
