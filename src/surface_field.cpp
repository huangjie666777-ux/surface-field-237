#include "surface_field237/surface_field.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#include "validation.h"

namespace surface_field237 {

namespace {

[[noreturn]] void fail(const std::string& msg) {
  throw std::invalid_argument("surface_field237: " + msg);
}

struct Aabb {
  Vec3 lo = Vec3::Constant(std::numeric_limits<double>::infinity());
  Vec3 hi = Vec3::Constant(-std::numeric_limits<double>::infinity());

  void extend(const Vec3& p) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  void extend(const Aabb& o) {
    lo = lo.cwiseMin(o.lo);
    hi = hi.cwiseMax(o.hi);
  }
  // 点到包围盒的欧氏距离下界。
  double distance_lower_bound(const Vec3& p) const {
    double s = 0.0;
    for (int k = 0; k < 3; ++k) {
      const double d =
          p[k] < lo[k] ? lo[k] - p[k] : (p[k] > hi[k] ? p[k] - hi[k] : 0.0);
      s += d * d;
    }
    return std::sqrt(s);
  }
};

// 精确点到三角形最近点（Christer Ericson, RTCD 5.1.5），
// 覆盖面内、边与顶点区域；返回最近点，重心坐标由调用方另行计算。
Vec3 closest_point_triangle(const Vec3& p, const Vec3& a, const Vec3& b,
                            const Vec3& c) {
  const Vec3 ab = b - a;
  const Vec3 ac = c - a;
  const Vec3 ap = p - a;
  const double d1 = ab.dot(ap);
  const double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0) return a;  // 顶点 A 区域

  const Vec3 bp = p - b;
  const double d3 = ab.dot(bp);
  const double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3) return b;  // 顶点 B 区域

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    return a + (d1 / (d1 - d3)) * ab;  // 边 AB

  const Vec3 cp = p - c;
  const double d5 = ab.dot(cp);
  const double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6) return c;  // 顶点 C 区域

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    return a + (d2 / (d2 - d6)) * ac;  // 边 AC

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);  // 边 BC

  // 面内区域。
  const double denom = 1.0 / (va + vb + vc);
  return a + ab * (vb * denom) + ac * (vc * denom);
}

// 由最近点（位于三角形平面内）求重心坐标，顺序与面索引一致。
std::array<double, 3> barycentric_of(const Vec3& p, const Vec3& a,
                                     const Vec3& b, const Vec3& c) {
  const Vec3 v0 = b - a;
  const Vec3 v1 = c - a;
  const Vec3 v2 = p - a;
  const double d00 = v0.dot(v0);
  const double d01 = v0.dot(v1);
  const double d11 = v1.dot(v1);
  const double d20 = v2.dot(v0);
  const double d21 = v2.dot(v1);
  const double denom = d00 * d11 - d01 * d01;  // 非退化面保证为正
  const double v = (d11 * d20 - d01 * d21) / denom;
  const double w = (d00 * d21 - d01 * d20) / denom;
  return {1.0 - v - w, v, w};
}

}  // namespace

struct SurfaceFieldQuery::Impl {
  std::vector<Vec3> vertices;
  std::vector<Tri> triangles;
  std::vector<double> scalars;

  struct Node {
    Aabb box;
    std::int32_t left = -1;   // 内部节点子节点下标
    std::int32_t right = -1;
    std::uint32_t begin = 0;  // 叶子：tri_order 区间 [begin, end)
    std::uint32_t end = 0;
    bool leaf = false;
  };
  std::vector<Node> nodes;
  std::vector<std::uint32_t> tri_order;  // 叶子引用的原面 ID 排列

  static constexpr std::size_t kLeafSize = 8;

  Impl(std::vector<Vec3> v, std::vector<Tri> t, std::vector<double> s)
      : vertices(std::move(v)),
        triangles(std::move(t)),
        scalars(std::move(s)) {
    build_bvh();
  }

  void build_bvh() {
    const std::size_t n = triangles.size();
    std::vector<Aabb> boxes(n);
    std::vector<Vec3> centroids(n);
    for (std::size_t i = 0; i < n; ++i) {
      const Tri& t = triangles[i];
      Aabb b;
      b.extend(vertices[t[0]]);
      b.extend(vertices[t[1]]);
      b.extend(vertices[t[2]]);
      boxes[i] = b;
      centroids[i] = (vertices[t[0]] + vertices[t[1]] + vertices[t[2]]) / 3.0;
    }
    tri_order.resize(n);
    std::iota(tri_order.begin(), tri_order.end(), 0u);
    nodes.reserve(2 * n);
    build_node(0, static_cast<std::uint32_t>(n), boxes, centroids);
  }

  std::uint32_t build_node(std::uint32_t begin, std::uint32_t end,
                           const std::vector<Aabb>& boxes,
                           const std::vector<Vec3>& centroids) {
    const std::uint32_t idx = static_cast<std::uint32_t>(nodes.size());
    nodes.emplace_back();
    Aabb box;
    for (std::uint32_t i = begin; i < end; ++i) box.extend(boxes[tri_order[i]]);
    nodes[idx].box = box;
    if (end - begin <= kLeafSize) {
      nodes[idx].leaf = true;
      nodes[idx].begin = begin;
      nodes[idx].end = end;
      return idx;
    }
    // 按包围盒质心最长轴中位切分。
    const Vec3 ext = box.hi - box.lo;
    int axis = 0;
    if (ext.y() > ext.x() && ext.y() >= ext.z()) axis = 1;
    else if (ext.z() > ext.x() && ext.z() >= ext.y()) axis = 2;
    const std::uint32_t mid = begin + (end - begin) / 2;
    std::nth_element(tri_order.begin() + begin, tri_order.begin() + mid,
                     tri_order.begin() + end,
                     [&](std::uint32_t x, std::uint32_t y) {
                       return centroids[x][axis] < centroids[y][axis];
                     });
    nodes[idx].left = build_node(begin, mid, boxes, centroids);
    nodes[idx].right = build_node(mid, end, boxes, centroids);
    return idx;
  }

  ProjectionResult query_one(const Vec3& p, double max_distance) const {
    ProjectionResult best;
    double best_dist = max_distance;  // 等界可接受：仅当下界更大才剪枝
    std::uint32_t best_face = std::numeric_limits<std::uint32_t>::max();
    bool found = false;

    // 栈式遍历，以包围盒距离下界剪枝；下界不超过当前最优才可能更优
    // （或等距但面 ID 更小），不得跳过。
    std::vector<std::uint32_t> stack;
    stack.push_back(0);
    while (!stack.empty()) {
      const Node& node = nodes[stack.back()];
      stack.pop_back();
      if (node.box.distance_lower_bound(p) > best_dist) continue;
      if (!node.leaf) {
        stack.push_back(static_cast<std::uint32_t>(node.left));
        stack.push_back(static_cast<std::uint32_t>(node.right));
        continue;
      }
      for (std::uint32_t i = node.begin; i < node.end; ++i) {
        const std::uint32_t fid = tri_order[i];
        const Tri& t = triangles[fid];
        const Vec3 q = closest_point_triangle(p, vertices[t[0]],
                                              vertices[t[1]], vertices[t[2]]);
        const double d = (p - q).norm();
        if (!found || d < best_dist || (d == best_dist && fid < best_face)) {
          found = true;
          best_dist = d;
          best_face = fid;
          best.projected = true;
          best.closest_point = q;
          best.face_id = fid;
          best.barycentric = barycentric_of(q, vertices[t[0]], vertices[t[1]],
                                            vertices[t[2]]);
          best.distance = d;
          best.scalar = best.barycentric[0] * scalars[t[0]] +
                        best.barycentric[1] * scalars[t[1]] +
                        best.barycentric[2] * scalars[t[2]];
        }
      }
    }
    if (!found) {
      best.projected = false;
      best.scalar = std::numeric_limits<double>::quiet_NaN();
    }
    return best;
  }
};

SurfaceFieldQuery::SurfaceFieldQuery(const std::vector<Vec3>& vertices,
                                     const std::vector<Tri>& triangles,
                                     const std::vector<double>& scalars) {
  if (scalars.size() != vertices.size())
    fail("scalar count does not match vertex count");
  for (std::size_t i = 0; i < scalars.size(); ++i)
    if (!std::isfinite(scalars[i]))
      fail("non-finite scalar at vertex " + std::to_string(i));
  // 网格合法性沿用减面校验规则（允许开口与多分量）。
  internal::build_validated_state(vertices, triangles);
  // 独立快照：构造后修改输入不影响查询。
  impl_ = std::make_shared<Impl>(vertices, triangles, scalars);
}

std::vector<ProjectionResult> SurfaceFieldQuery::project(
    const std::vector<Vec3>& points, double max_distance) const {
  if (!std::isfinite(max_distance) || max_distance < 0.0)
    fail("max projection distance must be non-negative and finite");
  for (std::size_t i = 0; i < points.size(); ++i) {
    const Vec3& p = points[i];
    if (!std::isfinite(p.x()) || !std::isfinite(p.y()) ||
        !std::isfinite(p.z()))
      fail("non-finite query point at index " + std::to_string(i));
  }
  std::vector<ProjectionResult> out;
  out.reserve(points.size());
  for (const Vec3& p : points) out.push_back(impl_->query_one(p, max_distance));
  return out;
}

std::size_t SurfaceFieldQuery::vertex_count() const {
  return impl_->vertices.size();
}

std::size_t SurfaceFieldQuery::triangle_count() const {
  return impl_->triangles.size();
}

SimplifyFieldResult simplify_with_field(const std::vector<Vec3>& vertices,
                                        const std::vector<Tri>& triangles,
                                        const std::vector<double>& scalars,
                                        std::size_t target_face_count,
                                        double max_distance) {
  // 共用同一校验、空间索引与插值实现；输入按值语义不被修改。
  const SurfaceFieldQuery query(vertices, triangles, scalars);
  SimplifyFieldResult result;
  result.mesh = simplify(vertices, triangles, target_face_count);
  result.projections = query.project(result.mesh.vertices, max_distance);
  return result;
}

}  // namespace surface_field237
