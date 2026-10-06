#include "bvh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "closest_point.h"

namespace surface_field237::internal {

namespace {

constexpr std::uint32_t kLeafSize = 4;

// 点到轴对齐包围盒的平方距离（盒外最近点距离下界，盒内为零）。
double box_sq_distance(const Vec3& p, const Vec3& lo, const Vec3& hi) {
  double sq = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double d = p[k] < lo[k] ? lo[k] - p[k]
                     : p[k] > hi[k] ? p[k] - hi[k]
                                    : 0.0;
    sq += d * d;
  }
  return sq;
}

}  // namespace

TriBVH::TriBVH(const std::vector<Vec3>& pos, const std::vector<Tri>& tris)
    : pos_(pos), tris_(tris) {
  order_.resize(tris_.size());
  for (std::uint32_t i = 0; i < tris_.size(); ++i) order_[i] = i;
  if (!order_.empty()) build(0, static_cast<std::uint32_t>(order_.size()));
}

std::uint32_t TriBVH::build(std::uint32_t begin, std::uint32_t end) {
  Node node;
  node.lo = Vec3::Constant(std::numeric_limits<double>::infinity());
  node.hi = Vec3::Constant(-std::numeric_limits<double>::infinity());
  Vec3 cen_lo = node.lo;
  Vec3 cen_hi = node.hi;
  for (std::uint32_t i = begin; i < end; ++i) {
    const Tri& t = tris_[order_[i]];
    Vec3 cen = Vec3::Zero();
    for (int k = 0; k < 3; ++k) {
      node.lo = node.lo.cwiseMin(pos_[t[k]]);
      node.hi = node.hi.cwiseMax(pos_[t[k]]);
      cen += pos_[t[k]];
    }
    cen /= 3.0;
    cen_lo = cen_lo.cwiseMin(cen);
    cen_hi = cen_hi.cwiseMax(cen);
  }

  const std::uint32_t index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back(node);

  if (end - begin <= kLeafSize) {
    nodes_[index].leaf = true;
    nodes_[index].begin = begin;
    nodes_[index].end = end;
    return index;
  }

  // 按质心最长轴中位切分。
  int axis = 0;
  double span = cen_hi.x() - cen_lo.x();
  if (cen_hi.y() - cen_lo.y() > span) {
    axis = 1;
    span = cen_hi.y() - cen_lo.y();
  }
  if (cen_hi.z() - cen_lo.z() > span) axis = 2;
  const std::uint32_t mid = begin + (end - begin) / 2;
  std::nth_element(order_.begin() + begin, order_.begin() + mid,
                   order_.begin() + end, [&](std::uint32_t x, std::uint32_t y) {
                     auto centroid = [&](std::uint32_t f) {
                       const Tri& t = tris_[f];
                       return (pos_[t[0]][axis] + pos_[t[1]][axis] +
                               pos_[t[2]][axis]) /
                              3.0;
                     };
                     return centroid(x) < centroid(y);
                   });

  const std::uint32_t left = build(begin, mid);
  const std::uint32_t right = build(mid, end);
  nodes_[index].left = static_cast<std::int32_t>(left);
  nodes_[index].right = static_cast<std::int32_t>(right);
  return index;
}

std::uint32_t TriBVH::nearest_face(const Vec3& p, double max_dist,
                                   double* out_dist) const {
  if (nodes_.empty()) return UINT32_MAX;

  double best_sq = max_dist * max_dist;  // 距离等于上限可接收
  std::uint32_t best_face = UINT32_MAX;

  // 手工栈遍历；仅当包围盒下界严格大于当前最优才剪枝，
  // 等距情形继续深入以保证面 ID 较小者优先。
  std::vector<std::uint32_t> stack;
  stack.push_back(0);
  while (!stack.empty()) {
    const Node& node = nodes_[stack.back()];
    stack.pop_back();
    if (box_sq_distance(p, node.lo, node.hi) > best_sq) continue;

    if (!node.leaf) {
      stack.push_back(static_cast<std::uint32_t>(node.left));
      stack.push_back(static_cast<std::uint32_t>(node.right));
      continue;
    }
    for (std::uint32_t i = node.begin; i < node.end; ++i) {
      const std::uint32_t fid = order_[i];
      const Tri& t = tris_[fid];
      const Vec3 q = closest_point_on_triangle(p, pos_[t[0]], pos_[t[1]],
                                               pos_[t[2]]);
      const double sq = (p - q).squaredNorm();
      if (sq < best_sq || (sq == best_sq && fid < best_face)) {
        best_sq = sq;
        best_face = fid;
      }
    }
  }

  if (best_face == UINT32_MAX) return UINT32_MAX;
  if (out_dist) *out_dist = std::sqrt(best_sq);
  return best_face;
}

}  // namespace surface_field237::internal
