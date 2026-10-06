#pragma once

#include <cstdint>
#include <vector>

#include "surface_field237/mesh_lod.h"

namespace surface_field237::internal {

// 三角形 AABB 二叉层次索引。构建一次，供多批最近面查询复用；
// 以包围盒距离下界剪枝，叶内做精确点到三角形距离计算。
class TriBVH {
 public:
  // pos/tris 由调用方持有，生命周期须覆盖本对象。
  TriBVH(const std::vector<Vec3>& pos, const std::vector<Tri>& tris);

  // 返回距离 p 最近（且不超过 max_dist，含等距边界）的原面 ID；
  // 等距时取面 ID 较小者。未找到返回 UINT32_MAX。
  // out_dist 非空时写回最近距离。
  std::uint32_t nearest_face(const Vec3& p, double max_dist,
                             double* out_dist) const;

 private:
  struct Node {
    Vec3 lo;
    Vec3 hi;
    std::int32_t left = -1;   // 内部节点子节点下标
    std::int32_t right = -1;
    std::uint32_t begin = 0;  // 叶子：tris 置换区间 [begin, end)
    std::uint32_t end = 0;
    bool leaf = false;
  };

  const std::vector<Vec3>& pos_;
  const std::vector<Tri>& tris_;
  std::vector<Node> nodes_;
  std::vector<std::uint32_t> order_;  // 叶子内三角形原面 ID 置换

  std::uint32_t build(std::uint32_t begin, std::uint32_t end);
};

}  // namespace surface_field237::internal
