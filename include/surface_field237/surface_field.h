#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "surface_field237/mesh_lod.h"

namespace surface_field237 {

// 单点投射结果。projected 为 false 时其余字段不含有效值
// （数值字段为 NaN，face_id 为 UINT32_MAX），不以零值冒充。
struct ProjectionResult {
  static constexpr double kInvalid =
      std::numeric_limits<double>::quiet_NaN();

  bool projected = false;
  Vec3 closest_point = Vec3::Constant(kInvalid);
  std::uint32_t face_id = UINT32_MAX;  // 原网格面 ID
  // 按该面顶点索引顺序的重心坐标。
  std::array<double, 3> barycentric{kInvalid, kInvalid, kInvalid};
  double distance = kInvalid;               // 欧氏距离
  double scalar = kInvalid;                 // 重心插值标量
};

// 由原扫描表面（顶点 + 有向三角索引 + 逐顶点标量）建立的可复用查询对象。
// 对象持有输入的独立快照，构造后修改输入数组不影响查询；
// 内部三角形 AABB 层次索引供多批查询复用。
class SurfaceField {
 public:
  // 网格合法性规则与 simplify 一致（允许开口与多分量）；
  // 标量须逐顶点、长度一致且全部有限。非法输入抛 std::invalid_argument。
  SurfaceField(std::vector<Vec3> vertices, std::vector<Tri> triangles,
               std::vector<double> scalars);
  ~SurfaceField();
  SurfaceField(SurfaceField&&) noexcept;
  SurfaceField& operator=(SurfaceField&&) noexcept;
  SurfaceField(const SurfaceField&) = delete;
  SurfaceField& operator=(const SurfaceField&) = delete;

  // 批量投射同坐标系查询点。max_distance 须非负有限；距离等于上限可接收，
  // 超限点标记 projected=false。等距时选原面 ID 较小者。
  // 任一点坐标非有限或上限非法时整批拒绝（抛 std::invalid_argument），
  // 对象状态不受影响。输出按输入顺序一一对应。
  std::vector<ProjectionResult> project(const std::vector<Vec3>& points,
                                        double max_distance) const;

  const std::vector<Vec3>& vertices() const;
  const std::vector<Tri>& triangles() const;
  const std::vector<double>& scalars() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// 携带投射结果的减面输出。projections 与 mesh.vertices 一一对应；
// 未投射顶点保留在网格中并以 projected=false 明确标记。
struct SimplifyFieldResult {
  SimplifyResult mesh;
  std::vector<ProjectionResult> projections;
};

// 减面后以最终输出顶点查询原表面标量场。输入按值语义处理，不被修改；
// 投射复用同一份空间索引与重心插值。
SimplifyFieldResult simplify_with_field(
    const std::vector<Vec3>& vertices, const std::vector<Tri>& triangles,
    const std::vector<double>& scalars, std::size_t target_face_count,
    double max_projection_distance);

}  // namespace surface_field237
