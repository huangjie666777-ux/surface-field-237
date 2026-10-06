#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "surface_field237/mesh_lod.h"

namespace surface_field237 {

// 单点投射结果。projected 为 false 时仅 projected 字段有效，
// scalar 为 NaN（不以零值冒充），其余字段不参与解释。
struct ProjectionResult {
  bool projected = false;
  Vec3 closest_point = Vec3::Zero();
  std::uint32_t face_id = 0;                    // 原网格面 ID
  std::array<double, 3> barycentric{0, 0, 0};   // 按该面索引顺序
  double distance = 0.0;
  double scalar = 0.0;
};

// 由源顶点、有向三角索引与逐顶点标量建立的可复用查询对象。
// 构造时校验并持有独立快照；内部维护三角形 AABB 层次索引，
// 供任意多批查询复用。线程安全：所有查询方法均为 const。
class SurfaceFieldQuery {
 public:
  // 非法输入（非有限坐标/标量、长度不匹配、网格不合法）抛
  // std::invalid_argument，网格合法性规则与 simplify 一致。
  SurfaceFieldQuery(const std::vector<Vec3>& vertices,
                    const std::vector<Tri>& triangles,
                    const std::vector<double>& scalars);

  // 批量投射：points 与网格同坐标系，max_distance 为非负有限值。
  // 任一点坐标非有限或 max_distance 非法时整批拒绝（抛
  // std::invalid_argument），对象状态不受影响。
  // 返回与输入同序；距离不超过 max_distance（含等界）才视为投射成功，
  // 等距候选取原面 ID 较小者。
  std::vector<ProjectionResult> project(const std::vector<Vec3>& points,
                                        double max_distance) const;

  std::size_t vertex_count() const;
  std::size_t triangle_count() const;

 private:
  struct Impl;
  std::shared_ptr<const Impl> impl_;
};

// 减面并携带投射结果：以原网格与标量场建立查询对象，
// 将减面输出的每个顶点投射回原表面。
struct SimplifyFieldResult {
  SimplifyResult mesh;
  // 与 mesh.vertices 同序；未投射顶点保留在网格中并在此明确标记。
  std::vector<ProjectionResult> projections;
};

SimplifyFieldResult simplify_with_field(const std::vector<Vec3>& vertices,
                                        const std::vector<Tri>& triangles,
                                        const std::vector<double>& scalars,
                                        std::size_t target_face_count,
                                        double max_distance);

}  // namespace surface_field237
