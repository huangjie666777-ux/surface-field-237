#pragma once

#include "mesh_state.h"

namespace surface_field237::internal {

// 由初始面单位法向平面建立每顶点二次误差矩阵。
void init_quadrics(MeshState& st);

// 端点矩阵相加后的最小误差位置；奇异时比较两端点与中点。
// 返回误差值，位置写入 out_pos。
double quadric_optimal(const Eigen::Matrix4d& q, const Vec3& pa, const Vec3& pb,
                       Vec3& out_pos);

}  // namespace surface_field237::internal

