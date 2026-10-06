#pragma once

#include <array>

#include "surface_field237/mesh_lod.h"

namespace surface_field237::internal {

// 点到非退化三角形的精确最近点，覆盖面内、边、顶点全部区域。
Vec3 closest_point_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b,
                               const Vec3& c);

// 三角形所在平面上一点的重心坐标，按 (a, b, c) 顺序返回。
std::array<double, 3> barycentric_coords(const Vec3& q, const Vec3& a,
                                         const Vec3& b, const Vec3& c);

}  // namespace surface_field237::internal
