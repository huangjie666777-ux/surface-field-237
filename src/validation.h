#pragma once

#include "mesh_state.h"

namespace surface_field237::internal {

// 校验输入并构建初始网格状态；非法输入抛 std::invalid_argument。
MeshState build_validated_state(const std::vector<Vec3>& vertices,
                                const std::vector<Tri>& triangles);

}  // namespace surface_field237::internal

