#pragma once

#include "mesh_state.h"

namespace surface_field237::internal {

struct SimplifyOutcome {
  std::size_t collapse_count = 0;
  bool target_reached = false;
};

// 按误差最小合法候选迭代折叠，直到面数不大于目标或无合法候选。
SimplifyOutcome run_simplification(MeshState& st, std::size_t target);

}  // namespace surface_field237::internal

