#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace surface_field237 {

using Vec3 = Eigen::Vector3d;
using Tri = std::array<std::uint32_t, 3>;

enum class StopReason {
  kAlreadyAtOrBelowTarget,
  kTargetReached,
  kNoValidCandidate,
};

std::string to_string(StopReason reason);

struct SimplifyResult {
  std::vector<Vec3> vertices;
  std::vector<Tri> triangles;
  std::size_t actual_face_count = 0;
  std::size_t collapse_count = 0;
  StopReason stop_reason = StopReason::kNoValidCandidate;
};

SimplifyResult simplify(const std::vector<Vec3>& vertices,
                        const std::vector<Tri>& triangles,
                        std::size_t target_face_count);

}  // namespace surface_field237
