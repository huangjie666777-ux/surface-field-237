#include "surface_field237/mesh_lod.h"

#include <stdexcept>
#include <unordered_map>

#include "qem.h"
#include "simplify.h"
#include "validation.h"

namespace surface_field237 {

std::string to_string(StopReason reason) {
  switch (reason) {
    case StopReason::kAlreadyAtOrBelowTarget:
      return "already_at_or_below_target";
    case StopReason::kTargetReached:
      return "target_reached";
    case StopReason::kNoValidCandidate:
      return "no_valid_candidate";
  }
  return "unknown";
}

SimplifyResult simplify(const std::vector<Vec3>& vertices,
                        const std::vector<Tri>& triangles,
                        std::size_t target_face_count) {
  if (target_face_count == 0)
    throw std::invalid_argument("surface_field237: target face count must be positive");

  SimplifyResult result;
  // 先校验输入，再判断目标是否已满足。
  internal::MeshState st = internal::build_validated_state(vertices, triangles);

  if (target_face_count >= triangles.size()) {
    // 目标不小于当前面数：原样返回，不改输入。
    result.vertices = vertices;
    result.triangles = triangles;
    result.actual_face_count = triangles.size();
    result.collapse_count = 0;
    result.stop_reason = StopReason::kAlreadyAtOrBelowTarget;
    return result;
  }

  internal::init_quadrics(st);
  const internal::SimplifyOutcome outcome =
      internal::run_simplification(st, target_face_count);

  // 紧凑输出：重映射活顶点与三角索引。
  std::unordered_map<std::uint32_t, std::uint32_t> remap;
  remap.reserve(st.pos.size());
  for (std::uint32_t i = 0; i < st.pos.size(); ++i) {
    if (!st.v_alive[i]) continue;
    remap[i] = static_cast<std::uint32_t>(result.vertices.size());
    result.vertices.push_back(st.pos[i]);
  }
  for (const auto& f : st.faces) {
    if (!f.alive) continue;
    result.triangles.push_back(Tri{remap[f.v[0]], remap[f.v[1]], remap[f.v[2]]});
  }
  result.actual_face_count = result.triangles.size();
  result.collapse_count = outcome.collapse_count;
  result.stop_reason = outcome.target_reached
                           ? StopReason::kTargetReached
                           : StopReason::kNoValidCandidate;
  return result;
}

}  // namespace surface_field237
