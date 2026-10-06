#include "surface_field237/surface_field.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "bvh.h"
#include "closest_point.h"
#include "validation.h"

namespace surface_field237 {

struct SurfaceField::Impl {
  std::vector<Vec3> vertices;
  std::vector<Tri> triangles;
  std::vector<double> scalars;
  internal::TriBVH bvh;

  Impl(std::vector<Vec3> v, std::vector<Tri> t, std::vector<double> s)
      : vertices(std::move(v)),
        triangles(std::move(t)),
        scalars(std::move(s)),
        bvh(vertices, triangles) {}
};

SurfaceField::SurfaceField(std::vector<Vec3> vertices,
                           std::vector<Tri> triangles,
                           std::vector<double> scalars) {
  // 网格沿用减面入口的合法性规则（允许开口与多分量）。
  internal::build_validated_state(vertices, triangles);
  if (scalars.size() != vertices.size())
    throw std::invalid_argument(
        "surface_field237: scalar count must match vertex count");
  for (std::size_t i = 0; i < scalars.size(); ++i)
    if (!std::isfinite(scalars[i]))
      throw std::invalid_argument("surface_field237: non-finite scalar at " +
                                  std::to_string(i));
  // 快照在构造后独立于调用方持有的数组。
  impl_ = std::make_unique<Impl>(std::move(vertices), std::move(triangles),
                                 std::move(scalars));
}

SurfaceField::~SurfaceField() = default;
SurfaceField::SurfaceField(SurfaceField&&) noexcept = default;
SurfaceField& SurfaceField::operator=(SurfaceField&&) noexcept = default;

const std::vector<Vec3>& SurfaceField::vertices() const {
  return impl_->vertices;
}
const std::vector<Tri>& SurfaceField::triangles() const {
  return impl_->triangles;
}
const std::vector<double>& SurfaceField::scalars() const {
  return impl_->scalars;
}

std::vector<ProjectionResult> SurfaceField::project(
    const std::vector<Vec3>& points, double max_distance) const {
  // 整批校验：任一非法即拒绝，对象状态不受影响。
  if (!std::isfinite(max_distance) || max_distance < 0.0)
    throw std::invalid_argument(
        "surface_field237: max projection distance must be non-negative and "
        "finite");
  for (std::size_t i = 0; i < points.size(); ++i) {
    const Vec3& p = points[i];
    if (!std::isfinite(p.x()) || !std::isfinite(p.y()) ||
        !std::isfinite(p.z()))
      throw std::invalid_argument(
          "surface_field237: non-finite query point at " + std::to_string(i));
  }

  std::vector<ProjectionResult> out(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    double dist = 0.0;
    const std::uint32_t fid =
        impl_->bvh.nearest_face(points[i], max_distance, &dist);
    if (fid == UINT32_MAX) continue;  // 超限：保持未投射标记与 NaN 字段

    const Tri& t = impl_->triangles[fid];
    const Vec3& a = impl_->vertices[t[0]];
    const Vec3& b = impl_->vertices[t[1]];
    const Vec3& c = impl_->vertices[t[2]];
    const Vec3 q = internal::closest_point_on_triangle(points[i], a, b, c);
    const auto w = internal::barycentric_coords(q, a, b, c);

    ProjectionResult& r = out[i];
    r.projected = true;
    r.closest_point = q;
    r.face_id = fid;
    r.barycentric = w;
    r.distance = dist;
    r.scalar = w[0] * impl_->scalars[t[0]] + w[1] * impl_->scalars[t[1]] +
               w[2] * impl_->scalars[t[2]];
  }
  return out;
}

SimplifyFieldResult simplify_with_field(
    const std::vector<Vec3>& vertices, const std::vector<Tri>& triangles,
    const std::vector<double>& scalars, std::size_t target_face_count,
    double max_projection_distance) {
  // 先建立查询对象（含校验与空间索引），再减面，最后共用同一索引投射。
  const SurfaceField field(vertices, triangles, scalars);
  SimplifyFieldResult out;
  out.mesh =
      simplify(field.vertices(), field.triangles(), target_face_count);
  out.projections = field.project(out.mesh.vertices, max_projection_distance);
  return out;
}

}  // namespace surface_field237
