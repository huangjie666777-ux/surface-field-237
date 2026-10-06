#include "qem.h"

namespace surface_field237::internal {

namespace {

double quadric_error(const Eigen::Matrix4d& q, const Vec3& p) {
  const Eigen::Vector4d v(p.x(), p.y(), p.z(), 1.0);
  return v.dot(q * v);
}

}  // namespace

void init_quadrics(MeshState& st) {
  st.quad.assign(st.pos.size(), Eigen::Matrix4d::Zero());
  for (const auto& f : st.faces) {
    if (!f.alive) continue;
    const Vec3& a = st.pos[f.v[0]];
    const Vec3& b = st.pos[f.v[1]];
    const Vec3& c = st.pos[f.v[2]];
    Vec3 n = (b - a).cross(c - a);
    n.normalize();  // 输入已校验非零面积
    const double d = -n.dot(a);
    const Eigen::Vector4d plane(n.x(), n.y(), n.z(), d);
    const Eigen::Matrix4d kp = plane * plane.transpose();
    for (int k = 0; k < 3; ++k) st.quad[f.v[k]] += kp;
  }
}

double quadric_optimal(const Eigen::Matrix4d& q, const Vec3& pa, const Vec3& pb,
                       Vec3& out_pos) {
  const Eigen::Matrix3d a = q.block<3, 3>(0, 0);
  const Eigen::Vector3d b = q.block<3, 1>(0, 3);
  Eigen::FullPivLU<Eigen::Matrix3d> lu(a);
  lu.setThreshold(1e-10);
  if (lu.isInvertible()) {
    const Vec3 p = lu.solve(-b);
    if (p.allFinite()) {
      out_pos = p;
      return quadric_error(q, p);
    }
  }
  // 奇异：比较两端点与中点。
  const Vec3 mid = 0.5 * (pa + pb);
  const double ea = quadric_error(q, pa);
  const double eb = quadric_error(q, pb);
  const double em = quadric_error(q, mid);
  if (ea <= eb && ea <= em) {
    out_pos = pa;
    return ea;
  }
  if (eb <= em) {
    out_pos = pb;
    return eb;
  }
  out_pos = mid;
  return em;
}

}  // namespace surface_field237::internal
