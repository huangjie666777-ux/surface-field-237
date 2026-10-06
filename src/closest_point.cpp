#include "closest_point.h"

namespace surface_field237::internal {

// Christer Ericson, Real-Time Collision Detection 5.1.5 区域判别法。
Vec3 closest_point_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b,
                               const Vec3& c) {
  const Vec3 ab = b - a;
  const Vec3 ac = c - a;
  const Vec3 ap = p - a;
  const double d1 = ab.dot(ap);
  const double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0) return a;  // 顶点 A 区域

  const Vec3 bp = p - b;
  const double d3 = ab.dot(bp);
  const double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3) return b;  // 顶点 B 区域

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    return a + ab * (d1 / (d1 - d3));  // 边 AB 区域

  const Vec3 cp = p - c;
  const double d5 = ab.dot(cp);
  const double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6) return c;  // 顶点 C 区域

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    return a + ac * (d2 / (d2 - d6));  // 边 AC 区域

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + (c - b) * w;  // 边 BC 区域
  }

  // 面内区域。
  const double denom = 1.0 / (va + vb + vc);
  return a + ab * (vb * denom) + ac * (vc * denom);
}

std::array<double, 3> barycentric_coords(const Vec3& q, const Vec3& a,
                                         const Vec3& b, const Vec3& c) {
  const Vec3 v0 = b - a;
  const Vec3 v1 = c - a;
  const Vec3 v2 = q - a;
  const double d00 = v0.dot(v0);
  const double d01 = v0.dot(v1);
  const double d11 = v1.dot(v1);
  const double d20 = v2.dot(v0);
  const double d21 = v2.dot(v1);
  const double denom = d00 * d11 - d01 * d01;  // 非退化三角形，不为零
  const double w1 = (d11 * d20 - d01 * d21) / denom;
  const double w2 = (d00 * d21 - d01 * d20) / denom;
  return {1.0 - w1 - w2, w1, w2};
}

}  // namespace surface_field237::internal
