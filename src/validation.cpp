#include "validation.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>

namespace surface_field237::internal {

namespace {

[[noreturn]] void fail(const std::string& msg) {
  throw std::invalid_argument("surface_field237: " + msg);
}

// 非流形顶点检测：顶点邻接面在其对边上形成的邻域图
// 必须是单条路径（开口扇）或单个环（闭合扇）。
void check_manifold_vertices(const MeshState& st) {
  for (std::uint32_t x = 0; x < st.pos.size(); ++x) {
    if (st.v_faces[x].empty()) continue;
    std::map<std::uint32_t, std::vector<std::uint32_t>> adj;
    for (std::uint32_t fid : st.v_faces[x]) {
      const auto& fv = st.faces[fid].v;
      std::uint32_t opp[2];
      int n = 0;
      for (int k = 0; k < 3; ++k)
        if (fv[k] != x) opp[n++] = fv[k];
      adj[opp[0]].push_back(opp[1]);
      adj[opp[1]].push_back(opp[0]);
    }
    for (const auto& [nb, lst] : adj)
      if (lst.size() > 2)
        fail("non-manifold vertex " + std::to_string(x));
    // 连通性：从任一邻居出发遍历应覆盖全部邻居。
    std::set<std::uint32_t> seen;
    std::vector<std::uint32_t> stack{adj.begin()->first};
    while (!stack.empty()) {
      std::uint32_t cur = stack.back();
      stack.pop_back();
      if (!seen.insert(cur).second) continue;
      for (std::uint32_t nb : adj[cur]) stack.push_back(nb);
    }
    if (seen.size() != adj.size())
      fail("non-manifold vertex " + std::to_string(x));
  }
}

}  // namespace

MeshState build_validated_state(const std::vector<Vec3>& vertices,
                                const std::vector<Tri>& triangles) {
  if (vertices.empty()) fail("empty vertex array");
  if (triangles.empty()) fail("empty triangle array");

  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const Vec3& p = vertices[i];
    if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
      fail("non-finite coordinate at vertex " + std::to_string(i));
  }

  MeshState st;
  st.pos = vertices;
  st.v_alive.assign(vertices.size(), true);
  st.locked.assign(vertices.size(), false);
  st.v_faces.resize(vertices.size());
  st.faces.resize(triangles.size());

  std::set<FaceKey> seen_faces;
  for (std::size_t i = 0; i < triangles.size(); ++i) {
    const Tri& t = triangles[i];
    for (int k = 0; k < 3; ++k) {
      if (t[k] >= vertices.size())
        fail("triangle " + std::to_string(i) + " index out of range");
    }
    if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
      fail("triangle " + std::to_string(i) + " has repeated indices");
    if (!seen_faces.insert(face_key(t)).second)
      fail("duplicate triangle " + std::to_string(i));
    if (tri_degenerate(vertices[t[0]], vertices[t[1]], vertices[t[2]]))
      fail("zero-area triangle " + std::to_string(i));
    st.faces[i].v = t;
  }

  // 组边并检查非流形边与绕向一致性。
  std::map<EdgeKey, std::vector<std::pair<std::uint32_t, bool>>> directed;
  for (std::size_t i = 0; i < triangles.size(); ++i) {
    const Tri& t = triangles[i];
    for (int k = 0; k < 3; ++k) {
      const std::uint32_t a = t[k];
      const std::uint32_t b = t[(k + 1) % 3];
      directed[edge_key(a, b)].push_back(
          {static_cast<std::uint32_t>(i), a < b});
    }
  }
  for (const auto& [key, lst] : directed) {
    if (lst.size() > 2)
      fail("non-manifold edge (" + std::to_string(key.first) + ", " +
           std::to_string(key.second) + ")");
    if (lst.size() == 2 && lst[0].second == lst[1].second)
      fail("inconsistent winding on edge (" + std::to_string(key.first) +
           ", " + std::to_string(key.second) + ")");
  }

  for (std::size_t i = 0; i < triangles.size(); ++i)
    st.add_face(static_cast<std::uint32_t>(i));

  check_manifold_vertices(st);

  // 锁定初始边界顶点。
  for (const auto& [key, rec] : st.edges) {
    if (rec.count == 1) {
      st.locked[key.first] = true;
      st.locked[key.second] = true;
    }
  }
  return st;
}

}  // namespace surface_field237::internal

