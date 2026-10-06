#include "simplify.h"

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "qem.h"

namespace surface_field237::internal {

namespace {

struct Candidate {
  double cost = 0.0;
  std::uint32_t u = 0;  // 稳定原顶点 ID，u < v
  std::uint32_t v = 0;
  Vec3 pos = Vec3::Zero();

  bool operator<(const Candidate& o) const {
    if (cost != o.cost) return cost < o.cost;
    if (u != o.u) return u < o.u;
    return v < o.v;
  }
};

using CandidateSet = std::set<Candidate>;
using CandidateIndex = std::map<EdgeKey, CandidateSet::iterator>;

// 拓扑 link 条件：两端点共同邻居必须恰好等于该边邻接面的对顶点。
bool link_condition_ok(const MeshState& st, std::uint32_t u, std::uint32_t v) {
  const std::set<std::uint32_t> nu = st.neighbors(u);
  const std::set<std::uint32_t> nv = st.neighbors(v);
  std::set<std::uint32_t> common;
  std::set_intersection(nu.begin(), nu.end(), nv.begin(), nv.end(),
                        std::inserter(common, common.begin()));
  std::set<std::uint32_t> expected;
  for (std::uint32_t w : st.opposite_vertices(u, v)) expected.insert(w);
  return common == expected;
}

// 折叠 u -> v（v 移动到 new_pos）的合法性检查。
bool collapse_legal(const MeshState& st, std::uint32_t u, std::uint32_t v,
                    const Vec3& new_pos) {
  if (st.locked[u] || st.locked[v]) return false;  // 边界顶点不移动不删除
  if (!link_condition_ok(st, u, v)) return false;

  std::set<std::uint32_t> edge_faces;
  auto eit = st.edges.find(edge_key(u, v));
  if (eit == st.edges.end()) return false;
  for (int i = 0; i < eit->second.count; ++i)
    edge_faces.insert(eit->second.faces[i]);

  // 不删除整分量：折叠后两端点须仍有剩余邻面。
  std::set<std::uint32_t> remaining = st.v_faces[u];
  remaining.insert(st.v_faces[v].begin(), st.v_faces[v].end());
  for (std::uint32_t fid : edge_faces) remaining.erase(fid);
  if (remaining.empty()) return false;

  // u 的剩余邻面：重映射 u->v 后不得退化、翻转或重复。
  for (std::uint32_t fid : st.v_faces[u]) {
    if (edge_faces.count(fid)) continue;
    Tri t = st.faces[fid].v;
    for (auto& x : t)
      if (x == u) x = v;
    const Vec3& p0 = (t[0] == v) ? new_pos : st.pos[t[0]];
    const Vec3& p1 = (t[1] == v) ? new_pos : st.pos[t[1]];
    const Vec3& p2 = (t[2] == v) ? new_pos : st.pos[t[2]];
    if (tri_degenerate(p0, p1, p2)) return false;
    const Vec3 old_n = (st.pos[st.faces[fid].v[1]] - st.pos[st.faces[fid].v[0]])
                           .cross(st.pos[st.faces[fid].v[2]] -
                                  st.pos[st.faces[fid].v[0]]);
    const Vec3 new_n = (p1 - p0).cross(p2 - p0);
    if (new_n.dot(old_n) <= 0.0) return false;  // 法线翻转
    if (st.face_keys.count(face_key(t))) return false;  // 重复面
  }
  // v 的邻面：位置移动后不得退化或翻转。
  for (std::uint32_t fid : st.v_faces[v]) {
    if (edge_faces.count(fid)) continue;
    const Tri& t = st.faces[fid].v;
    bool has_u = false;
    for (auto x : t)
      if (x == u) has_u = true;
    if (has_u) continue;  // 已在上面检查
    const Vec3& p0 = (t[0] == v) ? new_pos : st.pos[t[0]];
    const Vec3& p1 = (t[1] == v) ? new_pos : st.pos[t[1]];
    const Vec3& p2 = (t[2] == v) ? new_pos : st.pos[t[2]];
    if (tri_degenerate(p0, p1, p2)) return false;
    const Vec3 old_n = (st.pos[t[1]] - st.pos[t[0]]).cross(st.pos[t[2]] - st.pos[t[0]]);
    const Vec3 new_n = (p1 - p0).cross(p2 - p0);
    if (new_n.dot(old_n) <= 0.0) return false;
  }
  return true;
}

// 执行折叠 u -> v：删除边邻接面，u 的其余邻面重映射到 v，合并误差矩阵。
void apply_collapse(MeshState& st, std::uint32_t u, std::uint32_t v,
                    const Vec3& new_pos) {
  auto eit = st.edges.find(edge_key(u, v));
  std::vector<std::uint32_t> edge_faces;
  for (int i = 0; i < eit->second.count; ++i)
    edge_faces.push_back(eit->second.faces[i]);
  for (std::uint32_t fid : edge_faces) st.remove_face(fid);

  std::vector<std::uint32_t> remap(st.v_faces[u].begin(), st.v_faces[u].end());
  for (std::uint32_t fid : remap) {
    st.remove_face(fid);
    for (auto& x : st.faces[fid].v)
      if (x == u) x = v;
    st.add_face(fid);
  }

  st.pos[v] = new_pos;
  st.quad[v] += st.quad[u];  // 保留合并矩阵，不重新初始化
  st.v_alive[u] = false;
  st.v_faces[u].clear();
}

// 为边计算候选；合法则插入候选集。
void evaluate_edge(MeshState& st, const EdgeKey& key, CandidateSet& cands,
                   CandidateIndex& index) {
  const std::uint32_t a = key.first;
  const std::uint32_t b = key.second;
  if (st.locked[a] || st.locked[b]) return;
  Candidate c;
  c.u = a;
  c.v = b;
  c.cost = quadric_optimal(st.quad[a] + st.quad[b], st.pos[a], st.pos[b], c.pos);
  if (!collapse_legal(st, a, b, c.pos)) return;
  auto [it, inserted] = cands.insert(c);
  if (inserted) index[key] = it;
}

void erase_edge_candidate(const EdgeKey& key, CandidateSet& cands,
                          CandidateIndex& index) {
  auto it = index.find(key);
  if (it == index.end()) return;
  cands.erase(it->second);
  index.erase(it);
}

}  // namespace

SimplifyOutcome run_simplification(MeshState& st, std::size_t target) {
  SimplifyOutcome out;
  CandidateSet cands;
  CandidateIndex index;

  for (const auto& [key, rec] : st.edges)
    evaluate_edge(st, key, cands, index);

  while (st.alive_faces > target) {
    if (cands.empty()) return out;  // 无合法候选，保留当前网格
    const Candidate best = *cands.begin();

    // 折叠前清除所有触及 u/v 的候选，失效代价不得继续使用。
    for (const auto& [ekey, rec] : st.edges)
      if (ekey.first == best.u || ekey.second == best.u ||
          ekey.first == best.v || ekey.second == best.v)
        erase_edge_candidate(ekey, cands, index);

    apply_collapse(st, best.u, best.v, best.pos);
    ++out.collapse_count;

    // 仅刷新受影响邻域的候选。
    std::set<std::uint32_t> region = st.neighbors(best.v);
    region.insert(best.v);
    std::vector<EdgeKey> dirty;
    for (const auto& [ekey, rec] : st.edges)
      if (region.count(ekey.first) || region.count(ekey.second))
        dirty.push_back(ekey);
    for (const EdgeKey& ekey : dirty) erase_edge_candidate(ekey, cands, index);
    for (const EdgeKey& ekey : dirty) evaluate_edge(st, ekey, cands, index);
  }
  out.target_reached = true;
  return out;
}

}  // namespace surface_field237::internal
