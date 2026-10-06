# surface_field237

基于二次误差度量（QEM）边折叠的三角网格轻量化库，并支持原扫描表面
标量场投射。C++20 + Eigen 3.4.0，命名空间 `surface_field237`，
仅依赖 `third_party/eigen3` 头文件，不含前端或 HTTP 服务。

## 构建与运行

```sh
make lib      # 构建静态库 build/libsurface_field237.a
make test     # 构建并运行自测 tests/test_main.cpp
make example  # 构建并运行示例（开口曲面减面 + 标量场投射）
make check    # test + example
```

## 减面接口

```cpp
#include "surface_field237/mesh_lod.h"

surface_field237::SimplifyResult r = surface_field237::simplify(vertices, triangles, target);
```

- 输入：双精度三维顶点 `std::vector<Eigen::Vector3d>`、有向三角索引
  `std::vector<std::array<uint32_t,3>>`、正整数目标面数。输入按值语义处理，不被修改。
- 非法输入抛出 `std::invalid_argument`：非有限坐标、索引越界、重复面、
  零面积面、共享边绕向不一致、非流形边（>2 面共边）、非流形顶点
  （邻域扇不连通或分叉）、目标面数为 0。
- 输出 `SimplifyResult`：紧凑顶点与重映射三角索引、`actual_face_count`、
  `collapse_count`、`stop_reason`
  （`already_at_or_below_target` / `target_reached` / `no_valid_candidate`）。
- 目标面数不小于当前面数时原样返回（输入仍先校验）。
- 减面输出仅保留被存活面引用的顶点，不保留未引用顶点。

## 标量场投射接口

```cpp
#include "surface_field237/surface_field.h"

surface_field237::SurfaceField field(vertices, triangles, scalars);
auto hits = field.project(points, max_distance);

auto r = surface_field237::simplify_with_field(
    vertices, triangles, scalars, target, max_projection_distance);
```

- `SurfaceField` 由源顶点、有向三角索引和逐顶点标量建立可复用查询对象；
  网格合法性与减面入口同一套规则（允许开口与多分量），标量须长度匹配
  且全部有限。对象持有独立快照，构造后修改输入数组不影响查询。
- 内部构建三角形 AABB 二叉层次索引（质心最长轴中位切分，叶大小 4），
  以包围盒距离下界剪枝，叶内精确计算点到三角面的最近点（面内、边、
  顶点全区域），不以最近顶点代替、不做全表扫描；索引供多批查询复用。
- `project` 批量接收同坐标系查询点与非负有限上限，返回
  `ProjectionResult`：最近位置、原面 ID、按该面索引顺序的重心坐标、
  欧氏距离与重心插值标量，按输入顺序输出。等距选原面 ID 较小者；
  距离等于上限可接收；超限点 `projected=false`，数值字段为 NaN、
  `face_id` 为 `UINT32_MAX`，不以零值冒充。任一点非有限或上限非法时
  整批抛 `std::invalid_argument`，对象状态不受影响。
- `simplify_with_field` 减面后以最终输出顶点查询原表面，共用同一索引
  与插值；`projections` 与输出顶点一一对应，未投射顶点保留并以
  `projected=false` 明确标记。输入按值语义处理，不被修改。

## 容差

- 退化三角形判定：双倍面积 ≤ 1e-12 × 最长边²（相对容差），
  边长 1e-7 的有效三角形不会被误拒。
- 投射等距与上限判定使用精确平方距离比较：距离 ≤ 上限（含等号）即接收，
  等距时取原面 ID 较小者。
- 二次误差法方程求解阈值：FullPivLU 阈值 1e-10，奇异时比较两端点
  与中点取误差最小者。

## 算法要点

- 由初始面单位法向平面建立每顶点 4x4 二次误差矩阵；候选仅限当前网格边。
- 端点矩阵相加求最小误差位置（FullPivLU 解 3x3 法方程）；奇异时比较两端点
  与中点取误差最小者。折叠后保留合并矩阵，不重新初始化。
- 按误差最小的合法候选推进，代价并列时按稳定原顶点 ID 升序。
- 初始边界顶点全部锁定：不移动、不删除，边界坐标逐位保持；折叠不连接
  不同连通分量。
- 每次折叠检查拓扑 link 条件，删除边邻接面并重映射其余邻面；拒绝产生
  重复面、零面积面、法线翻转或删除整分量的折叠。
- 候选集在每次折叠后局部失效并重建（折叠前清除所有触及两端点的候选，
  折叠后仅重估受影响邻域），不会使用失效代价。
- 面数不大于目标即停止；无合法候选时保留当前网格并以
  `no_valid_candidate` 明确未达目标。
- 不提供全局自相交防护，不保证全局最优近似。

## 源码组织（按职责划分）

| 文件 | 职责 |
| --- | --- |
| `include/surface_field237/mesh_lod.h` | 公共调用接口与结果类型 |
| `src/mesh_lod.cpp` | 接口实现：校验编排、早退、输出紧凑化 |
| `src/validation.cpp` | 输入校验与初始拓扑构建、边界锁定 |
| `src/topology.cpp` | 网格状态与邻接维护（边/面/顶点邻接） |
| `src/qem.cpp` | 二次误差矩阵初始化与最优位置求解 |
| `src/simplify.cpp` | 候选评估、合法性检查、折叠与迭代推进 |
| `include/surface_field237/surface_field.h` | 标量场投射公共接口 |
| `src/closest_point.cpp` | 点到三角形精确最近点与重心坐标 |
| `src/bvh.cpp` | 三角形 AABB 层次索引构建与最近面查询 |
| `src/surface_field.cpp` | 查询对象、批量投射与携带投射的减面入口 |

## 示例输出

`examples/open_surface.cpp` 对 12x12 顶点的开口波浪网格减面到 60：

```
input : 144 vertices, 242 triangles
output: 53 vertices, 60 triangles (target 60)
collapses: 91
stop reason: target_reached
boundary vertices: 44 -> 44, coords preserved: yes
```
