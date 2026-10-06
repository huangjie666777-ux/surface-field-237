# surface_field237

基于二次误差度量（QEM）边折叠的三角网格轻量化库。C++20 + Eigen 3.4.0，
命名空间 `surface_field237`，仅依赖 `third_party/eigen3` 头文件。

## 构建与运行

```sh
make lib      # 构建静态库 build/libsurface_field237.a
make test     # 构建并运行自测 tests/test_main.cpp
make example  # 构建并运行开口曲面示例 examples/open_surface.cpp
make check    # test + example
```

## 接口

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

## 示例输出

`examples/open_surface.cpp` 对 12x12 顶点的开口波浪网格减面到 60：

```
input : 144 vertices, 242 triangles
output: 53 vertices, 60 triangles (target 60)
collapses: 91
stop reason: target_reached
boundary vertices: 44 -> 44, coords preserved: yes
```

