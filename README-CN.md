# Mini Geometric Control（迷你几何控制理论）

**从零开始、零依赖的 C 语言实现**，涵盖几何控制理论（Geometric Control Theory）——由 Brockett、Jurdjevic、Isidori、Nijmeijer、van der Schaft 等人奠基的、从微分几何视角研究非线性控制系统的理论框架。每个模块将经典几何控制形式化理论转化为可运行的 C 代码，在流形上的现代控制理论与计算实践之间架设桥梁。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-connection-principal-bundle](mini-connection-principal-bundle/) | Ehresmann 联络，联络 1-形式，曲率 2-形式，和乐群，平行移动，主 G-丛，格点规范理论（连接变量、Wilson 圈），规范变换 | MIT 18.755（李群与李代数），MIT 6.832 |
| [mini-differential-flatness](mini-differential-flatness/) | 微分平坦（Fliess–Levine–Martin–Rouchon 理论），平坦输出，相对阶，无需 ODE 积分的轨迹生成，Brunovsky 标准形，动态反馈线性化，多阶导数链 | MIT 6.832（欠驱动机器人），Stanford AA 203 |
| [mini-feedback-linearization-geo](mini-feedback-linearization-geo/) | 输入-输出 / 全状态反馈线性化，李导数与李括号，相对阶分析，Frobenius 定理与对合分布，零动态，非线性坐标变换，Isidori–Nijmeijer 框架 | MIT 6.241J（动力系统与控制），MIT 6.832，Stanford ENGR 205 |
| [mini-geometric-optimal-ctrl](mini-geometric-optimal-ctrl/) | 光滑流形上的控制仿射系统，切丛/余切丛，流形上的 Pontryagin 极大值原理，次黎曼几何，几何积分器（李群方法），可达集与可控性（Chow 定理） | MIT 6.832，Stanford AA 203，Caltech CDS 205 |
| [mini-hamiltonian-control](mini-hamiltonian-control/) | 辛流形上的哈密顿系统，泊松几何，Hamilton–Jacobi–Bellman（HJB）方程，能量整形（端口哈密顿系统），辛矩阵，相空间分析，庞特里亚金极大值原理（对偶形式） | MIT 2.151（高等动力学），Stanford AA 242 |
| [mini-lie-group-mechanics](mini-lie-group-mechanics/) | 力学中的李群（SO(3), SE(3), SU(2), SO(2), SE(2)），李代数结构常数，群作用与伴随/余伴随轨道，Euler–Poincaré 方程，Lie–Poisson 约化，李群上的几何积分 | MIT 2.032（动力学），MIT 2.151，Caltech CDS 205 |
| [mini-nonholonomic-planning](mini-nonholonomic-planning/) | 非完整约束与分布，对合性与 Frobenius 定理，链式系统转换（Brockett），运动学汽车 / 蛇形机器人 / 球-板系统，速度约束下的运动规划，Laumond–Murray 框架 | MIT 6.832，MIT 6.821（机器人操作） |
| [mini-symmetry-reduction](mini-symmetry-reduction/) | Marsden–Weinstein 辛约化，动量映射，余伴随轨道（Arnold），泊松约化，Euler–Poincaré 约化，自由/真确群作用，哈密顿约化，约化相空间 | MIT 18.755，Caltech CDS 205，MIT 2.151 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，与几何控制的原始文献对齐（Brockett, Jurdjevic, Isidori, Marsden, Agrachev & Sachkov 等）
- **从零实现** — 不使用现有库；每个概念从几何与控制论原语构建（李括号、微分形式、辛结构）

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-lie-group-mechanics
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-geometric-control/
├── mini-connection-principal-bundle/ # Ehresmann 联络，和乐，主丛上的曲率
├── mini-differential-flatness/       # 平坦输出，轨迹生成，Brunovsky 标准形
├── mini-feedback-linearization-geo/  # 基于李导数的反馈线性化，零动态
├── mini-geometric-optimal-ctrl/      # 流形上的 PMP，次黎曼几何，几何积分器
├── mini-hamiltonian-control/         # 辛几何，HJB 方程，能量整形
├── mini-lie-group-mechanics/         # SO(3)/SE(3) 运动学，Euler–Poincaré，Lie–Poisson 约化
├── mini-nonholonomic-planning/       # 分布，链式系统，约束下的运动规划
└── mini-symmetry-reduction/          # Marsden–Weinstein 约化，动量映射，余伴随轨道
```

## 许可证

MIT
