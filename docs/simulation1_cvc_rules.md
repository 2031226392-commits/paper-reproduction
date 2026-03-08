# Simulation 1 CVC 细化规则记录

## 默认实验设置（当前实现）

- 5-fold CV。
- CVC 参数：α = 0.05（单一默认值）、α' = 0.005、B = 200。
- 候选模型：含截距的 16 个子模型。
- 对比方法：cvc<sub>pmax</sub>、cvc<sub>frac</sub>、cv、bic、fy。

## 你关心的并列最简模型（tie）处理

记 A<sub>cv</sub> = {m : p̂<sub>cv,m</sub> ≥ α}。

1. 先找最简模型集合：
M<sub>min</sub> = {m in A<sub>cv</sub> : |J<sub>m</sub>| = min<sub>k in A<sub>cv</sub></sub> |J<sub>k</sub>|}。

2. 两种 CVC 口径：
- cvc<sub>pmax</sub>：从 M<sub>min</sub> 中选 p̂<sub>cv,m</sub> 最大的那个作为最终模型。
  若仍并列，则按更小 CV 风险、再按更小模型编号打破平局。
  正确率按 0/1 计。
- cvc<sub>frac</sub>：不再强行单选；若 |M<sub>min</sub>| = K，且真模型在 M<sub>min</sub> 中，
  该次重复记正确率 1/K，否则记 0。
  即“并列最简模型平均分摊正确率”。

## 正确率统计口径

- cv、bic、cvc<sub>pmax</sub>：每次重复是 0/1，最终对重复取均值。
- cvc<sub>frac</sub>：每次重复可取分数值（例如 0.25），最终对重复取均值。

## 备注

- 若某次 A<sub>cv</sub> 为空，按实现约定回退到 m<sub>cv</sub>，保证可选模型非空。
- fy 目前采用可落地近似实现：对每个子模型与全模型做嵌套 F 检验，
  不拒绝（p ≥ α）的模型进入 FY 置信集；并记录其集合大小。

## “无法拒绝模型数量”图的设定（CVC vs FY）

- 图片文件：`figures/sim1_alpha005_nonrejected_count.png`。
- 统计对象：每次重复中的“不可拒绝模型集合大小”。
- CVC 曲线：
  - 集合定义：A<sub>cv</sub> = {m : p̂<sub>cv,m</sub> ≥ α}；
  - 计数：|A<sub>cv</sub>|；
  - 图中使用方法标签 `cvc`（代码实现对应 cvc<sub>pmax</sub> 的同一 A<sub>cv</sub> 大小）。
- FY 曲线：
  - 以全模型为基准，对每个子模型做嵌套 F 检验；
  - 集合定义：A<sub>fy</sub> = {m : p<sub>fy,m</sub> ≥ α}；
  - 计数：|A<sub>fy</sub>|。
- 纵轴含义：在 100 次重复上取平均后的集合大小（Mean # Non-Rejected Models）。
- 分面：`model1/model2 × σ^2=1/4`，横轴是 `n = 40, 80, 160, 320, 640`。
