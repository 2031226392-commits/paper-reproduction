# Simulation 1 扩展：Nested Cross-Validation（NCV）选模方案

## 1. 重点阅读结论（来自 NCV.pdf 的核心思想）

结合 `NCV.pdf` Figure 7 与第 4.3 节（Holdout MSE identity / Nested CV），可提炼出与本项目最相关的原则：

- 外层 holdout（论文图中蓝色折）在内层调参过程中完全不可见。
- 内层仅使用训练部分做 CV/调参；外层只用于独立检验。
- 通过多个外层切分重复，得到更稳健的外层检验结果。

你给出的三折思路（3 折轮流 holdout，内层用剩余 2 折选参数，再到 holdout 检验）与论文“外层新鲜 holdout”思想一致。

## 2. 本项目约束下的 NCV 目标

- 任务：沿用 Simulation 1 的 16 个候选子模型做选模正确率比较。
- 外层：3-fold（每轮 1 折 holdout）。
- 内层：在剩余 2 折上做 2-fold CV 以选模型。
- 关键约束：每一轮内层选模不能使用该轮 holdout 的任何信息。

## 3. 方法定义（含你要求的 ncv1 / ncv2 / ncv3）

记外层三轮得到三组候选 `(m1, m2, m3)`，对应外层 holdout 误差 `(e1, e2, e3)`。

### 3.1 共同定义：模型“正则化强弱”

由于 Simulation 1 是离散子模型，不是连续 λ，定义：

- “正则化更强” = 模型更简约 = 变量数更少（`|J_m|` 更小）。
- 若变量数相同，按该层对应的 CV 误差更小优先；再同分按模型编号最小。

### 3.2 内层选择规则

- `inner_best`：内层 2-fold CV 平均误差最小。
- `inner_1se`：在内层最小误差 + 1SE 阈值内，选“正则化更强”的模型。

### 3.3 外层聚合规则

- `outer_best`：三轮外层误差 `e1,e2,e3` 中最小者对应模型。
- `outer_1se`：在 `min(e) + 1SE_outer` 范围内，选“正则化更强”的模型。
  这里 `1SE_outer` 用三轮误差的样本标准差除以 `sqrt(3)`。

### 3.4 三种 NCV 方法（项目命名）

- `ncv1`：`inner_1se + outer_best`
- `ncv2`：`inner_best + outer_1se`
- `ncv3`：`inner_1se + outer_1se`

## 4. 其他对照方法（同图比较）

- `cv`：现有 5-fold CV 最小误差选模（可复用已跑结果）。
- `cv_1se`：5-fold CV 的 1SE 规则（最小误差+1SE范围内取最简约模型）。
- `cvc`：沿用当前 Simulation 1 配置（5-fold，`α=0.05`，`α'=0.005`，`B=200`）。

说明：按你的要求，`cv` 和 `cvc` 可先复用已有结果；新增只需跑 `cv_1se` 与 `ncv1/2/3`。

## 5. 实验设置（沿用 Simulation 1）

- 两组真值 `β`：Model 1 / Model 2。
- `n = 40, 80, 160, 320, 640`。
- `σ² = 1, 4`。
- 候选模型 16 个（含截距）。
- 重复次数：先冒烟（每格 10 次）后全量（每格 100 次）。

## 6. 输出与绘图方案（正确率）

## 6.1 原始结果表

新增建议：`results/sim1_ncv_raw.csv`

字段建议：
- `method`（cv/cv_1se/cvc/ncv1/ncv2/ncv3）
- `model_name`, `sigma2`, `n`, `rep`
- `selected_model_id`, `selected_size`, `true_model_id`, `correct`

## 6.2 汇总表

新增建议：`results/sim1_ncv_summary.csv`

- 按 `method × model_name × sigma2 × n` 聚合：
  - `selection_rate = mean(correct)`
  - 可附 `se = sd(correct)/sqrt(reps)`

## 6.3 图片

新增建议：`figures/sim1_ncv_selection_rate.png`

- 2×2 分面（Model1/2 × σ²=1/4）。
- x 轴：`n`。
- y 轴：正确率。
- 线条：`cv`, `cv_1se`, `cvc`, `ncv1`, `ncv2`, `ncv3`。
- 可选加误差棒（±1SE）。

## 7. 可行性评估

- 统计可行：方法定义完整，且满足“内层调参与外层 holdout 隔离”。
- 计算可行：NCV 比普通 CV/CVC 更耗时，但在 16 个离散模型下可承受。
- 风险点：
  - 外层只有 3 轮，`outer_1se` 的标准误估计不稳定；需在报告中注明。
  - `cv/cvc` 若直接复用历史结果，需确保随机种子与本轮新增方法可比较；
    更严格做法是六种方法同一批数据重跑。

## 8. 建议执行顺序

1. 实现 `cv_1se` 与 `ncv1/2/3` 函数（先 n=40、rep=10 冒烟）。
2. 复用已有 `cv/cvc` 结果先画一版总览图验证管线。
3. 若要严格可比，六种方法在同一随机数据上全量重跑 100 次并覆盖图表。

