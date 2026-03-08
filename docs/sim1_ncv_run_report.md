# Simulation 1 NCV 扩展运行记录

## 执行口径

- 外层：5-fold（轮流 holdout）。
- 内层：4-fold（仅在外层训练子集内调参）。
- 方法：`cv`, `cv_1se`, `cvc`, `ncv1`, `ncv2`, `ncv3`。
- 其余设定：与 Simulation 1 主线一致（`β` 两组、`σ²=1/4`、`n=40..640`、Gaussian）。

## 运行命令

- 冒烟：`Rscript scripts/run_sim1_ncv.R --smoke=true`
- 全量：`Rscript scripts/run_sim1_ncv.R`

## 输出文件

- 冒烟原始：`results/sim1_ncv_smoke_raw.csv`
- 冒烟汇总：`results/sim1_ncv_smoke_summary.csv`
- 冒烟图：`figures/sim1_ncv_smoke_selection_rate.png`

- 全量原始：`results/sim1_ncv_raw.csv`
- 全量汇总：`results/sim1_ncv_summary.csv`
- 全量图：`figures/sim1_ncv_selection_rate.png`

## 备注

- `cv` 与 `cvc` 在该扩展里与新增方法同批次数据一起重跑，保证方法间可比性。
- `ncv1/ncv2/ncv3` 定义见：`docs/ncv_sim1_extension_plan.md`。
