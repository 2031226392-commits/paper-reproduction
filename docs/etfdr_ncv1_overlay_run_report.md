# ETFDR Figure 1 + NCV1 Overlay 运行记录

## 目标

在 ETFDR Figure 1 设定下，复用已有 `cv/cvc` 结果，仅新增 `ncv1` 与 `ncv1_1se` 选点，
并将 `cv/cvc/ncv1` 以及 `cv_1se/ncv1_1se` 竖线绘制到同一张图。

## 折数设置

- 外层：5-fold holdout。
- 内层：4-fold（仅外层训练数据内调参）。

## 运行命令

- 冒烟：`Rscript scripts/run_etfdr_ncv1_overlay.R --smoke=true`
- 全量：`Rscript scripts/run_etfdr_ncv1_overlay.R`

## 输出

- 冒烟图：`figures/etfdr_fig1_smoke_ncv1_overlay_figure.png`
- 全量图：`figures/etfdr_fig1_ncv1_overlay_figure.png`

- 冒烟选点：`results/etfdr_fig1_smoke_ncv1_overlay_ncv1_selected.csv`
- 全量选点：`results/etfdr_fig1_ncv1_overlay_ncv1_selected.csv`

- 冒烟全方法合并：`results/etfdr_fig1_smoke_ncv1_overlay_all_selected.csv`
- 全量全方法合并：`results/etfdr_fig1_ncv1_overlay_all_selected.csv`
