# ETFDR Figure 1 复现执行记录

## 运行命令

- 冒烟：`Rscript scripts/run_etfdr_fig1.R --smoke=true`
- 全量：`Rscript scripts/run_etfdr_fig1.R`

## 输出文件

- 冒烟主图：`figures/etfdr_fig1_smoke_figure1_repro.png`
- 冒烟对照图：`figures/etfdr_fig1_smoke_cv_cvc_fdr_mse.png`
- 全量主图：`figures/etfdr_fig1_figure1_repro.png`
- 全量对照图：`figures/etfdr_fig1_cv_cvc_fdr_mse.png`

- 冒烟曲线数据：`results/etfdr_fig1_smoke_curves.csv`
- 冒烟选择点：`results/etfdr_fig1_smoke_cv_cvc_selected.csv`
- 全量曲线数据：`results/etfdr_fig1_curves.csv`
- 全量选择点：`results/etfdr_fig1_cv_cvc_selected.csv`

## 备注

- 图 1 主图复现了“CV MSE / True FDR / FDR估计”三条曲线和 CV min / 1se 竖线。
- CVC vs CV 对照图在同类设置下比较了选点后的 `MSE` 与 `True FDR`。
