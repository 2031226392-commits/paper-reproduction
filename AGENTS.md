# AGENTS.md

## Project Goal

This repository reproduces the results of a statistical learning paper.

Target outputs:
- reproduced tables
- reproduced figures
- reproduced model metrics

---

## Environment

R version:
R 4.3 or later

Package management:
Use renv if available.

Setup:
install.packages("renv")
renv::restore()

If renv is not available, install required packages manually.

---

## Project Structure

R/          analysis and modeling code
scripts/    runnable pipeline scripts
data/       input datasets
results/    output tables, metrics, saved model objects
figures/    generated plots and figures

---

## Key Commands

Run full pipeline:

Rscript scripts/run_all.R

Fit model:

Rscript R/fit_model.R

Evaluate model:

Rscript R/evaluate_model.R

Generate figures:

Rscript R/make_figures.R

---

## Rules for Agents

- Prefer minimal and explainable changes
- Do not rewrite the entire analysis pipeline unless necessary
- Keep all assumptions in REPRO_GAP.md
- Save generated tables in results/
- Save generated figures in figures/
- If dataset is missing, stop and report clearly

---

## Expected Outputs

The final outputs should include:
- reproduced result tables
- reproduced figures
- a summary report in README_repro.md
- documented gaps in REPRO_GAP.md
