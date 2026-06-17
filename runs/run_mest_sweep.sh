#!/usr/bin/env bash
# M-estimator percentile sweep.
# One run at 5% WITH formula learners (Best_formula_GAM/GLM.txt) for comparison,
# then 1% / 5% / 10% / 20% WITHOUT formula learners (formula not tuned to this dataset).
# Runs sequentially — each job already uses mclapply internally.

set -euo pipefail
cd "$(dirname "$0")/.."   # always run from GRB_app root

INPUT="Data/superlearner_training_emcee_errcut_relative.csv"

run_sl() {
    local pct="$1"
    local use_formula="$2"
    local label="$3"
    local out_dir="runs/${label}/"
    mkdir -p "$out_dir"
    echo "============================================"
    echo "m_est_pct=${pct}  use_formula=${use_formula}  →  ${out_dir}"
    echo "============================================"
    Rscript superlearner.R \
        "$INPUT" TRUE FALSE TRUE FALSE 0.65 10 "$out_dir" "$pct" "$use_formula" \
        2>&1 | tee "${out_dir}run.log"
    echo "Done: ${out_dir}"
}

# 5% with formula learners (baseline comparison)
run_sl 0.05 TRUE  "mest_5pct_formula"

# 1% / 5% / 10% / 20% without formula learners
run_sl 0.01 FALSE "mest_1pct"
run_sl 0.05 FALSE "mest_5pct"
run_sl 0.10 FALSE "mest_10pct"
run_sl 0.20 FALSE "mest_20pct"

echo ""
echo "All 5 runs complete. Results in runs/mest_*/"
