#!/usr/bin/env bash
# Outlier-removal experiment sweep.
# All runs use: emcee errcut dataset, MICE, 5% default M-est, no formula learners.
# outlier_method (arg 11): hard | soft | infold | both

set -uo pipefail
cd /Users/peterlin/Desktop/test/GRB_app

INPUT="Data/superlearner_training_emcee_errcut_relative.csv"
SL="outlier_experiments/superlearner.R"

run_sl() {
  local m_est_pct=$1
  local outlier_method=$2
  local tag=$3
  local out="outlier_experiments/runs/${tag}"
  mkdir -p "$out"
  echo "$(date): Starting ${tag} (m_est_pct=${m_est_pct}, method=${outlier_method})..."
  Rscript "$SL" \
    "$INPUT" TRUE FALSE TRUE FALSE 0.65 10 "$out" "$m_est_pct" FALSE "$outlier_method" \
    2>&1 | tee "${out}/run.log"
  echo "$(date): Done ${tag}."
}

# Option A: soft obsWeights from M-est rlm, no hard cut (method=soft)
run_sl 0.05 soft   "optA_soft_weights"

# Option B: in-fold rlm wrapper, no pre-cut (method=infold, do_m_estimator still TRUE but skipped for non-hard)
run_sl 0.05 infold "optB_infold_rlm"

# Combined 5%: hard 5% pre-cut + in-fold rlm wrapper on top (method=both)
run_sl 0.05 both   "combined_5pct"

# Combined 2%: hard 2% pre-cut + in-fold rlm wrapper (method=both)
run_sl 0.02 both   "combined_2pct"

echo "$(date): All outlier experiments complete."
