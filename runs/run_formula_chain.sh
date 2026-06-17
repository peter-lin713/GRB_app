#!/usr/bin/env bash
# Post-sweep chain:
#   1. Wait for the M-estimator sweep to finish
#   2. Run GAM formula search on the emcee dataset
#   3. Install fresh formulas as Best_formula_GAM.txt + Best_formula_GLM.txt
#   4. Run a final SuperLearner with use_formula=TRUE at 5% M-estimator

set -euo pipefail
cd "$(dirname "$0")/.."

INPUT="Data/superlearner_training_emcee_errcut_relative.csv"
SENTINEL="runs/mest_20pct/run.log"

# ── Step 1: wait for sweep ────────────────────────────────────────────────────
echo "$(date): Waiting for sweep to complete (watching for ${SENTINEL})..."
until [ -f "$SENTINEL" ] && { grep -q "MC metrics" "$SENTINEL" 2>/dev/null || grep -q "Halted\|halted\|Error" "$SENTINEL" 2>/dev/null; }; do
  sleep 120
done
echo "$(date): Sweep complete (or terminated)."

# ── Step 2: GAM formula search ───────────────────────────────────────────────
# correlation_cutoff = top 1% by CV correlation
# RMSE_cutoff        = bottom 3% by CV RMSE
echo "$(date): Running formula search..."
mkdir -p runs/formula_search
Rscript Formula_Search_Aditya.R \
    "$INPUT" \
    "runs/formula_search/" \
    0.99 \
    0.03 \
    2>&1 | tee runs/formula_search/formula_search.log
echo "$(date): Formula search done."

# ── Step 3: install new formulas ─────────────────────────────────────────────
# Formula_Search_Aditya.R writes Best_formula.txt to the working dir.
# Copy it to both files that superlearner.R reads.
if [ ! -f Best_formula.txt ]; then
    echo "ERROR: Best_formula.txt not found — formula search may have failed." >&2
    exit 1
fi
cp Best_formula.txt Best_formula_GAM.txt
cp Best_formula.txt Best_formula_GLM.txt
echo "$(date): Installed Best_formula_GAM.txt and Best_formula_GLM.txt."

# ── Step 4: SuperLearner with fresh formulas ──────────────────────────────────
OUT="runs/mest_5pct_new_formula/"
mkdir -p "$OUT"
echo "$(date): Running SuperLearner with new formulas (5% M-estimator)..."
Rscript superlearner.R \
    "$INPUT" TRUE FALSE TRUE FALSE 0.65 10 "$OUT" 0.05 TRUE \
    2>&1 | tee "${OUT}run.log"
echo "$(date): All done. Results in ${OUT}"
