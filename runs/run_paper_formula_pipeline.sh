#!/usr/bin/env bash
set -euo pipefail
cd /Users/peterlin/Desktop/test/GRB_app

INPUT="Data/superlearner_training_emcee_errcut_relative.csv"
FORMULA_OUT="runs/formula_search_paper_cutoffs"
SL_OUT="runs/mest_5pct_paper_formula"

echo "$(date): Step 1 — formula search (0.999 / 0.02)..."
Rscript Formula_Search_Aditya.R "$INPUT" "$FORMULA_OUT/" 0.999 0.02 \
    2>&1 | tee "$FORMULA_OUT/formula_search.log"
echo "$(date): Formula search done."

if [ ! -f Best_formula.txt ]; then
    echo "ERROR: Best_formula.txt not produced." >&2; exit 1
fi
cp Best_formula.txt Best_formula_GAM.txt
cp Best_formula.txt Best_formula_GLM.txt
echo "$(date): Installed Best_formula_GAM.txt and Best_formula_GLM.txt."

echo "$(date): Step 2 — SuperLearner (5% M-estimator, use_formula=TRUE)..."
Rscript superlearner.R "$INPUT" TRUE FALSE TRUE FALSE 0.65 10 "$SL_OUT/" 0.05 TRUE \
    2>&1 | tee "$SL_OUT/run.log"
echo "$(date): All done. Results in $SL_OUT"
