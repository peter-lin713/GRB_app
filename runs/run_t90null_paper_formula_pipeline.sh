#!/usr/bin/env bash
set -euo pipefail
cd /Users/petelin/Desktop/test/GRB_app

INPUT="Data/superlearner_training_emcee_t90null.csv"
FORMULA_OUT="runs/formula_search_t90null_paper_cutoffs"
SL_OUT="runs/t90null_paper_formula_5pct"
BACKUP_DIR="/private/tmp/claude-501/-Users-petelin-Desktop-test-GRB-app/e05be0d8-9ae7-4510-8c50-2e6d5f9d1f7a/scratchpad"

echo "$(date): Step 1 -- formula search on T90-null data (0.999 / 0.02)..."
Rscript Formula_Search_Aditya.R "$INPUT" "$FORMULA_OUT/" 0.999 0.02 \
    2>&1 | tee "$FORMULA_OUT/formula_search.log"
echo "$(date): Formula search done."

if [ ! -f Best_formula.txt ]; then
    echo "ERROR: Best_formula.txt not produced." >&2; exit 1
fi
cp Best_formula.txt Best_formula_GAM.txt
cp Best_formula.txt Best_formula_GLM.txt
echo "$(date): Installed Best_formula_GAM.txt and Best_formula_GLM.txt (T90-null-derived)."

echo "$(date): Step 2 -- SuperLearner (5% M-estimator, use_formula=TRUE)..."
Rscript superlearner.R "$INPUT" TRUE FALSE TRUE FALSE 0.65 10 "$SL_OUT/" 0.05 TRUE \
    2>&1 | tee "$SL_OUT/run.log"
echo "$(date): All done. Results in $SL_OUT"

echo "$(date): Restoring paper-formula Best_formula_GAM.txt / Best_formula_GLM.txt..."
cp "$BACKUP_DIR/Best_formula_GAM_paper_backup.txt" Best_formula_GAM.txt
cp "$BACKUP_DIR/Best_formula_GLM_paper_backup.txt" Best_formula_GLM.txt
echo "$(date): Restored."
