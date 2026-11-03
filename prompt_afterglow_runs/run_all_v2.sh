#!/usr/bin/env bash
# Second series (runs 10-17, corrected GBM crosswalk; see build_gbm_v2_inputs.py).
# Run the prompt/afterglow pipeline runs strictly in series (one at a time,
# all cores to each). Per run:
#   1. formula search + M-estimator cut   GAM_supercomputer/pa_NN_<name>_formula_generation/
#   2. install that run's winning formula as Best_formula_GAM.txt / _GLM.txt
#   3. SuperLearner, flagship settings    runs/pa_NN_<name>/
#      (superlearner_CHECKPOINT_pre_daume_fix.R variant; MICE on, no upsampling,
#       no in-SL M-estimator, loop=10, formula learners on)
# Resumable: finished formula searches / SuperLearner runs are skipped, and the
# formula search itself resumes from its per-split checkpoints.
# The repo-root Best_formula_*.txt are backed up once and restored on exit.
#
# Usage (from repo root):  nohup caffeinate -i bash prompt_afterglow_runs/run_all_v2.sh &
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"
PA="$REPO/prompt_afterglow_runs"
LOG="$PA/progress_v2.log"
BK="$PA/best_formula_backup"

RUNS=(
  "10 prompt_nh_swiftpf prompt_nh"
  "11 prompt_nh_gbmpf2 prompt_nh"
  "12 prompt_gbmpf2 prompt"
  "13 all_gbmpf2 all"
  "14 all_gbmpf2cal all"
  "15 revmulti_all_gbmpf2 all"
  "16 all_gbmpf2_epeak0 all_epeak"
  "17 all_gbmpf2_epeakmice all_epeak"
)

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $*" | tee -a "$LOG"; }

mkdir -p "$BK"
[ -f "$BK/Best_formula_GAM.txt" ] || cp Best_formula_GAM.txt "$BK/"
[ -f "$BK/Best_formula_GLM.txt" ] || cp Best_formula_GLM.txt "$BK/"
restore() { cp "$BK/Best_formula_GAM.txt" "$BK/Best_formula_GLM.txt" "$REPO/"; log "Restored original Best_formula_GAM/GLM.txt"; }
trap restore EXIT

log "=== run_all.sh started (pid $$) ==="
for spec in "${RUNS[@]}"; do
  read -r NN NAME SET <<< "$spec"
  GDIR="GAM_supercomputer/pa_${NN}_${NAME}_formula_generation"
  ODIR="runs/pa_${NN}_${NAME}"

  # ---- 1. formula search + M-estimator cut
  if [ -f "$GDIR/final_outliers_removed.csv" ] && [ -f "$GDIR/Formula_for_outlier.txt" ]; then
    log "[$NN $NAME] formula search already complete -- skipping"
  else
    log "[$NN $NAME] formula search START"
    (cd "$GDIR" && Rscript GAM_analysis_8variables.R >> run_full_search.log 2>&1)
    if [ ! -f "$GDIR/final_outliers_removed.csv" ] || [ ! -f "$GDIR/Formula_for_outlier.txt" ]; then
      log "[$NN $NAME] formula search FAILED or incomplete -- see $GDIR/run_full_search.log; stopping"
      exit 1
    fi
    log "[$NN $NAME] formula search DONE | $(grep -E 'M-estimator outlier cut|Most-frequent winning' "$GDIR/run_full_search.log" | tail -2 | tr '\n' ' ')"
  fi

  # ---- 2. install the winning formula (same 3-line format as the originals)
  Rscript -e 'f <- paste(trimws(readLines(commandArgs(TRUE)[1])), collapse = " ");
              l <- rep(sprintf("[1] \"%s\"", f), 3);
              writeLines(l, "Best_formula_GAM.txt"); writeLines(l, "Best_formula_GLM.txt")' \
          "$GDIR/Formula_for_outlier.txt"
  log "[$NN $NAME] installed formula: $(head -1 Best_formula_GAM.txt)"

  # ---- 3. SuperLearner
  if [ -f "$ODIR/superlearner_model" ]; then
    log "[$NN $NAME] SuperLearner already complete -- skipping"
  else
    mkdir -p "$ODIR"
    cp "$GDIR/Formula_for_outlier.txt" "$ODIR/formula_used.txt"
    log "[$NN $NAME] SuperLearner START"
    Rscript "$PA/sl_scripts/superlearner_${SET}.R" "$GDIR/final_outliers_removed.csv" \
            TRUE FALSE FALSE FALSE 0.65 10 "$ODIR/" 0.05 TRUE > "$ODIR/run.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$ODIR/superlearner_model" ]; then
      log "[$NN $NAME] SuperLearner FAILED (exit $rc) -- see $ODIR/run.log; stopping"
      exit 1
    fi
    log "[$NN $NAME] SuperLearner DONE"
  fi
done
log "=== all v2 runs complete ==="
