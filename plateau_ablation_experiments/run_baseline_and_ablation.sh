#!/usr/bin/env bash
set -uo pipefail
cd /Users/petelin/Desktop/test/GRB_app

echo "$(date): [1/2] Baseline SuperLearner (no opt cols, no formula learners) starting"
mkdir -p plateau_ablation_experiments/runs/baseline_no_optcols
Rscript superlearner.R "plateau_ablation_experiments/baseline_no_optcols.csv" \
    TRUE FALSE FALSE FALSE 0.65 10 plateau_ablation_experiments/runs/baseline_no_optcols 0.05 FALSE hard \
    > plateau_ablation_experiments/runs/baseline_no_optcols/run.log 2>&1
echo "$(date): Baseline done (exit $?)"

echo "$(date): [2/2] Ablation SuperLearner (no opt cols, no plateau params, no formula learners) starting"
mkdir -p plateau_ablation_experiments/runs/ablation_no_plateau
Rscript plateau_ablation_experiments/superlearner_no_plateau.R "plateau_ablation_experiments/ablation_no_plateau.csv" \
    TRUE FALSE FALSE FALSE 0.65 10 plateau_ablation_experiments/runs/ablation_no_plateau 0.05 FALSE hard \
    > plateau_ablation_experiments/runs/ablation_no_plateau/run.log 2>&1
echo "$(date): Ablation done (exit $?)"

echo "$(date): ALL DONE."
