# GAM Analysis — Supercomputer Handoff

Two datasets, ready to run. Pick one, `cd` into its folder, and run:

```bash
Rscript GAM_analysis_8variables.R
```

- `emcee_v2/` — dataset built from the optical→X-ray MCMC calibration method
- `ot_v3/` — dataset built from the optimal-transport (OT) imputation method

Both use real T90 (no nulling). Each input file has **296 GRBs — outliers are
NOT pre-removed**. That's intentional: this matches how the original paper
did it (Dainotti et al. 2025, Sec. 4.2.3–4.3) —

1. Formula search runs on the full, uncut sample.
2. The single best-performing formula is picked.
3. *That* formula is then used to fit a robust regression (M-estimator) on
   the full sample, and the worst ~5% of GRBs by fit quality get dropped.

So the outlier cut depends on which formula wins — it can't be done ahead of
time. The script does step 3 automatically at the end of the run and writes
out the final cut dataset. See each folder's README for details.

`GAM_analysis_8variables.R` in this top folder is just a spare template copy —
you don't need it unless you're setting up a new dataset from scratch.
