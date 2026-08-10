# Formulas found this batch

| File | Dataset | How found | Result (r_z, Without 2σ) |
|---|---|---|---|
| `standing_emcee_r070_formula.txt` | single-variable emcee (original) | full 100-split search (earlier session) | 0.712 — reference/paper winner |
| `emcee_own_fresh_search_formula.txt` | single-variable emcee, re-searched with corrected MICE | full 100-split search | 0.639 (all) — worse than standing, not pursued further |
| `ransac_own_formula.txt` | RANSAC calibration | full 100-split search, corrected linear-z/bias-gate methodology | see main results table |
| `theilsen_own_formula.txt` | Theil-Sen calibration | full 100-split search, corrected methodology | 0.681 — slightly worse than reusing the standing formula (0.707) |
| `hardzero_domain_own_formula.txt` | hard-zero domain split **v1** (Alpha_opt/Beta_opt/etc, no regression fit) | full 100-split search; top 7 ranked formulas all hit exact M-estimator singularity (multiple `_opt` columns sharing one zero/nonzero pattern) -- this is the first one down the ranked list (8th, 2 votes) that fits `rlm()` cleanly | 0.683 |
| `multivariate_emcee_own_formula.txt` | multivariate emcee (MCMC) calibration, Dataset A (core cols only) | full 100-split search | 0.672 — worse than reusing the standing formula (0.705 with opt cols added, Dataset B) |
| `multivariate_theilsen_own_formula.txt` | multivariate Theil-Sen calibration (all 4 optical params -> each X-ray-scale plateau param, no opt cols) | full 100-split search | 0.692 — slightly worse than reusing the standing formula (0.707), same recurring pattern as single-variable Theil-Sen |
| `no_plateau_own_formula.txt` | no-plateau ablation (6-var: Alpha/Beta/log10Fa/log10Ta dropped entirely, no opt cols) | full 100-split search | 0.644 — largest accuracy cost of any scheme tested |
| `hardzero_domain_v2_own_formula.txt` | hard-zero domain split **v2** (shared GRBs' opt cols left NaN for MICE instead of hard-zeroed) | full 100-split search; M-estimator fallback triggered again (rank 1 singular, rank 2 fit cleanly) -- confirms the fallback logic built into the v2 script generalizes | 0.691 -- essentially ties v1 (0.683) |
| (none -- reuses `standing_emcee_r070_formula.txt`) | Theil-Sen multivariate + opt cols added (v2 convention) | no search -- deliberately reused the standing formula, since reusing it has consistently beaten fresh searches all batch | **0.708** -- best reused-formula result besides the two ~0.71 references |

`Best_formula_GAM.txt` at the repo root is restored to the standing r=0.70 formula by default; each queued chain installs its own dataset's formula only for the duration of that dataset's own SuperLearner run, then restores the standing one.
