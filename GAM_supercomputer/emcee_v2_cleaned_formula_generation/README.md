## emcee v2 formula generation — human-confirmed outlier cut

This folder is staged, not yet run. Workflow:

1. `superlearner_training_emcee_v2_recovered.csv` — 296 GRBs, post-combination +
   MICE + real-value recovery (today's fixes: single canonical optical source,
   live emcee-draw projection, MCMC refit every run). Not yet cut.
2. `candidate_outliers_emcee_v2.csv` — GRBs flagged by the automatic
   non-negotiable numeric cuts (Alpha>3, Beta>2, log10T90>6, PhotonIndex<0).
   Currently just 1 GRB (131117A, PhotonIndex<0) — the "log(err)/log(value)>1"
   criterion was tested and dropped from this list; it's numerically unstable
   whenever a variable's linear value sits near 1 (log10 -> 0), which produced
   nonsense flags (up to literal Inf) for PeakFlux/Alpha/Beta/PhotonIndex here.
3. Review the 4D fundamental-plane plots (Alpha, Beta, log10Fa, log10Ta) for
   additional physically-implausible GRBs the numeric cuts miss.
4. Edit `confirmed_outliers_to_drop.txt` (one GRB per line) — pre-populated
   with the candidate above; add or remove GRBs based on the plot review.
5. Run `GAM_analysis_8variables.R` (`MyLaptop = T`, ~4.3 hr for the 100-split
   search on this machine, matches the emcee/OT runs already done this
   session). It will:
   - remove every GRB in `confirmed_outliers_to_drop.txt` completely (not
     NA+impute), and write the final drop list + reasons to
     `removed_outliers_final.csv`;
   - run the standard LASSO -> O1/O2/Smoothed formula sweep -> 100-split
     search -> most-frequent-winner M-estimator cut, exactly like
     `emcee_v2_recovered_formula_generation/`;
   - produce `Formula_for_outlier.txt` and `final_outliers_removed.csv`, ready
     to feed into the SuperLearner stage the same way the other
     `*_formula_generation_output/` folders have been used this session.
