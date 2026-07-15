# Experiment ledger — GRB redshift CV r
Protocol: 5x5 paired folds (seed 31415) unless noted. Metric r = Pearson(OOF pred, log10(z+1)).
Baseline r=0.510. Permutation null −0.03±0.07.

| # | Hypothesis | Result | Verdict |
|---|---|---|---|
| 1 | library/metalearner/weights screen (8 cfgs) | best +0.009 | noise |
| 1b | paper-inspired: formula learners, top7, is_opt, upsampling | formula lib +0.029 | WIN |
| 2 | combos + dilution test + honest bias correction | BC hurts; optical rows dilute | insight |
| 3 | formula-lib combos | pf_top7 0.547 (+0.038) | WIN |
| 4 | recover real optical-catalog features | 0.572 (+0.025 vs R3) | BIG WIN |
| 5 | MICE 3-completion averaging | 0.579 (+0.004) | small win |
| 6 | expanded library (smooth GAMs, GP, band-GLM) | 0.570, no gain; band-GLM high weight | ceiling |
| 7 | corrected pipeline + fresh M-est cut + 10-fold | **0.589**; r_opt 0.42->0.51 | BIG WIN (new frame+protocol) |
| 8 | colors/diffs, PCA/kPCA, kNN, BART, Cubist, semi-MICE | semi-MICE +0.011 (4/5 reps); bagged 0.598 on old frame; others ~0; per-fold PCA broke | WIN: semi-MICE |
| 7b | partition bagging (retroactive on R7) | 0.589 -> 0.593 | free win |
| 9 | lit round on corrected frame: daume, PLS, err-noise, C-mixup, rank-gauss | daume 0.599 (+0.012, 3/3; r_opt 0.515->0.589!); PLS/noise/cmixup/rankgauss negative | WIN: daume |
| 10 | HEADLINE: corrected + semi-MICE + daume + fresh cut + avg3 + bagged 5x10f + full lib | **0.600** (r_xray 0.600, r_opt 0.554, linear-z 0.601) | NEW BEST; semi+daume gains overlap |
| 11 | screen: coral, wopt, coreg, poolpca, isoscore, knnenc | RUNNING | — |
| 11 | screen: coral +0.007 (r_opt 0.577), iso +0.006; wopt/coreg/poolpca/knnenc dead | 2 advance to ablation | mixed |
| 13 | ablation: daume +0.021 alone; CORAL+daume cannibalize (+0.011); iso adds 0 | daume-only confirmed optimal | insight |
| 14 | formula regeneration (fold-honest bestglm/bestgam selection) | -0.000/-0.001: no gain | clean negative (stale formulas not a bottleneck) |
| 15 | in-fold cut retest; inner V=10 | +0.003 / +0.002, sub-noise | park for headline retest |
| 16 | untied-variance Daumé c-sweep + V=20 | flat (best c=0.25 +0.004); scaling meaningless for unpenalized winners | null (mechanism understood) |
| 17 | additive learners: ppr, scam(monotone NH), gamboost | ppr -0.007, scam -0.005, gboost +0.003 | null/negative (ensemble saturated) |
| 18 | projection-draw bagging (emcee posterior) | param-only +0.015 (r_opt +0.05); full-scatter -0.05 | BIG WIN (param); scatter draws poison |
| 19 | imputation engines + robust losses | **missForest +0.020 (3/3 reps)**; sbgcop -0.004; huber -0.006 | BIG WIN: missForest |
| 19b | quad-EN (post-expansion daume) + nonneg-ridge stacking (lit rev 4) | RUNNING | — |
| 20 plan | FINAL: missForest(xseeds) x projection draws + daume + fresh cut + 10f bagged | — | — |

Lit rev 4 rulings: skip correlation losses (linear winners make it moot) & ordinal decomposition;
try quad-EN + ridge stacking; TabPFN cloud API flagged (needs user-provided key).

| 19b | quad-EN -0.005; nonneg-ridge stacking +0.006 | ridge stacking adopted | mixed |
| 20 | final consolidation (all-291 CV) | KILLED mid-run (protocol pivot) | superseded |
| 23 | PIPELINE CONVENTION eval (seed-42 233/58), winner config | **CV-233 r=0.661; holdout-58 r=0.459** (old frame same split: 0.324) | MILESTONE |
| 24 | max-stack (+infold, c=0.25, gamboost, iso, V=10) | RUNNING | — |

Protocol note: earlier all-291 CV numbers ran ~0.60-0.62 and are NOT the pipeline convention;
0.661 CV-233 is the apples-to-apples figure vs paper's 0.646 (pure-xray) and old runs' 0.50-0.60.
Holdout n=58 => SE~0.11; single-split holdout is a sanity check, not a precision measurement.

Observation: last 6 rounds all <= +0.004 — converging on ceiling ~0.60-0.61.
Side-find (r18 prep): emcee calibration slopes are shallow (logFa 0.42, logTa 0.36, Beta 0.05)
— projected optical Beta is ~constant; quantifies the cross-band information loss for Peter's OT work.

Adaptive lit review #2 key rulings: untied Daumé (Finkel-Manning) top pick; V>=20 for n<500
(Phillips 2023); m=20-40 justified; skip MTGP/lme4/Caruana-GES/SAM; missForest low-expectation check.
| 12 | Epeak enrichment (Fermi-GBM cross-match, 131/291 matched) | partial corr given top7 ~ 0 -> no added info; KILLED at screen | dead end (important negative) |

Note: TabPFN infeasible locally (torch>=2.5 has no Intel-mac wheels).
Round 9 plan: port round-8 winners onto round-7 corrected frame/protocol; add RFF,
pooled-PCA(unlabeled), OOF-kNN target encoding, pseudo-labeling, cut-pct sweep
(in-fold cut variant to keep eval set fixed).

Backlog (adaptive): TabPFN(reticulate), pseudo-labeling, pooled-PCA(unlabeled), RFF+glmnet,
CCA multi-view, isolation-forest/LOF cleaning, C-mixup/error-based augmentation, SIMEX,
residual boosting, gated per-band models, OOF-kNN target encoding, ridge interaction expansion,
stability selection, alt cut pcts on corrected frame, winsorized target, bagged stacking,
multi-seed averaging, Tukey-loss boosting, CORAL alignment for optical rows.
