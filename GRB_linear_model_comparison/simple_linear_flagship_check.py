#!/usr/bin/env python3
"""
=============================================================================
Does the flagship evaluation protocol inflate a SIMPLE linear model too?
-----------------------------------------------------------------------------
The flagship GRB redshift estimator reports r ~ 0.71 (linear z) using a
three-step evaluation protocol:
  (1) an M-estimator "catastrophic-outlier" cut, driven by the winning GLM
      formula, that drops the lowest-weight ~5% of bursts;
  (2) repeated 10-fold cross-validation of the model on what remains;
  (3) a post-hoc "2-sigma cone" clip that removes CV points lying > 2 sigma
      from the 1:1 line, then recomputes the correlation.

Question: is that ~0.71 something only the elaborate winning GLM achieves, or
does a PLAIN 7-feature additive linear model -- no interactions, no squared
terms -- reach nearly the same number under the IDENTICAL protocol?  If the
plain model also gets there, the inflation is the protocol removing generically
hard-to-predict bursts (plus the clip), not the specific model.

Reports, for the simple linear model (winning GLM shown for reference):
  honest_r      -- mean held-out r over 100 random 80/20 splits (NO cut, NO clip)
  flagship_cv   -- r from repeated 10-fold CV on the M-est-cut training set,
                   AVERAGED over 20 random 80/20 train/holdout splits
  flagship_clip -- the same CV predictions after the 2-sigma clip (the flagship
                   "without catastrophic outliers" number), also averaged

All correlations are Pearson r in LINEAR redshift, z = 10**prediction - 1.

Why average over splits: the flagship uses a single 80/20 split, but that one
split is arbitrary and swings the number by ~0.05-0.08 (which is why the exact
figure differs between R and Python -- different RNGs pick different splits).
Averaging over splits removes that dependence and is language-independent.
=============================================================================
"""
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

DATA = "grb_imputed_281.csv"   # 281 long GRBs, sits next to this script

# ---- Data -------------------------------------------------------------------
d = pd.read_csv(DATA, index_col=0)
for b in ["log10Fa", "log10T90", "log10PeakFlux", "log10NH", "PhotonIndex", "log10Ta", "Alpha"]:
    d[b + "Sqr"] = d[b] ** 2                              # squared terms, as the pipeline builds them
if "log10z" not in d.columns:
    d["log10z"] = np.log10(d["Redshift_crosscheck"] + 1)  # modelling target: log10(z+1)

# ---- Models -----------------------------------------------------------------
# Winning GLM (Eq. "winner"): defines the M-est cut AND is the flagship's
# dominant learner. In patsy, (a + b + ...)**2 = main effects + all pairwise
# interactions, exactly like R's (...)^2.
WINNER = ("log10z ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)**2 "
          "+ log10NH + PhotonIndex + log10Ta + Alpha "
          "+ log10NHSqr + log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr")

# SIMPLE model under test: the same 7 features, purely additive -- no
# interactions, no squares.
SIMPLE7 = "log10z ~ log10T90 + log10Fa + log10Ta + Alpha + log10NH + PhotonIndex + log10PeakFlux"

# ---- Step 1: the flagship M-estimator cut (always winning-formula-driven) ----
# Fit the winning formula with a Huber M-estimator; every burst gets an IRLS
# weight in (0, 1]. The flagship drops the lowest 5% by weight. This same cut is
# used for EVERY model we then evaluate -- exactly the flagship setup.
_y, _X = patsy.dmatrices(WINNER, d, return_type="dataframe")
_w = np.asarray(sm.RLM(_y, _X, M=sm.robust.norms.HuberT()).fit().weights)   # HuberT default c = 1.345
d_cut = d[_w > np.quantile(_w, 0.05)].copy()                                # ~266 bursts kept

# ---- Fast OLS helpers -------------------------------------------------------
# Building a patsy design matrix re-parses the formula, so do it ONCE per
# (formula, frame); then fit/predict on ROW SUBSETS with plain least squares.
# (The formulas use only row-wise terms, so subsetting rows of the full design
# is identical to rebuilding it on the subset.)
def design(formula, frame):
    y, X = patsy.dmatrices(formula, frame, return_type="matrix")
    return np.asarray(y).ravel(), np.asarray(X)

def r_lin(pred_log, obs_log):
    """Pearson r in LINEAR redshift, z = 10**x - 1."""
    return np.corrcoef(10 ** np.asarray(pred_log) - 1, 10 ** np.asarray(obs_log) - 1)[0, 1]

def _ols(Xtr, ytr, Xte):
    beta, *_ = np.linalg.lstsq(Xtr, ytr, rcond=None)
    return Xte @ beta

def _cv_clip(y, X, rng, reps, K):
    """Repeated K-fold CV on rows of (y, X) -> mean out-of-fold prediction ->
    2-sigma clip. Returns (CV r, CV+clip r) in linear z."""
    n = len(y); P = np.full((n, reps), np.nan)
    for _ in range(reps):
        fold = rng.permutation(np.arange(n) % K)
        for k in range(K):
            te = fold == k
            P[te, _] = _ols(X[~te], y[~te], X[te])
    cvp = P.mean(axis=1)
    L = np.std(y - cvp, ddof=1)                            # residual sigma (n-1, matches R's sd)
    keep = np.abs(cvp - y) < 2 * L                        # 2-sigma cone
    return r_lin(cvp, y), r_lin(cvp[keep], y[keep])

# ---- Metrics ----------------------------------------------------------------
def honest_r(formula, n_splits=100, seed=1):
    """Mean held-out r over random 80/20 splits on the FULL sample: refit per
    split, score the untouched 20%. No M-est cut, no clip."""
    y, X = design(formula, d)
    rng, n, rs = np.random.default_rng(seed), len(y), []
    for _ in range(n_splits):
        ho = rng.choice(n, size=int(0.20 * n), replace=False)
        te = np.zeros(n, bool); te[ho] = True
        rs.append(r_lin(_ols(X[~te], y[~te], X[te]), y[te]))
    return float(np.mean(rs))

def flagship_r(formula, n_splits=20, reps=5, K=10, seed=0):
    """Flagship protocol (M-est cut -> repeated K-fold CV -> 2-sigma clip)
    averaged over `n_splits` random 80/20 train/holdout splits of the cut set."""
    y, X = design(formula, d_cut)
    rng, n, cvs, clips = np.random.default_rng(seed), len(y), [], []
    for _ in range(n_splits):
        ho = rng.choice(n, size=int(0.20 * n), replace=False)
        tr = np.setdiff1d(np.arange(n), ho)               # ~213 training bursts
        cv, clip = _cv_clip(y[tr], X[tr], rng, reps, K)
        cvs.append(cv); clips.append(clip)
    return dict(cv=float(np.mean(cvs)), cv_sd=float(np.std(cvs)),
                clip=float(np.mean(clips)), clip_sd=float(np.std(clips)),
                n_train=n - int(0.20 * n))

# ---- Report -----------------------------------------------------------------
if __name__ == "__main__":
    print(f"M-est cut: {len(d_cut)} bursts kept of {len(d)} "
          f"(dropped lowest 5% by Huber weight)\n")

    s = flagship_r(SIMPLE7)
    print("Simple linear model (7 features, additive):")
    print(f"  honest held-out r (100x 80/20)        = {honest_r(SIMPLE7):.3f}")
    print(f"  flagship CV r          (avg 20 splits)= {s['cv']:.3f} +/- {s['cv_sd']:.3f}  (n_train={s['n_train']})")
    print(f"  flagship CV+2sigma-clip (avg 20 split)= {s['clip']:.3f} +/- {s['clip_sd']:.3f}\n")

    w = flagship_r(WINNER)
    print("Winning GLM (reference):")
    print(f"  honest held-out r (100x 80/20)        = {honest_r(WINNER):.3f}")
    print(f"  flagship CV r          (avg 20 splits)= {w['cv']:.3f} +/- {w['cv_sd']:.3f}")
    print(f"  flagship CV+2sigma-clip (avg 20 split)= {w['clip']:.3f} +/- {w['clip_sd']:.3f}\n")

    # Reproducibility: re-run the whole averaged estimate with different master
    # seeds; the mean should barely move.
    print("Simple linear, flagship CV+clip r (avg 20 splits) across master seeds:")
    for sd in (0, 1, 42, 123, 2024):
        rr = flagship_r(SIMPLE7, seed=sd)
        print(f"  seed={sd:<5d}  CV+clip r = {rr['clip']:.3f} +/- {rr['clip_sd']:.3f}")
