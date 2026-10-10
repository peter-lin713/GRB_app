#!/usr/bin/env python3
"""
Random-projection / random-feature gates: can many random directions separate the
outlier manifold from the bulk when axis-aligned features can't?

Two published ideas, tested honestly:
  (1) Random features (Rahimi-Recht "random kitchen sinks"): RBF random Fourier
      features + logistic regression. This APPROXIMATES an RBF-SVM -- which we
      already saw scores ~0.44 -- so it is a consistency check.
  (2) Supervised random-projection selection (the "find projections that separate
      the two classes" idea): draw M random directions, keep the J that best
      separate the classes, classify on those. The crucial subtlety is WHERE the
      selection happens:
        HONEST : directions scored & selected INSIDE each training fold, then the
                 held-out fold is scored. Deployable.
        LEAKED : directions selected using ALL the data (incl. test labels) once,
                 then CV only the final classifier. This is the overfitting trap --
                 with 15 positives and thousands of directions you can always find
                 some that separate the training data by chance.

If HONEST stays ~0.5 while LEAKED shoots up, the "separation" is fitted noise, not
recoverable structure -- the same target-leakage pattern as the M-est cut.
Target: the 15 M-est outliers.
"""
import warnings, numpy as np
warnings.filterwarnings("ignore")
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.kernel_approximation import RBFSampler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from simple_linear_flagship_check import d, d_cut

FEAT = ["log10T90","log10Fa","log10Ta","Alpha","Beta","Gamma",
        "log10Fluence","PhotonIndex","log10NH","log10PeakFlux"]
Xr = d[FEAT].to_numpy()
X  = (Xr - Xr.mean(0)) / Xr.std(0)                        # standardized once (features only)
lab = d.index.isin(d.index.difference(d_cut.index)).astype(int)
CV = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)

# ---- (1) random Fourier features + logistic ---------------------------------
print("(1) Random features (RBF random Fourier features) + logistic regression:")
for D in (100, 500, 2000):
    pipe = make_pipeline(StandardScaler(), RBFSampler(n_components=D, gamma=0.1, random_state=0),
                         LogisticRegression(class_weight="balanced", max_iter=5000))
    s = cross_val_score(pipe, Xr, lab, cv=CV, scoring="roc_auc")
    print(f"    D={D:5d} random features:  AUC = {s.mean():.3f} ± {s.std():.3f}")

# ---- (2) supervised random-projection selection -----------------------------
class RandProjSelect(BaseEstimator, ClassifierMixin):
    """Draw M random unit directions, keep the J that best separate the classes
    (by |point-biserial corr| with the label), classify on those projections."""
    def __init__(self, M=3000, J=25, seed=0):
        self.M, self.J, self.seed = M, J, seed
    def fit(self, X, y):
        rng = np.random.default_rng(self.seed)
        W = rng.standard_normal((X.shape[1], self.M)); W /= np.linalg.norm(W, axis=0)
        P = X @ W                                        # (n, M) projections
        yc = (y - y.mean())
        score = np.abs((P - P.mean(0)).T @ yc) / (P.std(0) * np.linalg.norm(yc) + 1e-9)
        self.keep_ = np.argsort(score)[::-1][:self.J]
        self.W_ = W[:, self.keep_]
        self.clf_ = LogisticRegression(class_weight="balanced", max_iter=5000).fit(X @ self.W_, y)
        self.classes_ = self.clf_.classes_
        return self
    def predict_proba(self, X): return self.clf_.predict_proba(X @ self.W_)
    def decision_function(self, X): return self.clf_.decision_function(X @ self.W_)

print("\n(2) Supervised random-projection selection (M=3000 directions, keep J=25):")
honest = cross_val_score(RandProjSelect(), X, lab, cv=CV, scoring="roc_auc")
print(f"    HONEST  (select directions inside each fold):  AUC = {honest.mean():.3f} ± {honest.std():.3f}")
# LEAKED: select directions using ALL data, then CV only the final logistic
rng = np.random.default_rng(0)
W = rng.standard_normal((X.shape[1], 3000)); W /= np.linalg.norm(W, axis=0)
P = X @ W; yc = lab - lab.mean()
score = np.abs((P - P.mean(0)).T @ yc) / (P.std(0) * np.linalg.norm(yc) + 1e-9)
keep = np.argsort(score)[::-1][:25]
Xproj = X @ W[:, keep]
leaked = cross_val_score(LogisticRegression(class_weight="balanced", max_iter=5000),
                         Xproj, lab, cv=CV, scoring="roc_auc")
print(f"    LEAKED  (select directions on ALL data first):  AUC = {leaked.mean():.3f} ± {leaked.std():.3f}")
print("\nreference ceiling (route on |global residual|, needs z): AUC ~1.00")
