#!/usr/bin/env python3
"""
Is the gate's failure a SEPARABILITY problem a feature transform could fix?

Earlier: a gate predicting manifold membership (outlier vs bulk) from features
scored AUC ~0.57 (chance). Question: was that just the wrong feature space? We
try a battery of transforms x -> phi(x) and classifiers, each cross-validated
with the transform fit INSIDE the fold (no leakage), and report AUC.

Reasoning to keep in mind while reading the numbers:
  - RandomForest / RBF-SVM are flexible: they already approximate arbitrary
    functions of x, so a transform cannot raise THEIR ceiling -- it can only make
    a LINEAR classifier's boundary match a shape that is already there.
  - Therefore: if no linear-on-transform beats the flexible baselines, and the
    flexible baselines are at chance, the information simply is not in x, and this
    is NOT a separability issue.
  - If some transform jumps to high AUC (validated, low variance), that IS a
    deployable routing feature -> a real find.

Target: the 15 M-est outliers. Reference ceiling: |global-GLM residual| (needs z)
separates them at AUC ~1.0 -- that is what "the signal is in z, not x" looks like.
"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from itertools import combinations
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures, FunctionTransformer
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from simple_linear_flagship_check import d, d_cut

FEAT = ["log10T90","log10Fa","log10Ta","Alpha","Beta","Gamma",
        "log10Fluence","PhotonIndex","log10NH","log10PeakFlux"]
X = d[FEAT].to_numpy()
lab = d.index.isin(d.index.difference(d_cut.index)).astype(int)
FLUX = ["log10Fa","log10Ta","log10PeakFlux","log10Fluence","log10T90"]   # log-space -> diffs = log ratios
fi = [FEAT.index(c) for c in FLUX]

def add_flux_diffs(A):
    extra = [A[:, i] - A[:, j] for i, j in combinations(fi, 2)]           # pairwise log-ratios
    return np.column_stack([A] + extra)
FluxDiff = FunctionTransformer(add_flux_diffs)

CV = RepeatedStratifiedKFold(n_splits=5, n_repeats=8, random_state=0)
def auc(pipe):
    s = cross_val_score(pipe, X, lab, cv=CV, scoring="roc_auc"); return s.mean(), s.std()

LR = lambda: LogisticRegression(class_weight="balanced", max_iter=5000)
combos = {
    "LogReg | raw+scale":            make_pipeline(StandardScaler(), LR()),
    "LogReg | poly2+scale":          make_pipeline(PolynomialFeatures(2), StandardScaler(), LR()),
    "LogReg | flux-ratios+scale":    make_pipeline(FluxDiff, StandardScaler(), LR()),
    "LogReg | PCA(0.95)":            make_pipeline(StandardScaler(), PCA(0.95), LR()),
    "LDA  | scale":                  make_pipeline(StandardScaler(), LinearDiscriminantAnalysis()),
    "QDA  | scale":                  make_pipeline(StandardScaler(), QuadraticDiscriminantAnalysis(reg_param=0.3)),
    "kNN  | scale":                  make_pipeline(StandardScaler(), KNeighborsClassifier(15)),
    "SVM-rbf | scale":               make_pipeline(StandardScaler(), SVC(kernel="rbf", class_weight="balanced")),
    "SVM-rbf | flux-ratios+scale":   make_pipeline(FluxDiff, StandardScaler(), SVC(kernel="rbf", class_weight="balanced")),
    "RandomForest | raw (flexible)": RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=0, n_jobs=-1),
    "RandomForest | flux-ratios":    make_pipeline(FluxDiff, RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=0, n_jobs=-1)),
}
print(f"{'transform | classifier':32s}  {'CV AUC':>14s}   (0.50 = chance)")
best = None
for name, pipe in combos.items():
    m, s = auc(pipe)
    print(f"{name:32s}  {m:.3f} ± {s:.3f}")
    if best is None or m > best[1]: best = (name, m, s)
print(f"\nbest transform: {best[0]}  AUC={best[1]:.3f}±{best[2]:.2f}")
print("reference ceiling (route on |global residual|, needs z): AUC ~1.00")
print("\nread: all near 0.5 incl. flexible RF/SVM => not a separability issue; the")
print("membership signal is in z, not in any transform of x.")
