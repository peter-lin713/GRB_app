#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"

EXPS="exp01_xray_only exp08_ot_baseline exp11_iterative_mice exp12_photonindex_relation \
      exp13_quantile_mapping exp14_kmeans_cluster exp15_gmm_softcluster exp16_copula \
      exp17_mlp_multitask exp18_svd_matrix_completion exp19_regression_residual_hybrid \
      exp20_stacked_ensemble"

echo "=== Generating datasets ==="
for exp in $EXPS; do
  echo "--- $exp ---"
  python3 "${exp}.py" 2>&1 | tail -6
done

echo
echo "=== Robustness check: 10 random seeds each ==="
python3 - << 'EOF'
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

FEATURES = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
            'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux']
RESPONSE = 'Redshift_crosscheck'

exps = "exp01_xray_only exp08_ot_baseline exp11_iterative_mice exp12_photonindex_relation \
      exp13_quantile_mapping exp14_kmeans_cluster exp15_gmm_softcluster exp16_copula \
      exp17_mlp_multitask exp18_svd_matrix_completion exp19_regression_residual_hybrid \
      exp20_stacked_ensemble".split()

results = {e: [] for e in exps}
ns = {}
for e in exps:
    df = pd.read_csv(f"results/{e}.csv")
    df = df[df[RESPONSE].notna()].copy()
    ns[e] = len(df)

for seed in range(10):
    for e in exps:
        df = pd.read_csv(f"results/{e}.csv")
        df = df[df[RESPONSE].notna()].copy()
        X = df[FEATURES].apply(pd.to_numeric, errors='coerce').values
        y = np.log10(pd.to_numeric(df[RESPONSE], errors='coerce').values + 1)
        model = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), RidgeCV(alphas=np.logspace(-3,3,25)))
        kf = KFold(n_splits=10, shuffle=True, random_state=seed)
        pred = cross_val_predict(model, X, y, cv=kf)
        r = np.corrcoef(y, pred)[0,1]
        results[e].append(r)

print(f"{'experiment':<32}{'n':>5}{'mean_r':>8}{'std':>8}{'min':>8}{'max':>8}")
summary = []
for e in exps:
    arr = np.array(results[e])
    summary.append((e, ns[e], arr.mean(), arr.std(), arr.min(), arr.max()))
summary.sort(key=lambda x: -x[2])
for e, n, m, s, lo, hi in summary:
    print(f"{e:<32}{n:>5}{m:>8.3f}{s:>8.3f}{lo:>8.3f}{hi:>8.3f}")

pd.DataFrame(summary, columns=['experiment','n','mean_r','std','min','max']).to_csv('results/round2_summary.csv', index=False)
EOF
