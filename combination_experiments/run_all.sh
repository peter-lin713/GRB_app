#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"

echo "=== Generating datasets ==="
for exp in exp01_xray_only exp02_median_impute exp03_knn_impute exp04_linreg_impute \
           exp05_ridge_impute exp06_rf_impute exp07_gp_impute exp08_ot_baseline \
           exp09_ot_plus_t90 exp10_ot_errweighted; do
  echo "--- $exp ---"
  python3 "${exp}.py" 2>&1 | tail -5
done

echo
echo "=== Quick ridge check on each ==="
echo "experiment,n,r_log,rmse_log,r_z" > results/summary.csv
for exp in exp01_xray_only exp02_median_impute exp03_knn_impute exp04_linreg_impute \
           exp05_ridge_impute exp06_rf_impute exp07_gp_impute exp08_ot_baseline \
           exp09_ot_plus_t90 exp10_ot_errweighted; do
  csv="results/${exp}.csv"
  if [ -f "$csv" ]; then
    out=$(python3 ../quick_ridge_check.py "$csv" 2>&1)
    n=$(echo "$out" | grep -oE "n = [0-9]+" | grep -oE "[0-9]+")
    r_log=$(echo "$out" | grep "log10(z+1)" | grep -oE "r = [0-9.-]+" | grep -oE "[0-9.-]+")
    rmse_log=$(echo "$out" | grep "log10(z+1)" | grep -oE "RMSE = [0-9.-]+" | grep -oE "[0-9.-]+$")
    r_z=$(echo "$out" | grep "z scale" | grep -oE "r = [0-9.-]+" | grep -oE "[0-9.-]+")
    echo "${exp},${n},${r_log},${rmse_log},${r_z}" >> results/summary.csv
    echo "$exp: n=$n r_log=$r_log r_z=$r_z"
  else
    echo "$exp: NO OUTPUT CSV"
  fi
done

echo
echo "=== Summary (sorted by r_log) ==="
python3 -c "
import pandas as pd
df = pd.read_csv('results/summary.csv')
df = df.sort_values('r_log', ascending=False)
print(df.to_string(index=False))
"
