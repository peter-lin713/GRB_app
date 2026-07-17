#!/usr/bin/env python3
"""exp02: naive baseline -- fill each optical-only GRB's 5 X-ray-only
features with the X-ray population's own median. No matching at all;
the floor any "smart" imputation method should beat."""
import os
from common import load_data, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

medians = {col: xray[col].median() for col in TARGET_COLS}
imputed = {col: [medians[col]] * len(opt_only_proj) for col in TARGET_COLS}
imputed_std = {col: [xray[col].std()] * len(opt_only_proj) for col in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp02_median_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp02_median_impute: {len(out)} GRBs -> {out_path}')
