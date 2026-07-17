#!/usr/bin/env python3
"""exp17: multi-output neural net (small MLP) predicting all 5 target
features jointly from the 4 shared params, fit on the 223 X-ray GRBs. Unlike
exp04-07 (one independent model per target), a shared hidden representation
can exploit correlations between the targets themselves (e.g. Fluence and
PeakFlux are related quantities) rather than treating each in isolation."""
import os
import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

X_train_raw = xray[FEATURE_COLS].values
X_query_raw = opt_only_proj[FEATURE_COLS].values
Y_train_raw = xray[TARGET_COLS].values

imp_x = SimpleImputer(strategy='median').fit(X_train_raw)
imp_y = SimpleImputer(strategy='median').fit(Y_train_raw)
X_train = imp_x.transform(X_train_raw)
X_query = imp_x.transform(X_query_raw)
Y_train = imp_y.transform(Y_train_raw)

x_scaler = StandardScaler().fit(X_train)
y_scaler = StandardScaler().fit(Y_train)
X_train_z = x_scaler.transform(X_train)
X_query_z = x_scaler.transform(X_query)
Y_train_z = y_scaler.transform(Y_train)

mlp = MLPRegressor(hidden_layer_sizes=(16, 8), activation='tanh', alpha=1.0,
                    max_iter=3000, random_state=42, early_stopping=True, validation_fraction=0.15)
mlp.fit(X_train_z, Y_train_z)
pred_z = mlp.predict(X_query_z)
pred = y_scaler.inverse_transform(pred_z)

# 5-fold CV residual std per target, as the uncertainty estimate
kf = KFold(n_splits=5, shuffle=True, random_state=42)
resid_std = np.zeros(len(TARGET_COLS))
for tr_idx, te_idx in kf.split(X_train_z):
    m = MLPRegressor(hidden_layer_sizes=(16, 8), activation='tanh', alpha=1.0,
                      max_iter=3000, random_state=42, early_stopping=True, validation_fraction=0.15)
    m.fit(X_train_z[tr_idx], Y_train_z[tr_idx])
    p = y_scaler.inverse_transform(m.predict(X_train_z[te_idx]))
    resid_std += (Y_train[te_idx] - p).std(axis=0) / kf.get_n_splits()

imputed = {col: pred[:, k] for k, col in enumerate(TARGET_COLS)}
imputed_std = {col: [resid_std[k]] * len(opt_only_proj) for k, col in enumerate(TARGET_COLS)}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp17_mlp_multitask.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp17_mlp_multitask: {len(out)} GRBs -> {out_path}')
