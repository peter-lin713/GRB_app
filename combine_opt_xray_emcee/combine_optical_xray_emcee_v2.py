#!/opt/anaconda3/envs/grb/bin/python3
"""
Combine optical and x-ray GRB data using MCMC (emcee) calibration -- v2.

Same emcee calibration approach as combine_optical_xray_emcee.py, but reads
the newer pre-processed/MICE-imputed/error-cut catalogs (which already carry
unit-fixed Fluence/NH, corrected T90, and Gamma/Fluence/NH/PeakFlux on the
optical side too) instead of the original optical_data.csv/xray_data.csv.
Writes to separate output/plot/cache paths so the original outputs are left
untouched.

Input files (both CSV with named columns):
  OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv
      GRB, z, class, log10Faopt, log10FaErr_opt, log10Taopt, logTaErr_opt,
      Alpha_opt, AlphaErr_opt, Beta_opt, betaErr_opt, PhotonIndex,
      PhotonIndexErr, Gamma, log10Fluence, log10FluenceErr, log10NH,
      log10PeakFlux, log10PeakFluxErr, log10T90
  Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv
      GRB_Name, Redshift_crosscheck, log10Fluence, log10PeakFlux,
      PhotonIndex, log10NH, Gamma, log10Fa_X, log10Ta_X, Alpha_X, Beta_X,
      log10TaErr, BetaErr_X, AlphaErr, T90Err, PhotonIndexErr, log10FaErr_X,
      logT90err, log10PeakFluxErr, logFluenceErr, log10T90

Unlike the original catalogs, neither file carries raw (linear) T90 anymore
-- only log10T90 -- so T90 is combined directly in log space instead of via
np.log10(final['T90']). Both files now also carry Gamma/Fluence/NH/PeakFlux
for every GRB (not just x-ray), so those fields are combined the same
xray-preferred/optical-fallback way as logFa/logTa/Alpha/Beta rather than
being left NaN for optical-only GRBs.
"""

import os
import numpy as np
import pandas as pd
import emcee
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
REPO_DIR    = os.path.dirname(SCRIPT_DIR)
DATA_DIR    = os.path.join(REPO_DIR, 'Data')
OUT_DIR     = REPO_DIR
OPT_FILE    = os.path.join(DATA_DIR, 'OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv')
XRAY_FILE   = os.path.join(DATA_DIR, 'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
OUT_FILE    = os.path.join(DATA_DIR, 'superlearner_training_emcee_v2.csv')
PLOT_DIR    = os.path.join(OUT_DIR, 'Plot_Output', 'emcee_calibration_v2')
CHAINS_FILE = os.path.join(DATA_DIR, 'emcee_chains_v2.npz')   # cached chains; delete to rerun MCMC

N_WALKERS  = 32
N_STEPS    = 3000
N_BURN     = 500
N_THIN     = 15
MIN_ERR    = 1e-6

# Load optical data

opt = pd.read_csv(OPT_FILE)
opt_sel = opt[['GRB', 'z', 'log10T90',
               'log10Faopt', 'log10FaErr_opt', 'log10Taopt', 'logTaErr_opt',
               'Alpha_opt', 'AlphaErr_opt', 'Beta_opt', 'betaErr_opt',
               'Gamma', 'log10Fluence', 'log10FluenceErr',
               'PhotonIndex', 'PhotonIndexErr',
               'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()
opt_sel = opt_sel.rename(columns={
    'log10Faopt':   'logFa_opt',  'log10FaErr_opt': 'logFaErr_opt',
    'log10Taopt':   'logTa_opt',  'logTaErr_opt':   'logTaErr_opt',
    'AlphaErr_opt': 'AlphaErr_opt',
    'betaErr_opt':  'BetaErr_opt',
})
for col in opt_sel.columns:
    if col != 'GRB':
        opt_sel[col] = pd.to_numeric(opt_sel[col], errors='coerce')

print(f'Optical data loaded:  {opt_sel.shape[0]} GRBs')

# Load x-ray data

xray = pd.read_csv(XRAY_FILE)
xray = xray.rename(columns={
    'GRB_Name':            'GRB',
    'Redshift_crosscheck': 'z',
    'log10Fa_X':    'logFa_x',   'log10FaErr_X': 'logFaErr_x',
    'log10Ta_X':    'logTa_x',   'log10TaErr':   'logTaErr_x',
    'Alpha_X':      'Alpha_x',   'AlphaErr':     'AlphaErr_x',
    'Beta_X':       'Beta_x',    'BetaErr_X':    'BetaErr_x',
    'logFluenceErr': 'log10FluenceErr',
})
xray['GRB'] = xray['GRB'].str.replace('GRB', '', regex=False)

# Normalize: x-ray IDs with no letter suffix get an 'A' appended so they
# match the optical catalog's consistent use of letter suffixes.
xray['GRB_norm'] = xray['GRB'].apply(
    lambda g: g if (g and g[-1].isalpha()) else g + 'A'
)

print(f'X-ray data loaded:    {xray.shape[0]} GRBs')

xray_sel = xray[[
    'GRB_norm', 'z', 'log10T90', 'T90Err', 'Gamma',
    'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
    'log10NH', 'log10PeakFlux', 'log10PeakFluxErr',
    'logFa_x', 'logFaErr_x', 'logTa_x', 'logTaErr_x',
    'Alpha_x', 'AlphaErr_x', 'Beta_x', 'BetaErr_x',
]].rename(columns={'GRB_norm': 'GRB'}).copy()

# Overlapping GRBs

overlap = set(xray_sel['GRB']) & set(opt_sel['GRB'])
print(f'Overlapping GRBs:     {len(overlap)}  '
      f'(exact suffix matches + naming-convention fixes)')

merged = pd.merge(
    xray_sel[['GRB', 'logFa_x', 'logFaErr_x', 'logTa_x', 'logTaErr_x',
              'Alpha_x', 'AlphaErr_x', 'Beta_x', 'BetaErr_x']],
    opt_sel,
    on='GRB',
)
print(f'Merged overlap shape: {merged.shape}')

# MCMC linear calibration fits
# Model: y = m*x + b,  total variance V_i = sy_i^2 + m^2*sx_i^2 + exp(2*log_s)
# Parameters theta = [m, b, log_s]

def log_likelihood(theta, x, y, sx, sy):
    m, b, log_s = theta
    V = sy**2 + m**2 * sx**2 + np.exp(2 * log_s)
    return -0.5 * np.sum((y - m * x - b)**2 / V + np.log(2 * np.pi * V))

def log_prior(theta):
    m, b, log_s = theta
    if -10 < m < 10 and -50 < b < 50 and -10 < log_s < 5:
        return 0.0
    return -np.inf

def log_posterior(theta, x, y, sx, sy):
    lp = log_prior(theta)
    return lp + log_likelihood(theta, x, y, sx, sy) if np.isfinite(lp) else -np.inf

params_map = {
    'logFa': ('logFa_opt',  'logFaErr_opt',  'logFa_x',  'logFaErr_x'),
    'logTa': ('logTa_opt',  'logTaErr_opt',  'logTa_x',  'logTaErr_x'),
    'Alpha': ('Alpha_opt',  'AlphaErr_opt',  'Alpha_x',  'AlphaErr_x'),
    'Beta':  ('Beta_opt',   'BetaErr_opt',   'Beta_x',   'BetaErr_x'),
}

# Per-parameter outlier cuts applied before MCMC fitting.
# Tuple: (x_max, y_max, include_err)
#   include_err=False -> remove if value > max on either axis
#   include_err=True  -> remove if value + error > max on either axis
param_cuts = {
    'logTa': (6.0, 6.0, False),
    'Alpha': (2.5, 2.5, True),
    'Beta':  (1.5, 1.5, True),
}

# Remove calibration GRBs where err/|value| > this on either axis.
CALIB_ERR_RATIO = 0.5

def apply_calib_cuts(sub, xc, xe, yc, ye, param):
    if param in param_cuts:
        x_max, y_max, use_err = param_cuts[param]
        before = len(sub)
        if use_err:
            sub = sub[(sub[xc] + sub[xe] < x_max) & (sub[yc] + sub[ye] < y_max)]
            label = 'val+err'
        else:
            sub = sub[(sub[xc] < x_max) & (sub[yc] < y_max)]
            label = 'val'
        if len(sub) < before:
            print(f'  {param}: dropped {before - len(sub)} value-outlier(s) ({label} >= {x_max})')
    before = len(sub)
    sub = sub[
        (sub[xe].abs() / sub[xc].abs().clip(lower=1e-10) <= CALIB_ERR_RATIO) &
        (sub[ye].abs() / sub[yc].abs().clip(lower=1e-10) <= CALIB_ERR_RATIO)
    ]
    if len(sub) < before:
        print(f'  {param}: dropped {before - len(sub)} error-ratio outlier(s) (err/|val| > {CALIB_ERR_RATIO})')
    return sub

mcmc_results = {}
print()

if os.path.exists(CHAINS_FILE):
    print(f'Loading cached MCMC chains from {CHAINS_FILE}')
    cache = np.load(CHAINS_FILE, allow_pickle=False)
    for param in params_map:
        flat = cache[param]
        m_med, b_med, ls_med = np.median(flat, axis=0)
        m_std, b_std, _      = np.std(flat, axis=0)
        sigma_int             = np.exp(ls_med)
        sub = merged[[params_map[param][0], params_map[param][2]]].apply(
            pd.to_numeric, errors='coerce').dropna()
        mcmc_results[param] = dict(
            m=m_med, b=b_med, sigma_int=sigma_int,
            m_err=m_std, b_err=b_std, flat=flat, n=len(sub)
        )
        print(f'{param}: m={m_med:.4f}±{m_std:.4f}  b={b_med:.4f}±{b_std:.4f}  '
              f'σ_int={sigma_int:.4f}  N={len(sub)}  [cached]')
else:
    chains_to_save = {}
    for param, (xc, xe, yc, ye) in params_map.items():
        sub = merged[[xc, xe, yc, ye]].apply(pd.to_numeric, errors='coerce').dropna()
        sub = apply_calib_cuts(sub, xc, xe, yc, ye, param)
        x_vals = sub[xc].values
        y_vals = sub[yc].values
        sx     = np.maximum(sub[xe].values, MIN_ERR)
        sy     = np.maximum(sub[ye].values, MIN_ERR)

        # OLS starting point
        cov = np.cov(x_vals, y_vals)
        m0  = cov[0, 1] / np.var(x_vals)
        b0  = np.mean(y_vals) - m0 * np.mean(x_vals)
        p0  = np.array([m0, b0, np.log(0.1)])
        pos = p0 + 1e-3 * np.random.default_rng(42).standard_normal((N_WALKERS, 3))

        sampler = emcee.EnsembleSampler(
            N_WALKERS, 3, log_posterior, args=(x_vals, y_vals, sx, sy)
        )
        sampler.run_mcmc(pos, N_STEPS, progress=True)

        flat = sampler.get_chain(discard=N_BURN, thin=N_THIN, flat=True)
        m_med,  b_med,  ls_med  = np.median(flat, axis=0)
        m_std,  b_std,  ls_std  = np.std(flat, axis=0)
        sigma_int = np.exp(ls_med)

        mcmc_results[param] = dict(
            m=m_med, b=b_med, sigma_int=sigma_int,
            m_err=m_std, b_err=b_std, flat=flat, n=len(x_vals)
        )
        chains_to_save[param] = flat
        print(f'{param}: m={m_med:.4f}±{m_std:.4f}  b={b_med:.4f}±{b_std:.4f}  '
              f'σ_int={sigma_int:.4f}  N={len(x_vals)}')

    np.savez(CHAINS_FILE, **chains_to_save)
    print(f'\nMCMC chains cached to {CHAINS_FILE}')

# Calibration diagnostic plots: optical vs x-ray for overlapping GRBs

os.makedirs(PLOT_DIR, exist_ok=True)
axis_labels = {
    'logFa': (r'log$_{10}(F_a)$ optical', r'log$_{10}(F_a)$ X-ray'),
    'logTa': (r'log$_{10}(T_a)$ optical', r'log$_{10}(T_a)$ X-ray'),
    'Alpha': (r'$\alpha$ optical',         r'$\alpha$ X-ray'),
    'Beta':  (r'$\beta$ optical',          r'$\beta$ X-ray'),
}

fig, axes = plt.subplots(2, 2, figsize=(10, 9))
fig.suptitle('Emcee calibration v2: optical vs X-ray (overlapping GRBs)\n'
             'Line = posterior median fit; band = 1σ posterior uncertainty',
             fontsize=11)

for ax, (param, (xc, xe, yc, ye)) in zip(axes.flat, params_map.items()):
    sub = merged[[xc, xe, yc, ye]].apply(pd.to_numeric, errors='coerce').dropna()
    sub = apply_calib_cuts(sub, xc, xe, yc, ye, param)
    x_d  = sub[xc].values
    y_d  = sub[yc].values
    sx   = sub[xe].values
    sy   = sub[ye].values

    r    = mcmc_results[param]
    flat = r['flat']
    m_med, b_med = r['m'], r['b']

    x_line = np.linspace(x_d.min(), x_d.max(), 200)
    y_line = m_med * x_line + b_med

    # 1σ band: envelope over all posterior-sample lines
    y_band = flat[:, 0:1] * x_line[None, :] + flat[:, 1:2]   # (n_samp, 200)
    y_lo_band = np.percentile(y_band, 16, axis=0)
    y_hi_band = np.percentile(y_band, 84, axis=0)

    ax.errorbar(x_d, y_d, xerr=sx, yerr=sy,
                fmt='o', ms=3, color='steelblue', alpha=0.6,
                elinewidth=0.7, capsize=2, label=f'Overlap (N={len(x_d)})')
    ax.plot(x_line, x_line, 'k--', lw=1, alpha=0.5, label='1:1')
    ax.plot(x_line, y_line, 'r-', lw=1.5,
            label=f'm={m_med:.2f}, b={b_med:.2f}')
    ax.fill_between(x_line, y_lo_band, y_hi_band,
                    color='red', alpha=0.15, label='1σ posterior band')

    xlab, ylab = axis_labels[param]
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(param, fontsize=10)
    ax.legend(fontsize=7, loc='upper left')
    ax.tick_params(labelsize=8)

plt.tight_layout()
out_path = os.path.join(PLOT_DIR, 'emcee_optical_vs_xray.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f'\nCalibration plot saved to {out_path}')

# Transform optical-only GRBs to x-ray scale
# Uses full posterior predictive propagation: for every emcee sample (m_i, b_i, sigma_i)
# we also perturb each optical measurement within its 1-sigma error, then take the
# median as the point estimate and half the 16th-84th percentile range as the 1-sigma
# uncertainty. This avoids the Gaussian approximation of the analytical error formula.

opt_tr = opt_sel[['GRB', 'z', 'log10T90', 'Gamma',
                   'log10Fluence', 'log10FluenceErr',
                   'PhotonIndex', 'PhotonIndexErr',
                   'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()

rng = np.random.default_rng(42)

for param, (xc, xe, _, _) in params_map.items():
    r    = mcmc_results[param]
    flat = r['flat']                        # (n_samples, 3)
    m_samp  = flat[:, 0]                    # posterior slope samples
    b_samp  = flat[:, 1]                    # posterior intercept samples
    s_samp  = np.exp(flat[:, 2])            # posterior intrinsic-scatter samples

    x_vals = pd.to_numeric(opt_sel[xc], errors='coerce').values   # (n_grbs,)
    x_errs = pd.to_numeric(opt_sel[xe], errors='coerce').values
    x_errs = np.where(np.isfinite(x_errs), np.maximum(x_errs, MIN_ERR), MIN_ERR)

    n_samp, n_grb = len(flat), len(x_vals)

    # Perturb each optical measurement within its 1-sigma error, one draw per sample.
    x_pert = x_vals[None, :] + rng.standard_normal((n_samp, n_grb)) * x_errs[None, :]

    # Predicted x-ray value: linear transform + intrinsic scatter draw.
    preds = (m_samp[:, None] * x_pert + b_samp[:, None]
             + rng.standard_normal((n_samp, n_grb)) * s_samp[:, None])

    opt_tr[f'{param}_x']    = np.nanmedian(preds, axis=0)
    lo = np.nanpercentile(preds, 16, axis=0)
    hi = np.nanpercentile(preds, 84, axis=0)
    opt_tr[f'{param}Err_x'] = (hi - lo) / 2.0      # 1-sigma half-width

    print(f'{param}: median 1σ err = {np.nanmedian(opt_tr[f"{param}Err_x"]):.4f}'
          f'  (σ_int={r["sigma_int"]:.4f})')

opt_tr = opt_tr.rename(columns={
    'logFa_x': 'logFa_x_opt', 'logFaErr_x': 'logFaErr_x_opt',
    'logTa_x': 'logTa_x_opt', 'logTaErr_x': 'logTaErr_x_opt',
    'Alpha_x': 'Alpha_x_opt', 'AlphaErr_x': 'AlphaErr_x_opt',
    'Beta_x':  'Beta_x_opt',  'BetaErr_x':  'BetaErr_x_opt',
})

# Outer merge; x-ray values take priority for overlapping GRBs

combined = pd.merge(
    opt_tr, xray_sel,
    on='GRB', how='outer', suffixes=('_from_opt', '_from_xray')
)

def resolve(df, xray_col, opt_col):
    xc = df.get(xray_col)
    oc = df.get(opt_col)
    if xc is None:
        return pd.to_numeric(oc, errors='coerce') if oc is not None else np.nan
    if oc is None:
        return pd.to_numeric(xc, errors='coerce')
    return pd.to_numeric(xc, errors='coerce').combine_first(pd.to_numeric(oc, errors='coerce'))

final = pd.DataFrame()
final['GRB'] = combined['GRB']

z_x = combined.get('z_from_xray', combined.get('z'))
z_o = combined.get('z_from_opt',  combined.get('z'))
final['Redshift_crosscheck'] = (
    pd.to_numeric(z_x, errors='coerce')
    .combine_first(pd.to_numeric(z_o, errors='coerce'))
)

# Neither catalog carries raw (linear) T90 anymore -- combine directly in log space.
final['log10T90'] = resolve(combined, 'log10T90_from_xray', 'log10T90_from_opt')

final['log10Fa']    = resolve(combined, 'logFa_x',    'logFa_x_opt')
final['log10FaErr'] = resolve(combined, 'logFaErr_x', 'logFaErr_x_opt')
final['log10Ta']    = resolve(combined, 'logTa_x',    'logTa_x_opt')
final['log10TaErr'] = resolve(combined, 'logTaErr_x', 'logTaErr_x_opt')
final['Alpha']      = resolve(combined, 'Alpha_x',    'Alpha_x_opt')
final['AlphaErr']   = resolve(combined, 'AlphaErr_x', 'AlphaErr_x_opt')
final['Beta']       = resolve(combined, 'Beta_x',     'Beta_x_opt')
final['BetaErr']    = resolve(combined, 'BetaErr_x',  'BetaErr_x_opt')

# Both catalogs now carry Gamma/Fluence/PhotonIndex/NH/PeakFlux (the old script only had
# these from the x-ray side) -- combine x-ray-preferred/optical-fallback like the params above
# instead of leaving optical-only GRBs NaN here.
for col in ['Gamma', 'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
            'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']:
    final[col] = resolve(combined, f'{col}_from_xray', f'{col}_from_opt')

final['T90Err'] = combined.get('T90Err', np.nan)   # only ever tracked on the x-ray side

# Finalize & save

col_order = [
    'GRB', 'Redshift_crosscheck', 'log10T90',
    'log10Fa', 'log10Ta', 'Alpha', 'Beta',
    'Gamma', 'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux',
    'T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
    'log10FluenceErr', 'PhotonIndexErr', 'log10PeakFluxErr',
]

out = final[col_order].copy()
out = out[out['Redshift_crosscheck'].notna()].copy()
out = out.set_index('GRB')

# Identify optical-only GRBs (projected, not in x-ray catalog)
xray_ids = set(xray_sel['GRB'])
opt_only_mask = ~out.index.isin(xray_ids)

# ---- Ratio error filter: err/|value| > 0.5 on original optical AND projected x-ray ----
# Original optical: avoids σ_int flattening all projected errors to the same value.
# Projected x-ray:  catches GRBs whose calibrated estimate itself is unreliable.
# A GRB is removed if ANY param fails the ratio cut on EITHER side.

RATIO_THRESH = 0.5

opt_only_grbs = out.index[opt_only_mask]
opt_orig = opt_sel[opt_sel['GRB'].isin(opt_only_grbs)].set_index('GRB')
opt_orig = opt_orig.apply(pd.to_numeric, errors='coerce')

# Original optical side
fa_cut_orig    = ((opt_orig['logFa_opt'].abs() > 0) & (opt_orig['logFaErr_opt'] / opt_orig['logFa_opt'].abs() > RATIO_THRESH)).reindex(opt_only_grbs, fill_value=False)
ta_cut_orig    = ((opt_orig['logTa_opt'].abs() > 0) & (opt_orig['logTaErr_opt'] / opt_orig['logTa_opt'].abs() > RATIO_THRESH)).reindex(opt_only_grbs, fill_value=False)
alpha_cut_orig = ((opt_orig['Alpha_opt'].abs() > 0) & (opt_orig['AlphaErr_opt'] / opt_orig['Alpha_opt'].abs() > RATIO_THRESH)).reindex(opt_only_grbs, fill_value=False)
beta_cut_orig  = ((opt_orig['Beta_opt'].abs()  > 0) & (opt_orig['BetaErr_opt']  / opt_orig['Beta_opt'].abs()  > RATIO_THRESH)).reindex(opt_only_grbs, fill_value=False)

# Projected x-ray side (from `out`)
out_opt = out.loc[opt_only_grbs].apply(pd.to_numeric, errors='coerce')
fa_cut_proj    = (out_opt['log10Fa'].abs() > 0) & (out_opt['log10FaErr'] / out_opt['log10Fa'].abs() > RATIO_THRESH)
ta_cut_proj    = (out_opt['log10Ta'].abs() > 0) & (out_opt['log10TaErr'] / out_opt['log10Ta'].abs() > RATIO_THRESH)
alpha_cut_proj = (out_opt['Alpha'].abs()   > 0) & (out_opt['AlphaErr']   / out_opt['Alpha'].abs()   > RATIO_THRESH)
beta_cut_proj  = (out_opt['Beta'].abs()    > 0) & (out_opt['BetaErr']    / out_opt['Beta'].abs()    > RATIO_THRESH)

# Combined: flag if EITHER side exceeds the threshold
fa_cut    = fa_cut_orig    | fa_cut_proj
ta_cut    = ta_cut_orig    | ta_cut_proj
alpha_cut = alpha_cut_orig | alpha_cut_proj
beta_cut  = beta_cut_orig  | beta_cut_proj

any_cut  = fa_cut | ta_cut | alpha_cut | beta_cut
drop_ids = any_cut[any_cut].index
out_cut  = out.drop(index=drop_ids)

cut_file = OUT_FILE.replace('.csv', '_errcut_relative.csv')
out_cut.to_csv(cut_file, index=True)
print(f'\nRatio cut (err/|value| > {RATIO_THRESH}) on original optical + projected x-ray:')
print(f'  Removed {len(drop_ids)} optical-only GRBs ({out_cut.shape[0]} remain)  → {cut_file}')
print(f'  logFa:  opt {fa_cut_orig.sum()}  proj {fa_cut_proj.sum()}  combined {fa_cut.sum()}')
print(f'  logTa:  opt {ta_cut_orig.sum()}  proj {ta_cut_proj.sum()}  combined {ta_cut.sum()}')
print(f'  Alpha:  opt {alpha_cut_orig.sum()}  proj {alpha_cut_proj.sum()}  combined {alpha_cut.sum()}')
print(f'  Beta:   opt {beta_cut_orig.sum()}  proj {beta_cut_proj.sum()}  combined {beta_cut.sum()}')

# ---- Diagnostic plot ----
fig2, axes2 = plt.subplots(2, 2, figsize=(10, 9))
fig2.suptitle(
    'Emcee calibration v2 + projected optical-only GRBs  (relative err cut: err/value > 1)\n'
    'Blue = overlap (calibration); green = projected kept; red = projected removed',
    fontsize=10)

proj_params_plot = {
    'logFa': ('log10Faopt', 'log10FaErr_opt', 'log10Fa',  'log10FaErr', fa_cut),
    'logTa': ('log10Taopt', 'logTaErr_opt',    'log10Ta',  'log10TaErr', ta_cut),
    'Alpha': ('Alpha_opt',  'AlphaErr_opt',    'Alpha',    'AlphaErr',   alpha_cut),
    'Beta':  ('Beta_opt',   'betaErr_opt',     'Beta',     'BetaErr',    beta_cut),
}

opt_only_grbs = out.index[opt_only_mask]
opt_proj_raw  = opt[opt['GRB'].isin(opt_only_grbs)].set_index('GRB')

for ax2, (param, (xc_o, xe_o, val_col, err_col, param_mask)) in zip(axes2.flat, proj_params_plot.items()):
    xc, xe, yc, ye = params_map[param]
    sub_cal = merged[[xc, xe, yc, ye]].apply(pd.to_numeric, errors='coerce').dropna()
    sub_cal = apply_calib_cuts(sub_cal, xc, xe, yc, ye, param)

    ax2.errorbar(sub_cal[xc], sub_cal[yc],
                 xerr=sub_cal[xe], yerr=sub_cal[ye],
                 fmt='o', ms=3, color='steelblue', alpha=0.5,
                 elinewidth=0.6, capsize=1.5, label=f'Overlap (N={len(sub_cal)})', zorder=2)

    r = mcmc_results[param]
    m_p, b_p = r['m'], r['b']
    x_all = pd.to_numeric(opt_proj_raw[xc_o], errors='coerce').dropna()
    x_lo  = x_all.min() if len(x_all) > 0 else sub_cal[xc].min()
    x_hi  = x_all.max() if len(x_all) > 0 else sub_cal[xc].max()
    x_line = np.linspace(x_lo, x_hi, 200)
    ax2.plot(x_line, m_p * x_line + b_p, 'k-', lw=1.5, zorder=3)
    ax2.plot(x_line, x_line, 'k--', lw=0.8, alpha=0.4)

    x_proj  = pd.to_numeric(opt_proj_raw[xc_o], errors='coerce')
    y_proj  = out.loc[out.index.isin(opt_proj_raw.index), val_col]
    ye_proj = out.loc[out.index.isin(opt_proj_raw.index), err_col]
    removed = param_mask.reindex(opt_proj_raw.index, fill_value=False)

    common = x_proj.index.intersection(y_proj.index)
    n_kept = n_rem = 0
    for grb in common:
        is_removed = bool(removed.get(grb, False))
        color = '#d62728' if is_removed else '#2ca02c'
        ax2.errorbar(float(x_proj[grb]), float(y_proj[grb]),
                     yerr=float(ye_proj[grb]),
                     fmt='^', ms=4, color=color, alpha=0.7,
                     elinewidth=0.7, capsize=2, zorder=1)
        if is_removed:
            n_rem += 1
        else:
            n_kept += 1

    from matplotlib.lines import Line2D
    handles, _ = ax2.get_legend_handles_labels()
    handles += [
        Line2D([0],[0], marker='^', color='w', markerfacecolor='#2ca02c', ms=6,
               label=f'Kept ({n_kept})'),
        Line2D([0],[0], marker='^', color='w', markerfacecolor='#d62728', ms=6,
               label=f'Removed ({n_rem})'),
    ]
    ax2.legend(handles=handles, fontsize=7, loc='upper left')

    xlab, ylab = axis_labels[param]
    ax2.set_xlabel(xlab, fontsize=9)
    ax2.set_ylabel(ylab, fontsize=9)
    ax2.set_title(param, fontsize=10)
    ax2.tick_params(labelsize=8)

plt.tight_layout()
plot_path = os.path.join(PLOT_DIR, 'emcee_projected_errcut_relative.png')
plt.savefig(plot_path, dpi=150)
plt.close()
print(f'Plot saved to {plot_path}')

# Save unfiltered base dataset
print(f'\nFinal dataset shape: {out.shape}')
print(f'Non-null counts:\n{out.notna().sum()}')
out.to_csv(OUT_FILE, index=True)
print(f'\nSaved to {OUT_FILE}')

# T90 completeness check
t90_missing_base = out['log10T90'].isna().sum()
t90_missing_cut  = out_cut['log10T90'].isna().sum()
print(f'\n=== T90 completeness check ===')
print(f'  Unfiltered ({out.shape[0]} GRBs):  log10T90 missing = {t90_missing_base}')
print(f'  Error-cut  ({out_cut.shape[0]} GRBs):  log10T90 missing = {t90_missing_cut}')
if t90_missing_base > 0 or t90_missing_cut > 0:
    print('  WARNING: missing T90 values detected!')
    print(out_cut[out_cut['log10T90'].isna()][['log10T90']].to_string())
else:
    print('  OK: all GRBs have log10T90')
