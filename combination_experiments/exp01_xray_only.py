#!/usr/bin/env python3
"""exp01: control -- no fusion at all. Train only on the 223 X-ray GRBs,
discard the 83 optical-only GRBs entirely. Tests whether adding the
optical-only GRBs (via any imputation method) is even worth doing."""
import os
import numpy as np
from common import load_data, xray_linear_errors, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
fluence_err_lin, peakflux_err_lin = xray_linear_errors(xray)

import pandas as pd
out = pd.DataFrame({
    'Redshift_crosscheck': xray['z'],
    'log10T90':            xray['log10T90'],
    'log10Fa':             xray['logFa_x'],
    'log10FaErr':          xray['logFaErr_x'],
    'log10Ta':             xray['logTa_x'],
    'log10TaErr':          xray['logTaErr_x'],
    'Alpha':               xray['Alpha_x'],
    'AlphaErr':            xray['AlphaErr_x'],
    'Beta':                xray['Beta_x'],
    'BetaErr':             xray['BetaErr_x'],
    'Gamma':               xray['Gamma'],
    'log10Fluence':        xray['log10Fluence'],
    'PhotonIndex':         xray['PhotonIndex'],
    'log10NH':             xray['log10NH'],
    'log10PeakFlux':       xray['log10PeakFlux'],
    'T90Err':              xray.get('T90Err', np.nan),
    'FluenceErr':          fluence_err_lin,
    'PhotonIndexErr':      xray.get('PhotonIndexErr', np.nan),
    'PeakFluxErr':         peakflux_err_lin,
}, index=xray.index)
out = out[pd.to_numeric(out['Redshift_crosscheck'], errors='coerce').notna()]
out.index.name = 'GRB'
out_path = os.path.join(SCRIPT_DIR, 'results', 'exp01_xray_only.csv')
out.to_csv(out_path)
print(f'exp01_xray_only: {len(out)} GRBs -> {out_path}')
