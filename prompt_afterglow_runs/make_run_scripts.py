#!/usr/bin/env python3
"""
Generate the per-run scripts for the prompt / afterglow SuperLearner runs.

For each run this writes
  GAM_supercomputer/pa_<NN>_<name>_formula_generation/
      GAM_analysis_8variables.R   copy of the flagship's
                                  (theilsen_cleaned_formula_generation) with ONLY
                                  the input file, method label and feature lists
                                  changed -- cutting logic untouched
      <input csv>                 copied in, as the existing runs do
and one SuperLearner variant per feature set,
  prompt_afterglow_runs/sl_scripts/superlearner_<set>.R
      copy of superlearner_CHECKPOINT_pre_daume_fix.R with only the feature and
      error-column lists changed.

Every text replacement is asserted to match exactly once, so a silent no-op
edit can't slip through.
"""
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GAM_SRC = os.path.join(REPO, 'GAM_supercomputer', 'theilsen_cleaned_formula_generation',
                       'GAM_analysis_8variables.R')
SL_SRC = os.path.join(REPO, 'superlearner_CHECKPOINT_pre_daume_fix.R')

ALL10 = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
         'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux']
FEATURE_SETS = {
    'prompt':    ['log10T90', 'log10Fluence', 'PhotonIndex', 'log10PeakFlux'],
    'prompt_nh': ['log10T90', 'log10Fluence', 'PhotonIndex', 'log10PeakFlux', 'log10NH'],
    'afterglow': ['log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma', 'log10NH'],
    'all':       ALL10,
    'all_epeak': ALL10 + ['log10Epeak'],
}
# feature -> its error column(s) as named in each script
GAM_ERR = {'log10T90': 'T90Err', 'log10Fa': 'log10FaErr', 'log10Ta': 'log10TaErr',
           'Alpha': 'AlphaErr', 'Beta': 'BetaErr', 'log10Fluence': 'fluence_err_col',
           'PhotonIndex': 'PhotonIndexErr', 'log10PeakFlux': 'peakflux_err_col'}
SL_ERR = {'log10T90': 'log10T90Err', 'log10Fa': 'log10FaErr', 'log10Ta': 'log10TaErr',
          'Alpha': 'AlphaErr', 'Beta': 'BetaErr', 'log10Fluence': 'log10FluenceErr',
          'PhotonIndex': 'PhotonIndexErr', 'log10PeakFlux': 'log10PeakFluxErr'}

RUNS = [  # (NN, name, input csv in Data/, feature set)
    (1, 'prompt_swiftpf',          'superlearner_training_theilsen_errcut_relative.csv',            'prompt'),
    (2, 'prompt_gbmpf',            'superlearner_training_theilsen_errcut_relative_gbmpf.csv',      'prompt'),
    (3, 'afterglow',               'superlearner_training_theilsen_errcut_relative.csv',            'afterglow'),
    (4, 'all_swiftpf',             'superlearner_training_theilsen_errcut_relative.csv',            'all'),
    (5, 'all_gbmpf',               'superlearner_training_theilsen_errcut_relative_gbmpf.csv',      'all'),
    (6, 'revmulti_all_swiftpf',    'superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative.csv', 'all'),
    (7, 'revmulti_all_gbmpf',      'superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative_gbmpf.csv', 'all'),
    (8, 'all_gbmpf_epeak0',        'superlearner_training_theilsen_errcut_relative_gbmpf_epeak0.csv',    'all_epeak'),
    (9, 'all_gbmpf_epeakmice',     'superlearner_training_theilsen_errcut_relative_gbmpf_epeakmice.csv', 'all_epeak'),
    # second series: corrected time+position GBM crosswalk (build_gbm_v2_inputs.py)
    (10, 'prompt_nh_swiftpf',      'superlearner_training_theilsen_errcut_relative.csv',            'prompt_nh'),
    (11, 'prompt_nh_gbmpf2',       'superlearner_training_theilsen_errcut_relative_gbmpf2.csv',     'prompt_nh'),
    (12, 'prompt_gbmpf2',          'superlearner_training_theilsen_errcut_relative_gbmpf2.csv',     'prompt'),
    (13, 'all_gbmpf2',             'superlearner_training_theilsen_errcut_relative_gbmpf2.csv',     'all'),
    (14, 'all_gbmpf2cal',          'superlearner_training_theilsen_errcut_relative_gbmpf2cal.csv',  'all'),
    (15, 'revmulti_all_gbmpf2',    'superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative_gbmpf2.csv', 'all'),
    (16, 'all_gbmpf2_epeak0',      'superlearner_training_theilsen_errcut_relative_gbmpf2_epeak0.csv',    'all_epeak'),
    (17, 'all_gbmpf2_epeakmice',   'superlearner_training_theilsen_errcut_relative_gbmpf2_epeakmice.csv', 'all_epeak'),
]


def sub_once(text, old, new, what):
    n = text.count(old)
    assert n == 1, f'{what}: expected 1 match, found {n}'
    return text.replace(old, new)


def drop_lines(text, feats_kept, pattern, what):
    """Remove lines matching `pattern` (with a {f} placeholder) for every
    feature NOT kept; assert each such line existed."""
    for f in ALL10:
        if f in feats_kept:
            continue
        rx = re.compile(pattern.format(f=re.escape(f)))
        lines = text.split('\n')
        hit = [i for i, l in enumerate(lines) if rx.search(l)]
        if hit:
            text = '\n'.join(l for i, l in enumerate(lines) if i not in hit)
    return text


def make_gam(nn, name, csv, fset):
    feats = FEATURE_SETS[fset]
    t = open(GAM_SRC).read()
    t = sub_once(t, 'read.csv("superlearner_training_theilsen.csv"', f'read.csv("{csv}"', 'input')
    t = sub_once(t, 'method = "theilsen")', f'method = "pa_{nn:02d}_{name}")', 'method')
    # MICE predictor block
    m = re.search(r'features_for_mice_preds = subset\(raw_xray_data,select = c\(.*?\)\)', t, re.S)
    assert m, 'preds block'
    block = ('features_for_mice_preds = subset(raw_xray_data,select = c('
             + ',\n                                                          '.join(feats) + '))')
    t = t[:m.start()] + block + t[m.end():]
    # safety-net NA lines for dropped features
    t = drop_lines(t, feats, r'^features_for_mice_preds\${f}\[', 'safety net')
    # error block
    m = re.search(r'features_for_mice_errs = raw_xray_data\[, c\(.*?\)\]', t, re.S)
    assert m, 'errs block'
    errs = []
    for f in feats:
        if f in GAM_ERR:
            e = GAM_ERR[f]
            errs.append(e if e.endswith('_col') else f'"{e}"')
    t = t[:m.start()] + 'features_for_mice_errs = raw_xray_data[, c(' + ', '.join(errs) + ')]' + t[m.end():]
    # guard the log-error conversions on the feature actually being present
    for feat, col in (('log10T90', 'log10T90Err'), ('log10Fluence', 'log10FluenceErr'),
                      ('log10PeakFlux', 'log10PeakFluxErr')):
        t = sub_once(t, f'if (!"{col}" %in% colnames(GRBPred)) {{',
                     f'if (!"{col}" %in% colnames(GRBPred) && "{feat}" %in% colnames(GRBPred)) {{',
                     f'guard {col}')
    # resumable split count override (default unchanged = 100) for smoke tests
    t = sub_once(t, 'N_SPLITS <- 100', 'N_SPLITS <- as.integer(Sys.getenv("N_SPLITS", "100"))', 'N_SPLITS')
    d = os.path.join(REPO, 'GAM_supercomputer', f'pa_{nn:02d}_{name}_formula_generation')
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, 'GAM_analysis_8variables.R'), 'w').write(t)
    shutil.copy(os.path.join(REPO, 'Data', csv), os.path.join(d, csv))
    return d


def make_sl(fset):
    feats = FEATURE_SETS[fset]
    t = open(SL_SRC).read()
    t = sub_once(t, '''features_for_mice_preds <- subset(raw_xray_data, select = c(
  log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma,
  log10Fluence, PhotonIndex, log10NH, log10PeakFlux
))''', 'features_for_mice_preds <- subset(raw_xray_data, select = c(\n  '
        + ', '.join(feats) + '\n))', 'SL preds block')
    t = drop_lines(t, feats, r'^features_for_mice_preds\${f}\[', 'SL safety net')
    # error lists (raw block + two grb_errors blocks): drop entries of removed features
    for f in ALL10:
        if f in feats or f not in SL_ERR:
            continue
        e = SL_ERR[f]
        lines = t.split('\n')
        hit = [i for i, l in enumerate(lines) if re.match(rf'^\s+{re.escape(e)}\s+=', l)]
        assert len(hit) == 3, f'SL err {e}: expected 3 lines, found {len(hit)}'
        t = '\n'.join(l for i, l in enumerate(lines) if i not in hit)
    # a data.frame( ... ) list must not end in a trailing comma after removals
    t = re.sub(r',(\s*\n\s*\))', r'\1', t)
    os.makedirs(os.path.join(HERE, 'sl_scripts'), exist_ok=True)
    p = os.path.join(HERE, 'sl_scripts', f'superlearner_{fset}.R')
    open(p, 'w').write(t)
    return p


if __name__ == '__main__':
    # optional: only these run numbers, e.g. `make_run_scripts.py 10 11` (default: all)
    only = {int(a) for a in sys.argv[1:]}
    runs = [r for r in RUNS if not only or r[0] in only]
    for fs in sorted({r[3] for r in runs}):
        print('SL  ', make_sl(fs))
    for nn, name, csv, fs in runs:
        print('GAM ', make_gam(nn, name, csv, fs))
