#' =============================================================================
#' superlearner.R -- GRB redshift prediction via a SuperLearner ensemble
#' =============================================================================
#'
#' Estimates the redshift (z) of gamma-ray bursts (GRBs) from their X-ray
#' afterglow and prompt-emission features, using a stacked SuperLearner ensemble
#' with Monte-Carlo (MC) propagation of the input measurement errors.
#'
#' Pipeline stages (in the order they appear below):
#'   1. Ingest & normalize   -- read the input CSV; map alternate column
#'                              conventions (X -> GRB, linear T90 -> log10T90);
#'                              key every row by its GRB id.
#'   2. Recover & clean      -- restore real prompt features for optical-projected
#'                              GRBs from Data/optical_data.csv; for pre-fused
#'                              frames (superlearner_training_emcee_v2*.csv, which
#'                              arrive complete) flag optical rows by X-ray-catalog
#'                              membership instead; drop short GRBs (log10T90 cut);
#'                              null physically-impossible values.
#'   3. Impute               -- predictors: missForest co-imputed with the unlabeled
#'                              generalization pool, repeated over imputation seeds x
#'                              emcee projection draws (predictions are averaged over
#'                              all frames; pre-fused frames use the 8 v2-chain draws
#'                              in Data/proj_draws_v2.csv); errors: MICE (midastouch,
#'                              m = 20).
#'   4. Feature engineering  -- order predictors by LASSO importance; squares of the
#'                              top-7; is_optical indicator + scaled per-domain
#'                              feature copies (Daume-style, c = 0.25).
#'   5. Outlier handling     -- optional M-estimator catastrophic-outlier cut
#'                              (m_estimator.R) and/or upsampling (upsampling.R).
#'   6. Split                -- seeded random 20% holdout / 80% train.
#'   7. Cross-validated fit  -- `loop` repetitions of k-fold CV SuperLearner,
#'                              parallelized across repetitions with forked
#'                              mclapply. Each fold also MC-perturbs its test
#'                              points through the fitted ensemble.
#'   8. Aggregate & report   -- pool MC draws, compute coverage / correlation /
#'                              RMSE, write OutputFiles/*.csv + Plot_Output/*.png,
#'                              and persist the final full-data model.
#'
#' Invocation:
#'   Rscript superlearner.R <input.csv> <do_mice> <upsampling> \
#'           <do_m_estimator> <custom_models> <weight_threshold> <loop>
#'   Booleans are the literal strings "true"/"false" (see the arg block below).
#'   Env vars SMOKE_TEST / SMOKE_N / ONE_FOLD enable fast structural runs;
#'   REMOVE_CAT_OUTLIERS=false disables the 2-sigma catastrophic-outlier re-run.
#'
#' Input contract: the CSV must provide a GRB id column, Redshift_crosscheck, the
#'   10 predictors (log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma, log10Fluence,
#'   PhotonIndex, log10NH, log10PeakFlux) and their 8 error columns (T90Err,
#'   log10FaErr, log10TaErr, AlphaErr, BetaErr, FluenceErr, PhotonIndexErr,
#'   PeakFluxErr).
#'
#' Outputs: OutputFiles/{grb_xray_imputed.csv, mc_predictions.csv, mc_raw.rds},
#'   Plot_Output/*.png, and the serialized model file `superlearner_model`.
#' =============================================================================

source("Load_Imports.R")          # package loads + shared helpers
source("Result_plot_maker.R")     # result_plotter()
source("mc_error_propagation.R")  # mc_predict(), summarize_mc(), inside_mc_cone()
source("mc_plots.R")              # make_mc_plots()


# ---- Run configuration ------------------------------------------------------
#' Configuration is sourced either from hard-coded debug defaults (`run_locally`)
#' or from the 7 positional CLI args passed by app.py. After this block the
#' following are defined:
#'   raw_xray_data    : data.frame  -- the input table, rows keyed by GRB id
#'   do_mice          : logical     -- impute missing values with MICE?
#'   upsampling       : logical     -- run upsampling.R before fitting?
#'   do_m_estimator   : logical     -- run the M-estimator outlier cut?
#'   custom_models    : logical     -- read the learner library from selected_models.txt?
#'   weight_threshold : numeric     -- M-estimator weight below which a GRB is dropped
#'   loop             : integer     -- number of CV repetitions
#'   m_est_pct        : numeric     -- percentile (0–1) of weights to drop; overrides
#'                                     the hardcoded 5% default in m_estimator.R
run_locally <- FALSE  # TRUE only for interactive debugging
if (run_locally) {
  raw_xray_data    <- read.csv("combined_data_with_redshift_V8.csv", header = TRUE, row.names = 1)
  do_mice          <- TRUE
  upsampling       <- FALSE
  do_m_estimator   <- TRUE
  custom_models    <- FALSE
  weight_threshold   <- 0.65
  m_est_pct          <- 0.05
  use_formula_learners <- FALSE
  loop               <- 5
} else {
  print("not running locally")
  args             <- commandArgs(trailingOnly = TRUE)
  input_file       <- args[1]                              # character: path to input CSV
  do_mice          <- as.logical(tolower(args[2]) == "true")
  upsampling       <- as.logical(tolower(args[3]) == "true")
  do_m_estimator   <- as.logical(tolower(args[4]) == "true")
  custom_models    <- as.logical(tolower(args[5]) == "true")
  weight_threshold <- as.numeric(args[6])
  loop             <- as.numeric(args[7])
  out_dir          <- if (length(args) >= 8) args[8] else "."
  m_est_pct        <- if (length(args) >= 9) as.numeric(args[9]) else 0.05
  use_formula_learners <- if (length(args) >= 10) as.logical(tolower(args[10]) == "true") else FALSE
  # outlier_method: "hard" (default, percentile cut), "soft" (obsWeights, no cut),
  #                 "infold" (in-fold rlm wrapper, no pre-cut), "both" (hard cut + infold wrapper)
  outlier_method   <- if (length(args) >= 11) tolower(args[11]) else "hard"
  raw_xray_data    <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

  #' Input normalization. Some files (e.g. x-ray_data.csv) use alternate column
  #' conventions: the GRB identifier sits in an unnamed first column (read in as
  #' "X") rather than a named "GRB" column, and T90 is on a linear scale rather
  #' than log10T90. Map them here so the rest of the pipeline -- which expects
  #' `GRB` + `log10T90` -- can ingest them. Both branches are no-ops for files
  #' that already conform (e.g. training_data.csv), so behavior is unchanged.
  if (!"GRB" %in% names(raw_xray_data) && "X" %in% names(raw_xray_data)) {
    raw_xray_data$GRB <- raw_xray_data$X
  }
  if (!"log10T90" %in% names(raw_xray_data) && "T90" %in% names(raw_xray_data)) {
    raw_xray_data$log10T90 <- log10(raw_xray_data$T90)
  }
  #' Some newer catalogs ship dex-scale errors directly (log10FluenceErr etc.)
  #' instead of the linear errors (FluenceErr etc.) that to_dex_err() below
  #' expects. Reconstruct the linear error so to_dex_err() round-trips back to
  #' the same dex value; a no-op when the linear column is already present.
  if (!"T90Err" %in% names(raw_xray_data) && "log10T90Err" %in% names(raw_xray_data)) {
    raw_xray_data$T90Err <- raw_xray_data$log10T90Err * 10^raw_xray_data$log10T90 * log(10)
  }
  if (!"FluenceErr" %in% names(raw_xray_data) && "log10FluenceErr" %in% names(raw_xray_data)) {
    raw_xray_data$FluenceErr <- raw_xray_data$log10FluenceErr * 10^raw_xray_data$log10Fluence * log(10)
  }
  if (!"PeakFluxErr" %in% names(raw_xray_data) && "log10PeakFluxErr" %in% names(raw_xray_data)) {
    raw_xray_data$PeakFluxErr <- raw_xray_data$log10PeakFluxErr * 10^raw_xray_data$log10PeakFlux * log(10)
  }

  stopifnot(!anyDuplicated(raw_xray_data$GRB))
  rownames(raw_xray_data) <- raw_xray_data$GRB
}
#' Do NOT re-read input_file after this point: the read above already loaded
#' raw_xray_data and set rownames <- GRB. A second read would crash in
#' run_locally mode (input_file undefined) and would reset rownames to 1..N,
#' discarding the GRB-keyed identity that grb_errors, the MC lookups and the GRB
#' name lookup all depend on.


# ---- Recover real prompt features for optical-projected GRBs ----------------
# Optical-only GRBs enter with 6/10 predictors missing; T90, Fluence, PhotonIndex,
# NH and PeakFlux exist in the optical catalog (unit conversions validated on the
# 88 overlap GRBs: Fluence in 1e-7 erg/cm2, NH in 1e21 cm-2). Gamma has no optical
# counterpart and stays NA for imputation. Native plateau values are kept for the
# emcee projection draws applied at the imputation stage.
pred_vars <- c("log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
               "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux")
raw_xray_data$is_optical <- 0
optical_native <- NULL
if (file.exists("Data/optical_data.csv")) {
  n_miss   <- rowSums(is.na(raw_xray_data[, pred_vars]))
  opt_rows <- which(n_miss >= 6)
  if (length(opt_rows) > 0) {
    opt_cat <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE)
    names(opt_cat)[1] <- "GRB"
    m  <- match(raw_xray_data$GRB[opt_rows], opt_cat$GRB)
    ok <- !is.na(m); opt_rows <- opt_rows[ok]; m <- m[ok]
    raw_xray_data$log10T90[opt_rows]      <- log10(opt_cat$T90[m])
    raw_xray_data$log10Fluence[opt_rows]  <- log10(opt_cat$Fluence[m]) - 7
    raw_xray_data$PhotonIndex[opt_rows]   <- opt_cat$PhotonIndex[m]
    raw_xray_data$log10NH[opt_rows]       <- log10(opt_cat$NH[m]) + 21
    raw_xray_data$log10PeakFlux[opt_rows] <- log10(ifelse(opt_cat$PeakFlux[m] > 0, opt_cat$PeakFlux[m], NA))
    raw_xray_data$Gamma[opt_rows]         <- NA
    raw_xray_data$is_optical[opt_rows]    <- 1
    optical_native <- data.frame(logFa = opt_cat$logFa[m], logTa = opt_cat$logT_a[m],
                                 Alpha = opt_cat$Alpha[m], Beta = opt_cat$Beta[m],
                                 row.names = raw_xray_data$GRB[opt_rows])
    cat("Recovered optical-catalog prompt features for", length(opt_rows), "GRBs\n")
  }
}

#' Pre-fused frames (e.g. Data/superlearner_training_emcee_v2*.csv) arrive with
#' the optical rows already projected to X-ray units and complete, so the
#' missingness test above finds nothing. Identify optical rows by X-ray-catalog
#' membership instead (the convention validated in
#' outlier_experiments/cv_experiments/campaign2/cv_maxstack_v2.R) and keep their
#' native optical plateau features so the emcee projection draws can re-project
#' them per posterior draw at the imputation stage. `prefused_frame` also
#' switches the draw source to the v2 emcee chains below.
prefused_frame <- FALSE
if (is.null(optical_native) &&
    all(pred_vars %in% names(raw_xray_data)) &&
    !any(rowSums(is.na(raw_xray_data[, pred_vars])) >= 6) &&
    file.exists("Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv") &&
    file.exists("Data/OnlyLGRBs_data_171_optical_corrected.csv")) {
  norm_grb <- function(s) {
    s <- trimws(gsub("GRB", "", as.character(s)))
    ifelse(grepl("[A-Za-z]$", s), s, paste0(s, "A"))
  }
  xr_ids <- norm_grb(read.csv("Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv",
                              stringsAsFactors = FALSE)[[1]])
  ids      <- norm_grb(raw_xray_data$GRB)
  opt_rows <- which(!(ids %in% xr_ids))
  if (length(opt_rows) > 0) {
    prefused_frame <- TRUE
    raw_xray_data$is_optical[opt_rows] <- 1
    opt_cat <- read.csv("Data/OnlyLGRBs_data_171_optical_corrected.csv", stringsAsFactors = FALSE)
    names(opt_cat)[1] <- "GRB"; opt_cat$GRB <- norm_grb(opt_cat$GRB)
    m  <- match(ids[opt_rows], opt_cat$GRB)
    ok <- !is.na(m)
    optical_native <- data.frame(logFa = opt_cat$log10Faopt[m[ok]], logTa = opt_cat$log10Taopt[m[ok]],
                                 Alpha = opt_cat$Alpha_opt[m[ok]], Beta = opt_cat$Beta_opt[m[ok]],
                                 row.names = rownames(raw_xray_data)[opt_rows[ok]])
    cat("Pre-fused frame:", length(opt_rows), "optical GRBs flagged by catalog membership;",
        sum(ok), "matched to the optical catalog for projection draws\n")
  }
}


# ---- Fast-run switches (env-controlled) -------------------------------------
#' smoke_test : logical -- run the whole pipeline on a small random subsample and
#'   a 2-learner library so an end-to-end check finishes in seconds.
#'     SMOKE_TEST=true SMOKE_N=120 Rscript superlearner.R <args...>
#' one_fold   : logical -- keep the full k-fold split (each model still trains on
#'   ~90% of GRBs) but run only fold 1, ~k x faster. Only that held-out fold gets
#'   CV + MC predictions; the rest stay as placeholder zeros, so the legacy CV
#'   correlation plot is not meaningful under ONE_FOLD (the MC outputs are).
#'     ONE_FOLD=true Rscript superlearner.R <args...>
smoke_test <- tolower(Sys.getenv("SMOKE_TEST", "false")) == "true"
smoke_n    <- as.integer(Sys.getenv("SMOKE_N", "120"))   # integer: subsample size
one_fold   <- tolower(Sys.getenv("ONE_FOLD", "false")) == "true"
if (smoke_test) {
  set.seed(42)
  keep <- sample(nrow(raw_xray_data), min(smoke_n, nrow(raw_xray_data)))
  raw_xray_data <- raw_xray_data[keep, ]
  cat(sprintf("SMOKE TEST: subsampled to %d GRBs\n", nrow(raw_xray_data)))
}

cat("Column names found:\n")
print(colnames(raw_xray_data))


#' Append a squared term for every column of a numeric data.frame.
#'
#' @param inputData data.frame of numeric predictors.
#' @return data.frame: `inputData` unchanged, plus one appended `<name>Sqr`
#'   column per input column (e.g. log10Fa -> log10Fa, log10FaSqr). Originals are
#'   never dropped, so no linear information is lost and there is no double
#'   squaring when called on an already-linear predictor set.
SqrTermGen <- function(inputData) {
  indVar <- colnames(inputData)
  for (i in seq_along(indVar)) {
    for (j in i:length(indVar)) {
      if (indVar[i] == indVar[j]) {  # square: the i == j case is the only term kept
        inputData[[paste0(indVar[i], "Sqr")]] <- inputData[, indVar[i]] * inputData[, indVar[j]]
      }
    }
  }
  inputData
}


# ---- Output locations & plot styling ----------------------------------------
sz   <- 0.8    # numeric: plot size multiplier (used by Result_plot_maker.R)
rez  <- 120    # integer: plot resolution in dpi (used by Result_plot_maker.R)
if (!exists("out_dir")) out_dir <- "."
addr     <- file.path(out_dir, "Results/")
PLOTaddr <- file.path(out_dir, "Plot_Output/")
if (!dir.exists(PLOTaddr)) dir.create(PLOTaddr, recursive = TRUE)
if (!dir.exists(addr))     dir.create(addr,     recursive = TRUE)


# ---- Clean & filter ---------------------------------------------------------
#' Keep long GRBs only: drop bursts with log10T90 <= 0.301 (T90 <= 2 s, the
#' short/long divide). Rows with a missing log10T90 are retained so MICE can
#' impute them rather than being silently discarded here.
raw_xray_data <- raw_xray_data[is.na(raw_xray_data$log10T90) | raw_xray_data$log10T90 > 0.301, ]

#' Predictor matrix fed to MICE: the 10 numeric features, independent of the
#' response (Redshift_crosscheck).
features_for_mice_preds <- subset(raw_xray_data, select = c(
  log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma,
  log10Fluence, PhotonIndex, log10NH, log10PeakFlux
))

#' Convert a linear-scale 1-sigma error to a dex (log10) error.
#'
#' For a quantity stored as log10(x), the error propagates as
#'   d(log10 x) = dx / (x * ln 10).
#' T90/Fluence/PeakFlux are stored as log10 features but ship with linear errors;
#' converting BEFORE imputation (a) puts every error column in the same
#' O(0.01-1) magnitude range, so mice's "constant"/near-zero-variance screen no
#' longer drops FluenceErr, and (b) yields exactly the dex quantities grb_errors
#' and the MC step consume, so no post-imputation conversion is needed.
#'
#' @param lin_err   numeric vector: linear-scale 1-sigma error (dx).
#' @param log10_val numeric vector: log10 of the measured value (log10 x).
#' @return numeric vector of dex-scale errors; non-finite results (missing or
#'   zero value/error) are set to NA so mice imputes them.
to_dex_err <- function(lin_err, log10_val) {
  out <- lin_err / (10^log10_val * log(10))
  out[!is.finite(out)] <- NA
  out
}

#' Error matrix fed to MICE. log10FaErr/log10TaErr are already dex;
#' AlphaErr/BetaErr/PhotonIndexErr are errors on linear spectral indices (no log)
#' and pass through unchanged; the three linear errors are converted via to_dex_err.
features_for_mice_errs <- data.frame(
  log10T90Err      = to_dex_err(raw_xray_data$T90Err,     raw_xray_data$log10T90),
  log10FaErr       = raw_xray_data$log10FaErr,
  log10TaErr       = raw_xray_data$log10TaErr,
  AlphaErr         = raw_xray_data$AlphaErr,
  BetaErr          = raw_xray_data$BetaErr,
  log10FluenceErr  = to_dex_err(raw_xray_data$FluenceErr, raw_xray_data$log10Fluence),
  PhotonIndexErr   = raw_xray_data$PhotonIndexErr,
  log10PeakFluxErr = to_dex_err(raw_xray_data$PeakFluxErr, raw_xray_data$log10PeakFlux),
  row.names        = rownames(raw_xray_data)
)

#' Null out non-finite / physically-impossible feature values so MICE re-imputes
#' them rather than fitting on bad data.
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA
features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20]       <- NA  # column density floor
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 3]              <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3]            <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3]            <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA


# ---- Impute ------------------------------------------------------------------
# Predictors: missForest, co-imputed with the unlabeled generalization pool so the
# imputation model learns the joint feature distribution from ~2x the rows (never
# sees redshift). Repeated over IMP_SEEDS; for optical rows the projected plateau
# features additionally vary over emcee posterior calibration draws. Downstream,
# predictions are averaged over all seed x draw frames. Errors: MICE as before.
if (do_mice) {
  require(missForest)
  set.seed(1)
  png(filename = file.path(out_dir, "MICE_missing_features.png"), width = 1000, height = 1000, res = 200)
  md.pattern(features_for_mice_preds, rotate.names = TRUE)
  title(main = "Missing-data pattern before imputation",
        sub  = "Rows = observed patterns (left count = # GRBs); blue = observed, red = missing; right/bottom = # missing")
  dev.off()

  gen_pool <- NULL
  if (file.exists("Data/TOTAL_GENERALIZATION_DATA_v4.csv")) {
    gen <- read.csv("Data/TOTAL_GENERALIZATION_DATA_v4.csv", row.names = 1, stringsAsFactors = FALSE)
    gen_pool <- data.frame(
      log10T90 = gen$logT90, log10Fa = gen$Fbest, log10Ta = gen$T_abest,
      Alpha = gen$Alpha, Beta = NA_real_, Gamma = gen$Gamma,
      log10Fluence = NA_real_,
      PhotonIndex = suppressWarnings(as.numeric(gen$photon_index)),
      log10NH = gen$logNH,
      log10PeakFlux = suppressWarnings(as.numeric(gen$logPeakFlux)))
    gen_pool <- gen_pool[, colnames(features_for_mice_preds)]
    gen_pool$log10NH[gen_pool$log10NH < 20]          <- NA
    gen_pool$PhotonIndex[gen_pool$PhotonIndex < 0]   <- NA
    gen_pool$Gamma[gen_pool$Gamma > 3]               <- NA
  }
  #' Projection-draw source. Pre-fused frames use the v2 emcee-chain draws
  #' (Data/proj_draws_v2.csv: 8 parameter-only draws from emcee_chains_v2.npz,
  #' the validated max-stack config — cv CV r 0.643 pre / 0.683 post outlier
  #' removal) with a single imputation seed, since those frames arrive complete
  #' and the seed sweep would only jitter a handful of leftover cells. Old-style
  #' frames keep the original 4-draw file and 3-seed sweep so their validated
  #' behavior is unchanged. Intrinsic-scatter draws are never used (adding
  #' scatter noise measurably degrades r).
  if (prefused_frame && file.exists("Data/proj_draws_v2.csv")) {
    proj_draws <- read.csv("Data/proj_draws_v2.csv", stringsAsFactors = FALSE)
    n_draws    <- 8L
    IMP_SEEDS  <- 12
    cat("Using v2 emcee projection draws (8 draws, 1 imputation seed)\n")
  } else {
    proj_draws <- if (file.exists("Data/emcee_projection_draws.csv"))
      read.csv("Data/emcee_projection_draws.csv", stringsAsFactors = FALSE) else NULL
    n_draws    <- 4L
    IMP_SEEDS  <- c(12, 77, 301)
  }
  draw_map   <- c(logFa = "log10Fa", logTa = "log10Ta", Alpha = "Alpha", Beta = "Beta")
  opt_here   <- intersect(rownames(features_for_mice_preds),
                          if (is.null(optical_native)) character(0) else rownames(optical_native))

  PROJ_DRAWS <- if (!is.null(proj_draws) && length(opt_here) > 0) seq_len(n_draws) else 0
  if (smoke_test) { IMP_SEEDS <- IMP_SEEDS[1]; PROJ_DRAWS <- PROJ_DRAWS[1] }

  feature_frames <- list()
  for (d in PROJ_DRAWS) {
    fp <- features_for_mice_preds
    if (d > 0) {
      for (p in names(draw_map)) {
        dr  <- proj_draws[proj_draws$param == p & proj_draws$draw == d - 1, ]
        val <- dr$m * optical_native[opt_here, p] + dr$b
        fin <- is.finite(val)
        fp[opt_here[fin], draw_map[[p]]] <- val[fin]
      }
    }
    for (s in IMP_SEEDS) {
      set.seed(s)
      imp <- missForest(rbind(fp, gen_pool))$ximp[seq_len(nrow(fp)), ]
      rownames(imp) <- rownames(fp)
      feature_frames[[length(feature_frames) + 1]] <- imp
    }
  }
  cat("Built", length(feature_frames), "imputation frames (seeds x projection draws)\n")
  features_for_mice_preds <- feature_frames[[1]]

  mice_model_errs        <- mice(data = features_for_mice_errs, m = 20, method = "midastouch", printFlag = FALSE)
  features_for_mice_errs <- complete(mice_model_errs, 20)

  GRBPred <- cbind(features_for_mice_preds, features_for_mice_errs)
} else {
  GRBPred <- na.omit(cbind(features_for_mice_preds, features_for_mice_errs))
  feature_frames <- list(GRBPred[, colnames(features_for_mice_preds)])
}

#' Assertion: GRBPred must be complete after imputation. anyNA() catches NaN too.
#' If this fires with do_mice = TRUE, MICE failed to fill some cell (e.g. an
#' all-NA column); with do_mice = FALSE, na.omit() did not drop every gap.
stopifnot("GRBPred still contains NA/NaN after imputation step" = !anyNA(GRBPred))
cat("Assertion OK: no NA in GRBPred after imputation (", nrow(GRBPred), "rows x", ncol(GRBPred), "cols )\n")

GRBPred$GRB <- raw_xray_data[rownames(GRBPred), "GRB"]

#' MC error frame, keyed by the same row identity GRBPred uses downstream so the
#' fold loop can look up errors directly with rownames(test_set). Error columns
#' are already dex (see features_for_mice_errs), so no conversion happens here.
grb_errors <- data.frame(
  log10FaErr       = GRBPred$log10FaErr,
  log10TaErr       = GRBPred$log10TaErr,
  PhotonIndexErr   = GRBPred$PhotonIndexErr,
  log10PeakFluxErr = GRBPred$log10PeakFluxErr,
  log10T90Err      = GRBPred$log10T90Err,
  AlphaErr         = GRBPred$AlphaErr,
  BetaErr          = GRBPred$BetaErr,
  log10FluenceErr  = GRBPred$log10FluenceErr,
  row.names        = rownames(GRBPred),
  stringsAsFactors = FALSE
)
cat("MC errors frame built for", nrow(grb_errors), "rows\n")


# ---- Feature engineering ----------------------------------------------------
source("lasso.R")  # defines `lassovar`: predictor names ordered by mean |LASSO coef|, strongest first
#' `lassovar` holds ALL predictor columns ordered by importance. To use the
#' leaner 7-feature model the Best_formula_*.txt formulas were tuned for, restore
#' the top-7 cut below; left disabled so every feature with an error column
#' (incl. Alpha/Beta/log10Fluence) stays in the model and its MC perturbation
#' actually propagates. Downstream code only references columns that exist, so
#' either choice is safe.
# lassovar <- head(lassovar, 7)

#' Build the squared-term design: square the (linear) LASSO-ordered predictors,
#' then re-attach the response/error columns SqrTermGen does not touch.
Variables          <- subset(GRBPred, select = lassovar)
Responses_and_Err  <- subset(GRBPred, select = !colnames(GRBPred) %in% colnames(Variables))
GRBPred            <- SqrTermGen(Variables)
GRBPred            <- cbind(GRBPred, Responses_and_Err)

#' Restore the response columns: Redshift_crosscheck (linear z) and its log10z =
#' log10(z + 1) transform, the SuperLearner target.
GRBPred$Redshift_crosscheck <- raw_xray_data$Redshift_crosscheck
GRBPred$log10z              <- log10(GRBPred$Redshift_crosscheck + 1)

out_files_dir <- file.path(out_dir, "OutputFiles")
if (!dir.exists(out_files_dir)) dir.create(out_files_dir, recursive = TRUE)
#' Snapshot the imputed/engineered frame; m_estimator.R reads it back.
if (do_mice) {
  write.csv(GRBPred, file.path(out_files_dir, "grb_xray_imputed.csv"))
} else {
  write.csv(GRBPred, file.path(out_files_dir, "grb_xray.csv"))
}


# ---- Outlier handling (optional) --------------------------------------------
if (upsampling) source("upsampling.R")

# "hard" and "both": run the standard M-estimator hard cut (drops bottom m_est_pct rows).
# "soft": compute rlm weights but keep all rows; weights fed to SuperLearner as obsWeights.
# "infold": skip pre-cut entirely; outlier removal happens inside each learner wrapper per fold.
sl_obs_weights <- NULL   # NULL => uniform weights in SuperLearner call

if (do_m_estimator && outlier_method %in% c("hard", "both")) {
  source("m_estimator.R")   # hard cut; overwrites GRBPred with GRBCut
} else if (do_m_estimator && outlier_method == "soft") {
  require(MASS)
  rlm_form_soft <- as.formula(paste("log10z ~", paste(lassovar, collapse = "+")))
  M_est_soft    <- MASS::rlm(rlm_form_soft,
                              data   = cbind(GRBPred, log10z = log10(GRBPred$Redshift_crosscheck + 1)),
                              method = "M", maxit = 50)
  sl_obs_weights <- setNames(M_est_soft$w, rownames(GRBPred))
  cat("Soft weights: min=", round(min(sl_obs_weights), 3),
      " median=", round(median(sl_obs_weights), 3),
      " (", sum(sl_obs_weights < quantile(sl_obs_weights, m_est_pct)), "GRBs below",
      m_est_pct * 100, "%-ile threshold)\n")
}
# "infold" and "both": define per-fold rlm detector + memoised cache + learner wrappers.
# Wrapped versions are registered in the global env so SuperLearner can resolve them.
if (outlier_method %in% c("infold", "both")) {
  require(MASS)
  require(digest)

  .detect_rlm <- function(Y, X, k = 3) {
    df  <- data.frame(Y = Y, X, check.names = FALSE)
    fit <- tryCatch(MASS::rlm(Y ~ ., data = df, method = "M", maxit = 50),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(TRUE, length(Y)))
    r <- residuals(fit)
    abs(r) <= k * mad(r, constant = 1.4826)
  }

  .store <- new.env(parent = emptyenv())
  .cached_detect <- function(Y, X) {
    key <- digest::digest(list(Y, X))
    hit <- .store[[key]]
    if (is.null(hit)) { hit <- .detect_rlm(Y, X); .store[[key]] <- hit }
    hit
  }

  .make_outlier_wrapper <- function(base_fn) {
    force(base_fn)
    function(Y, X, newX, family, obsWeights, id, ...) {
      keep <- .cached_detect(Y, as.data.frame(X))
      if (sum(keep) < max(20L, ncol(X) + 2L)) keep <- rep(TRUE, length(Y))
      cat(sprintf("  [infold] fold size %d → kept %d after rlm filter\n", length(Y), sum(keep)))
      base_fn(Y          = Y[keep],
              X          = X[keep, , drop = FALSE],
              newX       = newX,
              family     = family,
              obsWeights = obsWeights[keep],
              id         = id[keep],
              ...)
    }
  }
  cat("In-fold rlm wrapper defined.\n")
}

#' Assertion on the frame that actually feeds model fitting (post upsampling /
#' M-estimator). anyNA() covers NA and NaN; the is.finite check additionally
#' catches +/-Inf in the log-error columns, which would poison MC propagation.
.grb_num <- GRBPred[vapply(GRBPred, is.numeric, logical(1))]
stopifnot(
  "GRBPred contains NA/NaN after outlier removal"   = !anyNA(GRBPred),
  "GRBPred contains non-finite (Inf) numeric values" =
    all(vapply(.grb_num, function(col) all(is.finite(col)), logical(1)))
)
rm(.grb_num)
cat("Assertion OK: GRBPred has no NA/NaN/Inf (", nrow(GRBPred), "rows )\n")


# ---- Split ------------------------------------------------------------------
Responses <- subset(GRBPred, select = c("Redshift_crosscheck", "log10z"))

# Design matrix per imputation frame: all 10 linear predictors + squares of the
# top-7 by LASSO + is_optical + scaled per-domain copies (Daume augmentation,
# c = 0.25: penalized learners shrink domain-specific deviations harder than
# shared effects). One design per frame; frame 1 defines the modeling frame.
daume_c <- 0.25
build_design <- function(feats) {
  X <- feats[rownames(GRBPred), lassovar, drop = FALSE]
  for (v in head(lassovar, 7)) X[[paste0(v, "Sqr")]] <- X[[v]]^2
  X$is_optical <- raw_xray_data[rownames(X), "is_optical"]
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * (X$is_optical * daume_c)
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}
design_frames <- lapply(feature_frames, build_design)
GRBPred       <- cbind(design_frames[[1]], Responses)

#' Seeded random 20% holdout / 80% train. Selection is by row index then resolved
#' by rowname (safe because GRB rownames are unique). A random split avoids the
#' bias of a last-20% holdout if rows carry any ordering.
set.seed(42)
holdout_idx    <- sample(nrow(GRBPred), size = floor(0.20 * nrow(GRBPred)))
PredictionData <- GRBPred[holdout_idx, ]
TrainingData   <- GRBPred[!(rownames(GRBPred) %in% rownames(PredictionData)), ]
stopifnot(length(intersect(rownames(PredictionData), rownames(TrainingData))) == 0)

Response   <- TrainingData$log10z                                            # numeric vector: training target
Predictors <- subset(TrainingData, select = -c(log10z, Redshift_crosscheck)) # data.frame: training design matrix


# ---- Learner library --------------------------------------------------------
#' Custom SuperLearner wrappers. Each is sourced (base wrapper) then turned into
#' one learner per tuned formula via create.Learner. They must be defined in this
#' (global) environment so the forked CV workers inherit them and the final
#' full-data fit can resolve them -- a source() inside a fork does not propagate
#' back across the fork boundary. Each wrapper refits on every training fold (the
#' formula is fixed structure; coefficients are estimated per fold), so there is
#' no leakage.
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

formula_table_GAM <- read.table("Best_formula_GAM.txt")
bestGAM1          <- apply(as.matrix(formula_table_GAM[, 2]), 1, as.formula)  # list of formula
learner1          <- create.Learner("SL.mgcv_gam", tune = list(gam.model = c(bestGAM1)),
                                     detailed_names = FALSE, name_prefix = "gam")
libnames <- "_winner_"

# Winner library (paired-CV validated): first tuned GAM formula, band-interaction
# GLM with per-domain plateau slopes, all-feature GLM, glmnet, random forest,
# boosted GAM, mean. Stacked with non-negative ridge (regularized NNLS).
require(mboost)
band_form <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)),
                       detailed_names = FALSE, name_prefix = "glmB")
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df  <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = as.numeric(predict(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = newdata, type = "response"))
SL.gamboost <- function(Y, X, newX, family, obsWeights, ...) {
  df  <- data.frame(Y = Y, X, check.names = FALSE)
  num <- intersect(pred_vars, colnames(X))
  fm  <- as.formula(paste("Y ~", paste(c(sprintf("bbs(%s, df = 2)", num), "bols(is_optical)"), collapse = "+")))
  fit <- mboost::gamboost(fm, data = df, control = mboost::boost_control(mstop = 300, nu = 0.1))
  cvr <- try(mboost::cvrisk(fit, folds = mboost::cv(model.weights(fit), type = "kfold", B = 5),
                            papply = lapply), silent = TRUE)
  if (!inherits(cvr, "try-error")) fit <- fit[mboost::mstop(cvr)]
  out <- list(object = fit); class(out) <- "SL.gamboost"
  list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE))), fit = out)
}
predict.SL.gamboost <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE)))
method.NNRidge <- function() {
  out <- SuperLearner::method.NNLS()
  out$computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    fit <- glmnet::cv.glmnet(Z, Y, alpha = 0, lower.limits = 0, nfolds = 10)
    co  <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    co[co < 0] <- 0
    if (sum(co) > 0) co <- co / sum(co) else co <- rep(1 / ncol(Z), ncol(Z))
    list(cvRisk = apply(Z, 2, function(p) mean((p - Y)^2)), coef = co, optimizer = fit)
  }
  out
}
libs <- c(learner1$names[1], glmB$names, "SL.glm_all", "SL.glmnet",
          "SL.randomForest", "SL.gamboost", "SL.mean")
if (custom_models) {
  # User-selected learner names, written by app.py to selected_models.txt.
  eval(parse(text = readLines("selected_models.txt")))
  libnames <- "_custom_models_"
}

# For infold/both: wrap every learner with the rlm row-filter decorator and
# register wrapped versions in the global env under "<name>.olf" (outlier-filtered).
if (outlier_method %in% c("infold", "both")) {
  wrapped_libs <- character(0)
  for (.lname in libs) {
    .wname <- paste0(.lname, ".olf")
    .base  <- tryCatch(get(.lname), error = function(e) NULL)
    if (!is.null(.base) && is.function(.base)) {
      assign(.wname, .make_outlier_wrapper(.base), envir = globalenv())
      wrapped_libs <- c(wrapped_libs, .wname)
    } else {
      wrapped_libs <- c(wrapped_libs, .lname)   # keep unwrapped if not resolvable
    }
  }
  libs <- wrapped_libs
  cat("Wrapped", length(libs), "learners with in-fold rlm filter.\n")
}

#' Smoke mode: swap in two instant, always-fit learners so the end-to-end run
#' finishes in seconds while still exercising every pipeline stage.
if (smoke_test) {
  libs <- c("SL.mean", "SL.lm")
  cat("SMOKE TEST: learner library reduced to", paste(libs, collapse = ", "), "\n")
}

plotnames <- "correlation_plot"
TrainData <- TrainingData

#' Per-learner accumulators (one column per learner). Reset per repetition inside
#' run_one_loop; these top-level copies seed the column names.
Algo_coeff <- as.data.frame(matrix(nrow = 1, ncol = length(libs)))
colnames(Algo_coeff) <- libs
Algo_risk  <- as.data.frame(matrix(nrow = 1, ncol = length(libs)))
colnames(Algo_risk)  <- libs

cat("CV repetitions (loop):", loop, "\n")


# ---- Cross-validated fit ----------------------------------------------------
#' Final SuperLearner-input assertions: guard the design matrix / response the CV
#' loop repeatedly slices, and verify the MC error frame covers every training
#' row (an unmatched key yields an all-NA error row -> NaN MC predictions).
stopifnot(
  "Predictors contains NA/NaN"             = !anyNA(Predictors),
  "Response contains NA/NaN"               = !anyNA(Response),
  "Predictors/Response row count mismatch" = nrow(Predictors) == length(Response),
  "grb_errors missing some training rows"  = all(rownames(TrainingData) %in% rownames(grb_errors))
)
cat("Assertion OK: SuperLearner inputs valid;", nrow(Predictors), "rows,", ncol(Predictors), "predictors\n")

#' Custom learners must exist in THIS (parent) process: the forked mclapply
#' workers inherit the parent global env, but the parent also needs them for the
#' final full-data fit. (options(error = recover) cannot prompt under Rscript --
#' enable it only for interactive debugging.)
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

#' Run one independent CV repetition.
#'
#' Builds its own random k-fold split of the training data, fits SuperLearner per
#' fold (ensemble-weight CV with V = 5), MC-perturbs each fold's test points
#' through the fitted ensemble, and returns a self-contained result list. Reads
#' the globals `TrainData`, `PredictionData`, `libs`, `grb_errors`, `smoke_test`,
#' `one_fold`; mutates nothing global (the Algo_* accumulators are repetition-local
#' and intentionally shadow the top-level copies).
#'
#' @param j integer: repetition index (used only for logging).
#' @return list of 9 elements, by position:
#'   [[1]] numeric vector -- mean holdout prediction per PredictionData row (over folds)
#'   [[2]] numeric        -- mean ensemble coefficient per learner
#'   [[3]] numeric vector -- out-of-fold CV predictions (log10z) per training row
#'   [[4]] numeric        -- mean CV risk per learner
#'   [[5]] numeric vector -- per-fold test correlation (log10z)
#'   [[6]] character      -- suspect-GRB placeholder (unused)
#'   [[7]] numeric vector -- mgcv CV prediction placeholder (zeros, unused)
#'   [[8]] numeric vector -- per-fold test RMSE (log10z)
#'   [[9]] list           -- named GRB id -> numeric vector of MC-perturbed log10z preds
run_one_loop <- function(j) {
  cat(sprintf("[rep %d] starting\n", j))

  # Custom learners again inside the worker (belt-and-suspenders for the fork).
  source("Custom_SL/sl_mgcv_gam.R")
  source("Custom_SL/sl_custom_glm.R")
  source("Custom_SL/sl_custom_bayesglm.R")
  source("Custom_SL/sl_xgboost_safe.R")

  # Repetition-local accumulators (must not share mutable state across workers).
  Algo_coeff <- as.data.frame(matrix(nrow = 1, ncol = length(libs)))
  colnames(Algo_coeff) <- libs
  Algo_risk  <- as.data.frame(matrix(nrow = 1, ncol = length(libs)))
  colnames(Algo_risk)  <- libs

  # Strip response columns from the design matrix; keep the linear-z target.
  responses          <- c("Redshift_crosscheck", "invz", "log10z")
  all_data_scale_wo  <- subset(TrainData, select = !(colnames(TrainData) %in% responses))
  PredictionData_j   <- PredictionData[, colnames(all_data_scale_wo)]  # avoid mutating global
  allZ_wo            <- TrainData$Redshift_crosscheck
  invallZ_wo         <- TrainData$log10z

  # Per-frame design matrices for this repetition (frame 1 == all_data_scale_wo).
  X_frames <- lapply(design_frames, function(D) D[rownames(TrainData), colnames(all_data_scale_wo), drop = FALSE])
  P_frames <- lapply(design_frames, function(D) D[rownames(PredictionData), colnames(all_data_scale_wo), drop = FALSE])

  nwo            <- nrow(all_data_scale_wo)
  results_cv     <- data.frame(Predicted = numeric(nwo), Observed = numeric(nwo))
  results_cv_log <- data.frame(Predicted = numeric(nwo), Observed = numeric(nwo))
  results_cv_mgcv <- data.frame(Predicted = numeric(nwo), Observed = numeric(nwo))

  k_folds <- if (smoke_test) 5 else 10    # integer: outer CV folds
  mc_n    <- if (smoke_test) 20 else 100  # integer: MC draws per test point
  folds   <- createFolds(allZ_wo, k = k_folds)

  CVpred               <- matrix(0, nrow = nrow(PredictionData_j), ncol = 0)
  correlation_log_test <- numeric()
  rmse_log_test        <- numeric()
  sus_GRBs             <- character()
  mc_by_grb            <- list()  # GRB id -> n-vector of MC-perturbed predicted log10(z+1)

  fold_iters <- if (one_fold) 1 else seq_along(folds)
  if (one_fold) cat("ONE_FOLD: running only fold 1 of", length(folds), "\n")
  for (i in fold_iters) {
    cat(sprintf("[rep %d] Iteration %d out of %d\n", j, i, length(folds)))

    train_set <- as.data.frame(all_data_scale_wo[-folds[[i]], ])
    test_set  <- as.data.frame(all_data_scale_wo[ folds[[i]], ])
    invZtrain <- invallZ_wo[-folds[[i]]]
    invZtest  <- invallZ_wo[ folds[[i]]]

    # Per-fold validity check on the exact matrices handed to SuperLearner.
    stopifnot(
      "fold train_set has NA/NaN"      = !anyNA(train_set),
      "fold test_set has NA/NaN"       = !anyNA(test_set),
      "fold invZtrain has NA/NaN"      = !anyNA(invZtrain),
      "fold invZtest has NA/NaN"       = !anyNA(invZtest),
      "fold X/Y row mismatch"          = nrow(train_set) == length(invZtrain),
      "train/test column set mismatch" = identical(colnames(train_set), colnames(test_set))
    )

    # In-fold robust filter: drop training rows an rlm fit flags (>3*MAD residual).
    keep_tr <- rep(TRUE, nrow(train_set))
    rlm_fit <- tryCatch(MASS::rlm(Y ~ ., method = "M", maxit = 50,
                                  data = data.frame(Y = invZtrain, train_set[, head(lassovar, 7)])),
                        error = function(e) NULL)
    if (!is.null(rlm_fit)) {
      rres <- residuals(rlm_fit)
      k    <- abs(rres) <= 3 * mad(rres)
      if (sum(k) >= 30) keep_tr <- k
    }

    # For soft-weights: slice the pre-computed rlm weights to this fold's training rows.
    fold_obs_weights <- if (!is.null(sl_obs_weights)) {
      w <- sl_obs_weights[rownames(train_set)][keep_tr]
      w / mean(w)   # renormalise so mean == 1 (SuperLearner convention)
    } else {
      rep(1, sum(keep_tr))
    }

    # Fit once per imputation/projection frame; the fold prediction is the average.
    fits <- vector("list", length(X_frames))
    prm  <- matrix(NA_real_, nrow(test_set), length(X_frames))
    for (f in seq_along(X_frames)) {
      suppressMessages(capture.output(
        fits[[f]] <- SuperLearner(
          Y          = invZtrain[keep_tr],
          X          = X_frames[[f]][-folds[[i]], , drop = FALSE][keep_tr, , drop = FALSE],
          family     = gaussian(),
          newX       = X_frames[[f]][ folds[[i]], , drop = FALSE],
          SL.library = libs,
          obsWeights = fold_obs_weights,
          method     = method.NNRidge(),
          cvControl  = list(V = 10),  # ensemble-weight CV folds
          verbose    = FALSE
        ),
        file = nullfile()
      ))
      prm[, f] <- as.numeric(fits[[f]]$SL.predict)
    }
    s9 <- fits[[1]]   # frame-1 fit drives MC propagation and weight accounting
    pr <- matrix(rowMeans(prm), ncol = 1)

    #' MC error propagation for this fold's test points: use the fold-trained s9
    #' (so MC respects the out-of-sample regime) and grb_errors keyed the same
    #' way as rownames(test_set). Every test GRB must resolve to a real error row.
    test_grbs   <- rownames(test_set)
    test_errors <- grb_errors[test_grbs, , drop = FALSE]
    stopifnot(
      "test_errors row mismatch"                = nrow(test_errors) == length(test_grbs),
      "test_errors has NA (unmatched GRB key?)" = !anyNA(test_errors)
    )
    mc_matrix <- mc_predict(s9, test_set, test_errors, n = mc_n)
    for (g in test_grbs) {
      mc_by_grb[[g]] <- as.numeric(mc_matrix[g, ])
    }

    results_cv$Predicted[folds[[i]]] <- pr[, 1]
    results_cv$Observed[folds[[i]]]  <- invZtest

    correlation_log_test <- c(correlation_log_test, cor(pr[, 1], invZtest))
    rmse_log_test        <- c(rmse_log_test, sqrt(mean((pr[, 1] - invZtest)^2)))

    Algo_coeff <- rbind(Algo_coeff, coef(s9))
    Algo_risk  <- rbind(Algo_risk,  s9$cvRisk)
    hp <- rowMeans(sapply(seq_along(fits), function(f) as.numeric(predict(fits[[f]], P_frames[[f]])$pred)))
    CVpred     <- cbind(hp, CVpred)
    gc()
  }

  Algo_coeff <- na.omit(Algo_coeff)
  Algo_risk  <- na.omit(Algo_risk)

  list(
    rowMeans(CVpred),
    colMeans(Algo_coeff),
    results_cv$Predicted,
    colMeans(Algo_risk),
    correlation_log_test,
    sus_GRBs,
    results_cv_mgcv$Predicted,
    rmse_log_test,
    mc_by_grb
  )
}

#' Run the repetitions in parallel. Fork count is capped at the physical core
#' count and at `loop`, leaving spare HT threads so internally-threaded learners
#' (e.g. ranger) don't oversubscribe. L'Ecuyer-CMRG + mclapply per-fork seeding
#' give each repetition an independent RNG stream, so fold splits and MC draws
#' are never duplicated across repetitions.
n_cores <- max(1L, min(loop, parallel::detectCores(logical = FALSE)))
cat("Parallelizing", loop, "CV repetition(s) across", n_cores, "core(s)\n")
RNGkind("L'Ecuyer-CMRG")
set.seed(20240603)
CVmodel <- parallel::mclapply(seq_len(loop), run_one_loop, mc.cores = n_cores)

#' mclapply never aborts on a worker error -- it returns a try-error / NULL in
#' that slot. Fail loudly so a half-populated CVmodel can't feed the plots.
.failed <- !vapply(CVmodel, is.list, logical(1))
if (any(.failed)) {
  stop("Parallel CV failed in ", sum(.failed), " of ", loop,
       " repetition(s). First failure:\n",
       paste(utils::capture.output(print(CVmodel[[which(.failed)[1]]])), collapse = "\n"))
}


# ---- Aggregate per-repetition CV metrics ------------------------------------
z_e <- 1  # numeric: z + 1 offset used in the log10(z+1) transform

correl       <- as.vector(0)
LinearCorrel <- as.vector(0)
preds        <- matrix(nrow = nrow(Predictors), ncol = loop)  # per-rep CV predictions, averaged later
co           <- matrix(nrow = loop, ncol = length(libs))      # per-rep ensemble coefficients
AlgoRisks    <- matrix(nrow = loop, ncol = length(libs))      # per-rep CV risks
linearrms    <- as.vector(0)

for (j in 1:loop) {
  preds[, j]      <- CVmodel[[j]][[3]]
  co[j, ]         <- CVmodel[[j]][[2]]
  AlgoRisks[j, ]  <- CVmodel[[j]][[4]]
  correl[j]       <- cor(preds[, j], Response)
  LinearCorrel[j] <- cor(10^preds[, j] - z_e, TrainingData$Redshift_crosscheck)
  linearrms[j]    <- sqrt(mean((TrainingData$Redshift_crosscheck - (10^preds[, j] - 1))^2))
}
cat(sprintf("CV r (all training GRBs): log10(z+1)=%.4f linear-z=%.4f\n",
            cor(rowMeans(preds), Response),
            cor(10^rowMeans(preds) - 1, TrainingData$Redshift_crosscheck)))


# ---- Catastrophic-outlier removal (retrained second pass) --------------------
# Drop training GRBs whose pooled CV residual exceeds 2 sigma, then re-run the CV
# so all reported metrics come from retrained models -- never from post-hoc
# filtering of the same predictions. Disable with REMOVE_CAT_OUTLIERS=false.
remove_cat_outliers <- tolower(Sys.getenv("REMOVE_CAT_OUTLIERS", "true")) == "true"
if (remove_cat_outliers && !one_fold) {
  cv_resid <- rowMeans(preds) - Response
  keep_cat <- abs(cv_resid) <= 2 * sd(cv_resid)
  if (any(!keep_cat)) {
    cat("Catastrophic outliers removed (2-sigma), re-running CV:", sum(!keep_cat), "GRBs:",
        paste(rownames(TrainingData)[!keep_cat], collapse = ", "), "\n")
    TrainingData <- TrainingData[keep_cat, ]
    TrainData    <- TrainingData
    Response     <- TrainingData$log10z
    Predictors   <- subset(TrainingData, select = -c(log10z, Redshift_crosscheck))
    CVmodel  <- parallel::mclapply(seq_len(loop), run_one_loop, mc.cores = n_cores)
    .failed  <- !vapply(CVmodel, is.list, logical(1))
    if (any(.failed)) stop("Parallel CV failed in ", sum(.failed), " repetition(s) after outlier removal.")
    preds     <- matrix(nrow = nrow(Predictors), ncol = loop)
    co        <- matrix(nrow = loop, ncol = length(libs))
    AlgoRisks <- matrix(nrow = loop, ncol = length(libs))
    for (j in 1:loop) {
      preds[, j]      <- CVmodel[[j]][[3]]
      co[j, ]         <- CVmodel[[j]][[2]]
      AlgoRisks[j, ]  <- CVmodel[[j]][[4]]
      correl[j]       <- cor(preds[, j], Response)
      LinearCorrel[j] <- cor(10^preds[, j] - z_e, TrainingData$Redshift_crosscheck)
      linearrms[j]    <- sqrt(mean((TrainingData$Redshift_crosscheck - (10^preds[, j] - 1))^2))
    }
    cat(sprintf("CV r after outlier removal: log10(z+1)=%.4f linear-z=%.4f\n",
                cor(rowMeans(preds), Response),
                cor(10^rowMeans(preds) - 1, TrainingData$Redshift_crosscheck)))
  }
}


# ---- MC post-processing -----------------------------------------------------
#' Flatten MC samples across all repetitions. Each GRB appears in exactly one
#' fold per repetition, so its final MC vector has length loop * mc_n.
cat("Flattening MC samples across", loop, "repetitions...\n")
mc_raw <- list()
for (j in 1:loop) {
  rep_mc <- CVmodel[[j]][[9]]
  for (g in names(rep_mc)) {
    mc_raw[[g]] <- c(mc_raw[[g]], rep_mc[[g]])
  }
}

mc_grbs       <- names(mc_raw)
mc_matrix_all <- do.call(rbind, lapply(mc_grbs, function(g) mc_raw[[g]]))
rownames(mc_matrix_all) <- mc_grbs
cat("MC matrix:", nrow(mc_matrix_all), "rows x", ncol(mc_matrix_all), "samples\n")

mc_summary <- summarize_mc(mc_matrix_all, level = 0.68)

# Align truths to the MC summary's row identity.
y_true_named <- setNames(Response, rownames(TrainingData))                       # numeric: true log10(z+1)
z_true_named <- setNames(TrainingData$Redshift_crosscheck, rownames(TrainingData))  # numeric: true z

inside_flag <- inside_mc_cone(mc_summary, y_true_named)
cat(sprintf("MC coverage: %.1f%% of rows have truth inside 95%% MC interval\n", 100 * mean(inside_flag)))

#' All-GRB vs inside-cone metrics (predicted = MC mean, matching the plots),
#' reported in both log10(z+1) and linear z so we can see how much error the
#' GRBs whose truth falls outside their 95% MC band carry.
.mc_metrics <- function(obs, pred) {
  ok  <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  c(n = length(obs), r = cor(obs, pred), rmse = sqrt(mean((obs - pred)^2)))
}
.yt <- y_true_named[rownames(mc_summary)]
.zt <- z_true_named[rownames(mc_summary)]
for (.scope in c("all", "inside")) {
  .sel  <- if (.scope == "all") rep(TRUE, nrow(mc_summary)) else inside_flag
  .mlog <- .mc_metrics(.yt[.sel], mc_summary$log10z_mean[.sel])
  .mlin <- .mc_metrics(.zt[.sel], mc_summary$z_mean[.sel])
  cat(sprintf("MC metrics [%-6s] log10(z+1): n=%d r=%.3f RMSE=%.3f | z: r=%.3f RMSE=%.3f\n",
              .scope, .mlog["n"], .mlog["r"], .mlog["rmse"], .mlin["r"], .mlin["rmse"]))
}
rm(.mc_metrics, .yt, .zt, .scope, .sel, .mlog, .mlin)


# ---- Write outputs & plots --------------------------------------------------
#' Map row identity -> GRB name via raw_xray_data (same rownames); GRBPred no
#' longer carries a GRB column after the predictor rebuild.
grb_name_lookup <- setNames(raw_xray_data$GRB, rownames(raw_xray_data))
cv_pred_mean    <- setNames(rowMeans(preds), rownames(TrainingData))

mc_output <- data.frame(
  row_id         = rownames(mc_summary),
  GRB            = grb_name_lookup[rownames(mc_summary)],
  y_true         = y_true_named[rownames(mc_summary)],
  z_true         = z_true_named[rownames(mc_summary)],
  y_pred_cv_mean = cv_pred_mean[rownames(mc_summary)],
  log10z_mc_mean = mc_summary$log10z_mean,
  log10z_mc_sd   = mc_summary$log10z_sd,
  log10z_lower   = mc_summary$log10z_lower,
  log10z_upper   = mc_summary$log10z_upper,
  z_mc_mean      = mc_summary$z_mean,
  z_mc_sd        = mc_summary$z_sd,
  z_lower        = mc_summary$z_lower,
  z_upper        = mc_summary$z_upper,
  inside_cone    = inside_flag,
  row.names      = NULL,
  stringsAsFactors = FALSE
)
write.csv(mc_output, file.path(out_files_dir, "mc_predictions.csv"), row.names = FALSE)
saveRDS(mc_raw, file.path(out_files_dir, "mc_raw.rds"))
cat("Wrote OutputFiles/mc_predictions.csv and OutputFiles/mc_raw.rds\n")

make_mc_plots(mc_summary, y_true_named, z_true_named, out_dir = PLOTaddr)
cat("Wrote MC plots to", PLOTaddr, "\n")

# Ensemble-weight bar charts (sorted to PNG; unsorted to the default device).
png(filename = paste(PLOTaddr, "model_compare_plot.png"), res = 500, width = 3000, height = 3000)
par(mar = c(5, 10, 4, 2))
barplot(sort(colMeans(co)), names.arg = libs[order(colMeans(co))], horiz = TRUE, las = 1,
        cex.main = 0.95,
        main = "SuperLearner ensemble weights\n(mean over CV repetitions)",
        xlab = "Mean ensemble coefficient (larger = more influence)")
dev.off()

barplot(colMeans(co), names.arg = libs, horiz = TRUE, las = 1,
        cex.main = 0.95,
        main = "SuperLearner ensemble weights (unsorted)",
        xlab = "Mean ensemble coefficient (larger = more influence)")
par(mar = c(5, 4, 4, 2))

# Correlation plot over ALL training GRBs (with catastrophic outliers).
plotnames <- paste0("_with_catOutl_", "correlation_plot")
results <- result_plotter(rownames(TrainingData), rowMeans(preds), Response,
                          apply(preds, 1, max), apply(preds, 1, min))

# Persist the final model trained on the full training set.
InsideCone <- read.csv(paste(addr, "Results_wo_catout", plotnames, ".csv", sep = ""), row.names = 1)
sl_model   <- SuperLearner(Y = Response, X = Predictors, family = gaussian(),
                           SL.library = libs, method = method.NNRidge(),
                           cvControl = list(V = 10), verbose = FALSE)
saveRDS(sl_model, file = file.path(out_dir, "superlearner_model"))

# Correlation plot for data inside the 2-sigma cone (catastrophic outliers removed).
plotnames   <- paste0("_without_catOutl_", "correlation_plot")
Good_results <- result_plotter(rownames(InsideCone), InsideCone$InvZphot, InsideCone$InvZspec,
                               InsideCone$pred_max, InsideCone$pred_min)
