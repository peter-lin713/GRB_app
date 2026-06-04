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
#'   2. Clean & filter       -- drop short / instrumental GRBs (log10T90 cut) and
#'                              null out physically-impossible feature values.
#'   3. Impute (MICE)        -- multiply-impute missing predictors and (dex-scaled)
#'                              measurement errors with `mice` (midastouch, m = 20).
#'   4. Feature engineering  -- order predictors by LASSO importance; append a
#'                              squared term per predictor (SqrTermGen).
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
#'   Env vars SMOKE_TEST / SMOKE_N / ONE_FOLD enable fast structural runs.
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
run_locally <- FALSE  # TRUE only for interactive debugging
if (run_locally) {
  raw_xray_data    <- read.csv("combined_data_with_redshift_V8.csv", header = TRUE, row.names = 1)
  do_mice          <- TRUE
  upsampling       <- FALSE
  do_m_estimator   <- TRUE
  custom_models    <- FALSE
  weight_threshold <- 0.65
  loop             <- 5
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

  stopifnot(!anyDuplicated(raw_xray_data$GRB))
  rownames(raw_xray_data) <- raw_xray_data$GRB
}
#' Do NOT re-read input_file after this point: the read above already loaded
#' raw_xray_data and set rownames <- GRB. A second read would crash in
#' run_locally mode (input_file undefined) and would reset rownames to 1..N,
#' discarding the GRB-keyed identity that grb_errors, the MC lookups and the GRB
#' name lookup all depend on.


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
addr     <- "Results/"       # character: directory for result CSVs
PLOTaddr <- "Plot_Output/"   # character: directory for plot PNGs
if (!dir.exists(PLOTaddr)) dir.create(PLOTaddr)
if (!dir.exists(addr))     dir.create(addr)


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


# ---- Impute (MICE) ----------------------------------------------------------
#' With do_mice = TRUE, multiply-impute predictors and (dex) errors separately
#' (midastouch, m = 20, taking the 20th completed set) and combine into GRBPred.
#' With do_mice = FALSE, drop every incomplete row instead. `GRBPred` is the
#' modeling frame from here on: a data.frame of imputed features + errors keyed
#' by GRB id.
if (do_mice) {
  set.seed(1)
  png(filename = "MICE_missing_features.png", width = 1000, height = 1000, res = 200)
  md.pattern(features_for_mice_preds, rotate.names = TRUE)
  title(main = "Missing-data pattern before MICE imputation",
        sub  = "Rows = observed patterns (left count = # GRBs); blue = observed, red = missing; right/bottom = # missing")
  dev.off()

  mice_model_preds        <- mice(data = features_for_mice_preds, m = 20, method = "midastouch", printFlag = FALSE)
  features_for_mice_preds <- complete(mice_model_preds, 20)

  mice_model_errs        <- mice(data = features_for_mice_errs, m = 20, method = "midastouch", printFlag = FALSE)
  features_for_mice_errs <- complete(mice_model_errs, 20)

  GRBPred <- cbind(features_for_mice_preds, features_for_mice_errs)
} else {
  GRBPred <- na.omit(cbind(features_for_mice_preds, features_for_mice_errs))
}

#' Assertion: GRBPred must be complete after imputation. anyNA() catches NaN too.
#' If this fires with do_mice = TRUE, MICE failed to fill some cell (e.g. an
#' all-NA column); with do_mice = FALSE, na.omit() did not drop every gap.
stopifnot("GRBPred still contains NA/NaN after imputation step" = !anyNA(GRBPred))
cat("Assertion OK: no NA in GRBPred after imputation (", nrow(GRBPred), "rows x", ncol(GRBPred), "cols )\n")

GRBPred$GRB <- raw_xray_data$GRB

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

if (!dir.exists("OutputFiles")) dir.create("OutputFiles")
#' Snapshot the imputed/engineered frame; m_estimator.R reads it back.
if (do_mice) {
  write.csv(GRBPred, "OutputFiles/grb_xray_imputed.csv")
} else {
  write.csv(GRBPred, "OutputFiles/grb_xray.csv")
}


# ---- Outlier handling (optional) --------------------------------------------
if (upsampling)      source("upsampling.R")     # rebalances GRBPred by redshift
if (do_m_estimator)  source("m_estimator.R")    # robust-regression catastrophic-outlier cut

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

#' Re-derive the modeling matrix as (linear LASSO predictors + their squares) and
#' attach the responses. lassovar carries only linear names, so SqrTermGen
#' re-derives the squares cleanly (no "...SqrSqr").
O1Predictors <- subset(GRBPred, select = lassovar)
O2Predictors <- SqrTermGen(O1Predictors)
GRBPred      <- cbind(O2Predictors, Responses)

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
source("Custom_SL/sl_custom_bayesglm.R")

formula_table_GAM <- read.table("Best_formula_GAM.txt")
bestGAM1          <- apply(as.matrix(formula_table_GAM[, 2]), 1, as.formula)  # list of formula
learner1          <- create.Learner("SL.mgcv_gam", tune = list(gam.model = c(bestGAM1)),
                                     detailed_names = FALSE, name_prefix = "gam")

formula_table_GLM <- read.table("Best_formula_GLM.txt")
best_lm3          <- apply(as.matrix(formula_table_GLM[, 2]), 1, as.formula)  # list of formula
sl_glm1           <- create.Learner("SL.custom_glm", tune = list(glm.model = c(best_lm3)),
                                     detailed_names = FALSE, name_prefix = "cglm")

# bayesglm reuses the GLM formula set; needs the 'arm' package.
sl_bglm <- create.Learner("SL.custom_bayesglm", tune = list(bglm.model = c(best_lm3)),
                          detailed_names = FALSE, name_prefix = "bglm")

custom_learner_names <- c(learner1$names, sl_glm1$names, sl_bglm$names)  # character: tuned learner names
libnames             <- "_OG_ALL_CUSTOM_"

# caret random forest learner.
tune_caret    <- list(method = c("rf"), tuneLength = 1, verboseIter = FALSE)
caret_learner <- create.Learner("SL.caret", tune = tune_caret, detailed_names = TRUE)

analyze_all <- FALSE
#' analyze_all / custom_models set `libs` and `libnames`, but BOTH are overwritten
#' by the unconditional `libs <- ...` assignment below; they are kept for the
#' app's selected-models path and historical "_ALL_" runs. The active library is
#' always the one assigned at the end of this block.
if (analyze_all) {
  libs <- c(
    "SL.caret.rpart", "SL.earth", "SL.ipredbagg",
    "SL.mean", "SL.nnet", "SL.randomForest", "SL.ranger",
    "SL.rpart", "SL.step", "SL.step.forward",
    "SL.step.interaction", "SL.stepAIC", "SL.xgboost"
  )
  libnames <- "_ALL_"
}
if (custom_models) {
  # User-selected learner names, written by app.py to selected_models.txt.
  libs_line <- readLines("selected_models.txt")
  eval(parse(text = libs_line))
  libnames <- "_custom_models_"
}

#' Generic SuperLearner learners. NOTE: the original pipeline appended the custom
#' learners via `create.Learner(...)$library`, but that field is NULL on
#' create.Learner objects (the names live in `$names`), and `c(x, NULL)` silently
#' dropped them -- so the original results came from the generic learners alone.
#' We add the custom learners explicitly via `$names`.
generic_libs <- c(
  "SL.glmnet", "SL.xgboost_safe",
  "SL.caret.rpart", "SL.earth", "SL.ipredbagg",
  "SL.mean", "SL.nnet", "SL.randomForest", "SL.ranger",
  "SL.rpart", "SL.step", "SL.step.forward",
  "SL.step.interaction", "SL.stepAIC"
)
#' Active learner library: custom GAM + custom GLM + bayesglm + caret RF (by
#' $names) plus the generics. character vector.
libs <- c(custom_learner_names, caret_learner$names, generic_libs)

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

all_lasso_vars <- character()
gam_vars <- c(
  "log10FaSqr", "log10Fa", "log10PeakFlux",
  "log10NHSqr", "log10NH",
  "PhotonIndex", "PhotonIndexSqr",
  "log10Ta", "log10TaSqr",
  "Gamma", "GammaSqr",
  "Alpha", "AlphaSqr"
)
cat("The number of NAs in GRBPred is:", sum(is.na(GRBPred)), "\n")
cat("GAM vars missing from GRBPred:", paste(setdiff(gam_vars, colnames(GRBPred)), collapse = ", "), "\n")
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
source("Custom_SL/sl_custom_bayesglm.R")
source("Custom_SL/sl_xgboost_safe.R")

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

    suppressMessages(capture.output(
      s9 <- SuperLearner(
        Y          = invZtrain,
        X          = train_set,
        family     = gaussian(),
        newX       = test_set,
        SL.library = libs,
        cvControl  = list(V = 5),  # ensemble-weight CV folds
        verbose    = FALSE
      ),
      file = nullfile()
    ))
    pr <- s9$SL.predict

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
    CVpred     <- cbind(predict(s9, PredictionData_j)$pred, CVpred)
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

mc_summary <- summarize_mc(mc_matrix_all, level = 0.95)

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
write.csv(mc_output, "OutputFiles/mc_predictions.csv", row.names = FALSE)
saveRDS(mc_raw, "OutputFiles/mc_raw.rds")
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
                           SL.library = libs, cvControl = list(V = 5), verbose = FALSE)
saveRDS(sl_model, file = "superlearner_model")

# Correlation plot for data inside the 2-sigma cone (catastrophic outliers removed).
plotnames   <- paste0("_without_catOutl_", "correlation_plot")
Good_results <- result_plotter(rownames(InsideCone), InsideCone$InvZphot, InsideCone$InvZspec,
                               InsideCone$pred_max, InsideCone$pred_min)
