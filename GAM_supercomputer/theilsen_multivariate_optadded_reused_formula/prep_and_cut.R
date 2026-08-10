#' Prep script for "Theil-Sen multivariate calibration + opt cols added,
#' reusing the standing (more generalizable) formula" -- no fresh formula
#' search. Mirrors hardzero_domain_v2_formula_generation/GAM_analysis_8variables.R's
#' data-prep and M-estimator-cut logic exactly (same features_for_mice_preds/errs
#' construction, same non-negotiable-outlier removal, same MICE settings), but
#' skips the 100-split combinatorial search entirely and fits the M-estimator
#' directly with the standing formula, since that's the formula this run reuses.
require(mice)
require(MASS)
require(dplyr)

STANDING_FORMULA <- "Response ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)^2 + log10NH + PhotonIndex + log10Ta + Alpha + log10NHSqr + log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr"

raw_xray_data <- read.csv("superlearner_training_multivariate_theilsen_optadded.csv", header = TRUE, row.names = 1)

if ("T90" %in% colnames(raw_xray_data)) {
  raw_xray_data <- raw_xray_data[raw_xray_data$T90 > 2, ]
  raw_xray_data$log10T90 <- log10(raw_xray_data$T90)
} else {
  raw_xray_data <- raw_xray_data[raw_xray_data$log10T90 > log10(2), ]
}

features_for_mice_preds <- subset(raw_xray_data, select = c(log10T90,
                                                             log10Fa,
                                                             log10Ta,
                                                             Alpha,
                                                             Beta,
                                                             log10Fa_opt,
                                                             log10Ta_opt,
                                                             Alpha_opt,
                                                             Beta_opt,
                                                             Gamma,
                                                             log10Fluence,
                                                             PhotonIndex,
                                                             log10NH,
                                                             log10PeakFlux))

source("../../generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = ".", method = "theilsen_multivariate_optadded")
to_drop    <- .outliers$to_drop
reason_for <- setNames(.outliers$candidates$reasons, .outliers$candidates$GRB)

grb_ids      <- rownames(raw_xray_data)
dropped_mask <- grb_ids %in% to_drop
if (any(dropped_mask)) {
  dropped_log <- data.frame(
    GRB    = grb_ids[dropped_mask],
    reason = ifelse(grb_ids[dropped_mask] %in% names(reason_for),
                     reason_for[grb_ids[dropped_mask]], "manual/visual review")
  )
  write.csv(dropped_log, "removed_outliers_final.csv", row.names = FALSE)
  cat("Removed", sum(dropped_mask), "GRB(s) confirmed physically infeasible -- see removed_outliers_final.csv\n")
  raw_xray_data           <- raw_xray_data[!dropped_mask, ]
  features_for_mice_preds <- features_for_mice_preds[!dropped_mask, ]
} else {
  cat("confirmed_outliers_to_drop.txt is empty or matched nothing -- no GRBs removed.\n")
}

features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20]       <- NA
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 2]              <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3]            <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3]            <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

fluence_err_col  <- if ("FluenceErr"  %in% colnames(raw_xray_data)) "FluenceErr"  else "log10FluenceErr"
peakflux_err_col <- if ("PeakFluxErr" %in% colnames(raw_xray_data)) "PeakFluxErr" else "log10PeakFluxErr"

features_for_mice_errs <- raw_xray_data[, c("T90Err", "log10FaErr", "log10TaErr",
                                             "AlphaErr", "BetaErr",
                                             "log10FaErr_opt", "log10TaErr_opt",
                                             "AlphaErr_opt", "BetaErr_opt",
                                             fluence_err_col,
                                             "PhotonIndexErr", peakflux_err_col)]

# seed=1 (the convention used everywhere else) hits a degenerate midastouch
# resample on features_for_mice_preds for this dataset's specific collinearity
# structure (near-perfect correlation among the _opt columns for the ~218 rows
# where they aren't NaN -- rank is still full, 14/14, but some internal
# bootstrap subsample apparently isn't) -- "system is computationally
# singular". Seeds 2/3/42 all complete the full m=20,maxit=20 preds run
# cleanly; using seed=2.
set.seed(2)
mice_model_preds <- mice(data = features_for_mice_preds, m = 20, maxit = 20,
                         method = 'midastouch', printFlag = FALSE)
features_for_mice_preds <- complete(mice_model_preds, 20)

# features_for_mice_errs hits the SAME singularity with midastouch under every
# seed tried (1,2,3,4,5,6,7,8,9,10,42,99,123,555,2024,777 -- not a seed-luck
# issue): the four *Err_opt columns are identically 0 for all 134 X-ray-only
# rows (sd=0 in that complete-case block) while simultaneously NaN for the 78
# shared rows -- combined with T90Err's separate 73-row NaN block (missing
# only for optical-only rows), some internal midastouch donor-regression
# subproblem is left with a rank-deficient design no matter the seed. This
# block of columns is a diagnostic-only artifact anyway: none of the *Err
# columns feed the M-estimator formula below (Formula_for_outlier only
# references value columns), and final_outliers_removed.csv is written from
# raw_xray_data (not from this MICE-completed frame) exactly like every other
# formula-generation script -- so switching just this call to 'pmm' (mice's
# other standard, more numerically robust method) has no effect on the actual
# outlier cut or downstream SuperLearner input, only on the informational
# grb_xray_imputed.csv snapshot.
mice_model_errs <- mice(data = features_for_mice_errs, m = 20, maxit = 20,
                        method = 'pmm', printFlag = FALSE)
features_for_mice_errs <- complete(mice_model_errs, 20)

GRBPred <- cbind(features_for_mice_preds, features_for_mice_errs)

if (!"log10T90Err" %in% colnames(GRBPred)) {
  T90 <- 10^GRBPred$log10T90
  GRBPred$log10T90Err <- GRBPred$T90Err / (T90 * log(10))
}
if (!"log10FluenceErr" %in% colnames(GRBPred)) {
  Fluence <- 10^GRBPred$log10Fluence
  GRBPred$log10FluenceErr <- GRBPred$FluenceErr / (Fluence * log(10))
}
if (!"log10PeakFluxErr" %in% colnames(GRBPred)) {
  PeakFlux <- 10^GRBPred$log10PeakFlux
  GRBPred$log10PeakFluxErr <- GRBPred$PeakFluxErr / (PeakFlux * log(10))
}

GRBPred$Redshift_crosscheck <- raw_xray_data$Redshift_crosscheck
GRBPred$log10z <- log10(raw_xray_data$Redshift_crosscheck + 1)

write.csv(GRBPred, "grb_xray_imputed.csv")

TrainData <- GRBPred

SqrTermGen <- function(inputData) {
  indVar <- colnames(inputData)
  for (i in 1:length(indVar)) {
    for (j in i:(length(indVar))) {
      if (indVar[i] == indVar[j]) {
        inputData[[paste(indVar[i], "Sqr", sep = "")]] <- inputData[, indVar[i]] * inputData[, indVar[j]]
      }
    }
  }
  return(inputData)
}

# All base predictors kept (no LASSO top-7 restriction -- irrelevant here since
# we already know which formula we're fitting; restricting to a data-driven
# top-7 risked silently dropping a variable the standing formula needs).
Predictor <- subset(TrainData, select = colnames(features_for_mice_preds))
SqrData   <- SqrTermGen(Predictor)
SqrData$Response <- TrainData$log10z

Formula_for_outlier <- as.formula(STANDING_FORMULA)
M_est <- MASS::rlm(Formula_for_outlier, data = SqrData, method = "M", maxit = 50)
print(Formula_for_outlier)
writeLines(deparse(Formula_for_outlier), "Formula_for_outlier.txt")

weights <- M_est$w
weight_threshold <- quantile(weights, 0.05)

kept_rows    <- rownames(SqrData)[weights > weight_threshold]
removed_rows <- rownames(SqrData)[weights <= weight_threshold]

cat("M-estimator outlier cut:", length(removed_rows), "of", nrow(SqrData),
    "GRBs removed (weight <=", round(weight_threshold, 4), ")\n")
writeLines(removed_rows, "removed_outliers.txt")

final_cut <- raw_xray_data[rownames(raw_xray_data) %in% kept_rows, ]
write.csv(final_cut, "final_outliers_removed.csv")
cat("Final outlier-removed data saved:", nrow(final_cut), "GRBs -> final_outliers_removed.csv\n")
