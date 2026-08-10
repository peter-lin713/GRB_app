#' Prep script: take one of the standard 10-feature (+is_optical) combination
#' datasets, DROP is_optical entirely (so stock superlearner.R's automatic
#' Daume domain-adaptation block -- `if ("is_optical" %in% colnames(...))` --
#' never fires), then run the standard non-negotiable-outlier removal + MICE +
#' M-estimator cut using the STANDING (reused, well-generalized) formula
#' directly, no fresh formula search. Mirrors the earlier "Theil-Sen WITHOUT
#' is_optical" ablation, generalized to any of the same-schema datasets.
#'
#' Usage: Rscript prep_no_isoptical.R <input_csv> <method_tag> <out_dir>
require(mice)
require(MASS)
require(dplyr)

STANDING_FORMULA <- "Response ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)^2 + log10NH + PhotonIndex + log10Ta + Alpha + log10NHSqr + log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr"

args       <- commandArgs(trailingOnly = TRUE)
input_csv  <- normalizePath(args[1])
method_tag <- args[2]
out_dir    <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

raw_xray_data <- read.csv(input_csv, header = TRUE, row.names = 1)
# is_optical stays in raw_xray_data through the non-negotiable-outlier check
# below (it uses is_optical for its own diagnostic logging) and gets dropped
# right after -- never enters features_for_mice_preds, SqrData, or the final
# cut, so stock superlearner.R's Daume block never fires downstream.

if ("T90" %in% colnames(raw_xray_data)) {
  raw_xray_data <- raw_xray_data[raw_xray_data$T90 > 2, ]
  raw_xray_data$log10T90 <- log10(raw_xray_data$T90)
} else {
  raw_xray_data <- raw_xray_data[raw_xray_data$log10T90 > log10(2), ]
}

features_for_mice_preds <- subset(raw_xray_data, select = c(log10T90, log10Fa, log10Ta,
                                                             Alpha, Beta, Gamma,
                                                             log10Fluence, PhotonIndex,
                                                             log10NH, log10PeakFlux))

source("../../generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = ".", method = method_tag)
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

raw_xray_data$is_optical <- NULL   # the whole point of this script -- drop now,
                                    # after the outlier check used it, before
                                    # anything downstream (MICE, M-estimator,
                                    # final_cut) can see it.

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
                                             fluence_err_col, "PhotonIndexErr", peakflux_err_col)]

set.seed(1)
mice_model_preds <- mice(data = features_for_mice_preds, m = 20, maxit = 20,
                         method = 'midastouch', printFlag = FALSE)
features_for_mice_preds <- complete(mice_model_preds, 20)

mice_model_errs <- mice(data = features_for_mice_errs, m = 20, maxit = 20,
                        method = 'midastouch', printFlag = FALSE)
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
