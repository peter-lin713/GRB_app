#!/usr/bin/env Rscript
#' build_native_opt_dataset.R -- ADDS the real, native optical-domain
#' measurements (Alpha_opt/Beta_opt/log10Faopt/log10Taopt from the optical
#' catalog, renamed *_native here) as four EXTRA candidate columns, alongside
#' -- not replacing -- the existing emcee-projected Alpha/Beta/log10Fa/log10Ta.
#' Only the 73 optically-sourced GRBs (is_optical == 1) have a real value for
#' these; the 223 X-ray-native GRBs get NA, imputed downstream by MICE the
#' same way any other missing predictor is (75% missingness on this column,
#' so treat that imputation as weak/exploratory, not a strong signal).
#'
#' This lets LASSO's own ranking decide whether the native-optical version of
#' a variable carries more predictive signal than the standard (emcee-
#' calibrated) version for the GRBs that have both -- rather than us
#' presupposing an answer by overwriting one with the other.

setwd("/Users/petelin/Desktop/test/GRB_app")

d   <- read.csv("Data/superlearner_training_emcee_v2_recovered.csv", row.names = 1)
opt <- read.csv("Data/OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv", stringsAsFactors = FALSE)

norm_grb <- function(g) {
  g <- gsub("GRB", "", as.character(g))
  ifelse(grepl("[A-Za-z]$", g), g, paste0(g, "A"))
}

opt_rows <- which(d$is_optical == 1)
m  <- match(norm_grb(rownames(d)[opt_rows]), norm_grb(opt$GRB))
ok <- !is.na(m)
cat("is_optical rows:", length(opt_rows), "| matched to optical catalog:", sum(ok), "\n")
opt_rows <- opt_rows[ok]; m <- m[ok]

d$log10Fa_native   <- NA_real_
d$log10Ta_native   <- NA_real_
d$Alpha_native     <- NA_real_
d$Beta_native      <- NA_real_
d$log10FaErr_native <- NA_real_
d$log10TaErr_native <- NA_real_
d$AlphaErr_native   <- NA_real_
d$BetaErr_native    <- NA_real_

d$log10Fa_native[opt_rows]    <- opt$log10Faopt[m]
d$log10Ta_native[opt_rows]    <- opt$log10Taopt[m]
d$Alpha_native[opt_rows]      <- opt$Alpha_opt[m]
d$Beta_native[opt_rows]       <- opt$Beta_opt[m]
d$log10FaErr_native[opt_rows] <- opt$log10FaErr_opt[m]
d$log10TaErr_native[opt_rows] <- opt$logTaErr_opt[m]
d$AlphaErr_native[opt_rows]   <- opt$AlphaErr_opt[m]
d$BetaErr_native[opt_rows]    <- opt$betaErr_opt[m]

write.csv(d, "Data/superlearner_training_emcee_v2_native_opt.csv")
cat("Wrote Data/superlearner_training_emcee_v2_native_opt.csv (", nrow(d), "GRBs x", ncol(d),
    "cols; 4 native-optical columns added, non-NA for", length(opt_rows), "of", nrow(d), "rows)\n")
