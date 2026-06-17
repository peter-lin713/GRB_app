#!/usr/bin/env Rscript
# Trace sample size at each filtering step without running SuperLearner.
suppressPackageStartupMessages({
  library(mice)
})

INPUT <- "Data/superlearner_training_emcee_errcut_relative.csv"

raw <- read.csv(INPUT, header = TRUE, stringsAsFactors = FALSE)
if (!"GRB" %in% names(raw) && "X" %in% names(raw)) raw$GRB <- raw$X
if (!"log10T90" %in% names(raw) && "T90" %in% names(raw)) raw$log10T90 <- log10(raw$T90)
rownames(raw) <- raw$GRB

cat("1. After reading CSV:                    ", nrow(raw), "\n")

# Short GRB filter (keeps NAs for MICE)
short_grbs <- !is.na(raw$log10T90) & raw$log10T90 <= 0.301
cat("   Short GRBs (T90 <= 2s, known):        ", sum(short_grbs), "\n")
raw <- raw[is.na(raw$log10T90) | raw$log10T90 > 0.301, ]
cat("2. After T90 > 2s filter (NAs kept):     ", nrow(raw), "\n")

# Feature outlier nulling
feats <- subset(raw, select = c(log10T90, log10Fa, log10Ta, Alpha, Beta,
                                 Gamma, log10Fluence, PhotonIndex, log10NH, log10PeakFlux))
before_na <- colSums(is.na(feats))
feats$log10PeakFlux[is.infinite(feats$log10PeakFlux)] <- NA
feats$log10NH[feats$log10NH < 20]          <- NA
feats$Beta[feats$Beta > 3]                 <- NA
feats$Gamma[feats$Gamma > 3]               <- NA
feats$Alpha[feats$Alpha > 3]               <- NA
feats$PhotonIndex[feats$PhotonIndex < 0]   <- NA
after_na <- colSums(is.na(feats))
cat("3. Feature value nulling (turned to NA, not removed):\n")
for (col in names(after_na)) {
  added <- after_na[col] - before_na[col]
  if (added > 0) cat("     ", col, ":", added, "values nulled\n")
}
cat("   GRBs with >=1 NA in features:         ",
    sum(apply(feats, 1, anyNA)), "\n")
cat("   Complete cases (no MICE needed):       ",
    sum(complete.cases(feats)), "\n")

# MICE
cat("4. Running MICE (m=20, midastouch)...\n")
set.seed(1)
mice_out <- mice(data = feats, m = 20, method = "midastouch", printFlag = FALSE)
feats_imp <- complete(mice_out, 20)
cat("5. After MICE imputation:                ", nrow(feats_imp), "\n")

# Attach response
feats_imp$log10z        <- log10(raw$Redshift_crosscheck + 1)
feats_imp$Redshift_crosscheck <- raw$Redshift_crosscheck
cat("   Missing response (log10z):             ",
    sum(is.na(feats_imp$log10z)), "\n")
feats_imp <- feats_imp[!is.na(feats_imp$log10z), ]
cat("6. After dropping missing response:      ", nrow(feats_imp), "\n")

# M-estimator (source it)
GRBPred <- feats_imp
GRBPred$log10z <- log10(raw[rownames(GRBPred), "Redshift_crosscheck"] + 1)
TrainingData <- GRBPred
m_est_pct <- 0.05
PLOTaddr <- "runs/diag_plots/"
addr     <- "runs/diag_results/"
out_files_dir <- "runs/diag_output/"
dir.create(PLOTaddr, recursive = TRUE, showWarnings = FALSE)
dir.create(addr,     recursive = TRUE, showWarnings = FALSE)
dir.create(out_files_dir, recursive = TRUE, showWarnings = FALSE)
source("m_estimator.R")
cat("7. After M-estimator (", m_est_pct*100, "% removal): ", nrow(GRBCut), "\n")

# LASSO feature selection
source("lasso.R")
cat("   LASSO selected vars:", paste(lassovar, collapse=", "), "\n")

# Sqr term generation + response attachment
suppressPackageStartupMessages(library(SuperLearner))
Variables <- subset(GRBCut, select = lassovar)
GRBSqr    <- SqrTermGen(Variables)
GRBSqr$Response <- GRBCut$log10z
GRBSqr    <- GRBSqr[is.finite(GRBSqr$Response), ]
cat("8. After SqrTermGen + finite response:   ", nrow(GRBSqr), "\n")

# What mc_plots sees as training set
cat("   GRBs printed in mc_plots.R / 10fCV:  ", nrow(GRBSqr), "\n")
cat("   (without catOutl removes GRB090429B type GRBs separately in plotting)\n")
