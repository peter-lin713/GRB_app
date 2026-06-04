#' Load_Imports.R — package bootstrap for the GRB redshift pipeline.
#'
#' Sourced once at the top of the pipeline. Points the session at a CRAN
#' mirror, installs any of the required packages that are missing, then
#' attaches them all. Idempotent: re-sourcing only installs what is absent.
#'
#' Package roles in the pipeline:
#'   doParallel    parallel backend (legacy; current code uses parallel::mclapply)
#'   mice          multiple imputation of missing predictors
#'   VIM           missing-data visualisation (mice diagnostics)
#'   ggplot2       result / MC diagnostic plots
#'   lattice       predicted-vs-observed plots
#'   stringr       string manipulation
#'   dplyr         data wrangling / sorting
#'   MASS          rlm() M-estimator outlier cut; normal-fit utilities
#'   randomForest  SuperLearner base learner
#'   earth         SuperLearner base learner (MARS)
#'   glmnet        LASSO feature selection / ElasticNet
#'   SuperLearner  ensemble meta-learner
#'   mgcv          GAM base learner (custom SL.mgcv_gam wrapper)
#'   xgboost       gradient-boosted-tree base learner
#'   gbm           gradient-boosting base learner
#'   caret         training utilities
#'   party         conditional-inference trees/forests
#'   arm           bayesglm() base learner
#'   dgof, kSamples  goodness-of-fit / k-sample tests
#'   latticeExtra  lattice plot extensions
#'   biglasso      out-of-core LASSO
#'   Matrix        sparse-matrix support (glmnet dependency)
#'   e1071         SVM / misc ML utilities

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Packages required across the whole pipeline.
packages <- c("doParallel", "mice", "VIM", "ggplot2", "lattice", "stringr", "dplyr", "MASS", "randomForest", "earth", "glmnet", "SuperLearner", "mgcv", "xgboost", "gbm", "caret", "party", "arm", "dgof", "kSamples", "latticeExtra", "biglasso", "Matrix", "e1071")

# Install any packages not yet present, then attach them all.
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))
