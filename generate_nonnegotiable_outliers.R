#!/usr/bin/env Rscript
#' generate_nonnegotiable_outliers.R -- single source of truth for the
#' non-negotiable-cut outlier list. Defines update_nonnegotiable_outliers(),
#' callable two ways:
#'   1. Sourced from a formula-generation script (preferred -- keeps the list
#'      in sync with whatever data that script is actually using):
#'        source("generate_nonnegotiable_outliers.R")
#'        update_nonnegotiable_outliers(raw_xray_data, out_dir = ".", method = "emcee_v2")
#'      `raw_xray_data` must already be scope-filtered (log10T90 > log10(2)).
#'   2. As a standalone script (reads and scope-filters the CSV itself):
#'        Rscript generate_nonnegotiable_outliers.R <recovered_csv> <out_dir> <method_label>
#'
#' err>value is restricted to Alpha, Beta, log10Fa, log10Ta only -- these are
#' the four parameters emcee's MCMC calibration actually fits FOR THE
#' OPTICALLY-SOURCED SUBSET (is_optical == 1). T90, PhotonIndex, Fluence, and
#' NH are directly-measured catalog quantities for every GRB, so a large
#' relative error there reflects genuine measurement uncertainty, not an
#' unreliable projection fit -- applying err>value to them would conflate two
#' different kinds of "uncertain."
#'
#' Note this check runs on ALL GRBs, not just the is_optical==1 subset it was
#' motivated by: for X-ray-native GRBs (is_optical == 0), Alpha/Beta/Fa/Ta are
#' also direct measurements, not emcee-fit values, so a flag there means
#' "poorly-constrained measurement of this parameter," not "unreliable emcee
#' projection" -- a related but distinct kind of unreliability, still on the
#' same four parameters. In practice this matters: as of the current dataset,
#' all 4 non-negotiable-cut flags are on is_optical==0 GRBs, none on the
#' optically-projected subset the check was originally designed for.
#' PeakFlux was excluded earlier for an unrelated reason (unstable/unexplained,
#' see the emcee/OT PeakFlux err>value investigation).
#'
#' candidate_outliers_<method>.csv is fully auto-computed and OVERWRITTEN every
#' call -- never hand-edit it. confirmed_outliers_to_drop.txt is the actual
#' drop list: this function UNIONS the freshly-detected automatic candidates
#' into whatever's already there, so manual additions from visual review of
#' the 4D plots persist across reruns, but the automatic part can never go
#' stale (the bug this whole refactor exists to prevent: a stale hand-curated
#' file that nothing kept in sync with the current data).

update_nonnegotiable_outliers <- function(d, out_dir = ".", method) {
  hard <- d$Alpha > 3 | d$Beta > 2 | d$log10T90 > 6 | d$PhotonIndex < 0

  err_ratio <- function(v, e) e / abs(v)
  alpha_err <- err_ratio(d$Alpha,   d$AlphaErr)   > 0.5
  beta_err  <- err_ratio(d$Beta,    d$BetaErr)    > 0.5
  fa_err    <- err_ratio(d$log10Fa, d$log10FaErr) > 0.5
  ta_err    <- err_ratio(d$log10Ta, d$log10TaErr) > 0.5

  flagged <- hard | alpha_err | beta_err | fa_err | ta_err
  flagged[is.na(flagged)] <- FALSE

  reason_for <- function(i) {
    r <- character(0)
    if (isTRUE(d$Alpha[i] > 3))       r <- c(r, "Alpha>3")
    if (isTRUE(d$Beta[i] > 2))        r <- c(r, "Beta>2")
    if (isTRUE(d$log10T90[i] > 6))    r <- c(r, "log10T90>6")
    if (isTRUE(d$PhotonIndex[i] < 0)) r <- c(r, "PhotonIndex<0")
    if (isTRUE(alpha_err[i])) r <- c(r, "err>val:Alpha")
    if (isTRUE(beta_err[i]))  r <- c(r, "err>val:Beta")
    if (isTRUE(fa_err[i]))    r <- c(r, "err>val:log10Fa")
    if (isTRUE(ta_err[i]))    r <- c(r, "err>val:log10Ta")
    paste(r, collapse = "; ")
  }

  idx <- which(flagged)
  out <- data.frame(
    GRB         = rownames(d)[idx],
    is_optical  = d$is_optical[idx],
    reasons     = vapply(idx, reason_for, character(1)),
    Alpha       = d$Alpha[idx],
    Beta        = d$Beta[idx],
    log10T90    = d$log10T90[idx],
    PhotonIndex = d$PhotonIndex[idx],
    log10Fa     = d$log10Fa[idx],
    log10Ta     = d$log10Ta[idx]
  )
  out <- out[order(out$GRB), ]

  candidate_file <- file.path(out_dir, sprintf("candidate_outliers_%s.csv", method))
  confirm_file   <- file.path(out_dir, "confirmed_outliers_to_drop.txt")

  write.csv(out, candidate_file, row.names = FALSE)

  manual_lines    <- if (file.exists(confirm_file)) readLines(confirm_file) else character(0)
  manual_comments <- manual_lines[grepl("^#", manual_lines) | !nzchar(trimws(manual_lines))]
  manual_ids      <- trimws(manual_lines[!grepl("^#", manual_lines) & nzchar(trimws(manual_lines))])
  to_drop <- union(manual_ids, out$GRB)

  writeLines(c(
    if (length(manual_comments)) manual_comments else c(
      sprintf("# Final drop list (%s). err>value restricted to Alpha/Beta/log10Fa/log10Ta only", method),
      "# -- the 4 parameters emcee's MCMC fits for optically-sourced GRBs; for",
      "# X-ray-native GRBs these are direct measurements instead, so a flag there",
      "# means a poorly-constrained measurement, not an unreliable emcee fit.",
      "# Checked on all GRBs regardless. PeakFlux and other directly-measured",
      "# columns (T90, PhotonIndex, Fluence, NH) excluded from err>value on",
      "# principle: a large relative error there is catalog measurement",
      "# uncertainty, not an unreliable fit, for every GRB. Hard value cuts",
      "# (Alpha>3, Beta>2, log10T90>6, PhotonIndex<0) still apply regardless.",
      "# Auto-candidates are unioned in fresh every run (generate_nonnegotiable_",
      "# outliers.R); add GRBs below by hand after visual review of the 4D",
      "# plots to extend this list. One GRB per line.",
      "#"),
    to_drop
  ), confirm_file)

  cat(sprintf("[%s] N=%d | flagged=%d (hard=%d, err>value Alpha/Beta/Fa/Ta=%d) | drop list: %d manual + %d automatic = %d unique\n",
              method, nrow(d), length(idx), sum(hard), sum(alpha_err | beta_err | fa_err | ta_err),
              length(manual_ids), length(out$GRB), length(to_drop)))
  cat("  ->", paste(to_drop, collapse = ", "), "\n")
  cat("Wrote:", candidate_file, "and", confirm_file, "\n")

  invisible(list(candidates = out, to_drop = to_drop))
}

# ---- CLI entry point --------------------------------------------------------
# Only runs when this file itself is the one Rscript was invoked on (i.e. NOT
# when another script source()s it -- in that case --file= points at the
# OUTER script, not this one).
.this_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.is_main <- length(.this_file_arg) > 0 &&
  basename(sub("^--file=", "", .this_file_arg)) == "generate_nonnegotiable_outliers.R"

if (.is_main) {
  args <- commandArgs(trailingOnly = TRUE)
  recovered_csv <- args[1]
  out_dir       <- args[2]
  method        <- args[3]

  d <- read.csv(recovered_csv, row.names = 1)
  d <- d[d$log10T90 > log10(2), ]  # scope filter: long GRBs only, matches GAM_analysis_8variables.R
  update_nonnegotiable_outliers(d, out_dir, method)
}
