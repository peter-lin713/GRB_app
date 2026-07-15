# Round 3: data-side scan. (a) projection-draw count (bagging variance), (b) OT-no-z
# features as extra cols, (c) measurement-error cols as features, (d) physical ratios.
source("../rounds/harness_prep.R")

# ---- (a) frames-count sweep: build 24 draws once, subset ----
frames24 <- build_frames(24L)
cat("\n=== frames-count sweep (base config) ===\n")
for (nf in c(4, 8, 16, 24)) {
  r <- run_config(default_cfg(), train_all, frames24[seq_len(nf)], paste0("nf", nf))
}

# ---- build augmented frames: merge OT-no-z + keep error cols + ratios ----
norm_grb2 <- function(s){s<-trimws(gsub("GRB","",as.character(s))); ifelse(grepl("[A-Za-z]$",s),s,paste0(s,"A"))}
ot <- read.csv("Data/superlearner_training_ot_v3_noz_errcut_relative.csv", stringsAsFactors=FALSE)
names(ot)[1]<-"GRB"; ot$GRB<-norm_grb2(ot$GRB)
ot_feats <- c("log10Fa","log10Ta","Alpha","Beta")
err_feats <- c("log10FaErr","log10TaErr","AlphaErr","BetaErr","log10FluenceErr","PhotonIndexErr","log10PeakFluxErr")
mo2 <- match(f$GRB, ot$GRB)
augment <- function(fr) {
  for (v in ot_feats) fr[[paste0(v,"_ot")]] <- ot[[v]][mo2]
  # ratios / physical combos
  fr$FaTa   <- fr$log10Fa - fr$log10Ta      # plateau lum-duration
  fr$FluT90 <- fr$log10Fluence - fr$log10T90
  fr$AlBe   <- fr$Alpha - fr$Beta
  fr
}
frames_aug <- lapply(frames, augment)

base <- default_cfg()
ot_vars    <- c(VARS, paste0(ot_feats,"_ot"))
err_vars   <- c(VARS, err_feats)
ratio_vars <- c(VARS, "FaTa","FluT90","AlBe")
cfgs <- list(
  base      = base,
  ot_extra  = modifyList(base, list(vars = ot_vars)),
  err_feats = modifyList(base, list(vars = err_vars)),
  ratios    = modifyList(base, list(vars = ratio_vars)),
  ot_err_rat= modifyList(base, list(vars = c(VARS, paste0(ot_feats,"_ot"), err_feats, "FaTa","FluT90","AlBe")))
)
cat("\n=== data-augmentation configs (on 4 base frames) ===\n")
t0 <- Sys.time()
res <- run_set(cfgs, train_all, frames_aug)
cat("elapsed:", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
cat("ROUND3_DONE\n")
