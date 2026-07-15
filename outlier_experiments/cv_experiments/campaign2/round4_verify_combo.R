# Round 4: tighter verification of the two R3 positives (more draws, physical ratios)
# and whether they stack. N_REP=4 for stronger paired sign test.
N_REP <- 4L; K <- 5L; V_INNER <- 5L; N_FRAMES <- 4L
source("../rounds/harness_prep.R")
frames16 <- build_frames(16L)

norm_grb2 <- function(s){s<-trimws(gsub("GRB","",as.character(s))); ifelse(grepl("[A-Za-z]$",s),s,paste0(s,"A"))}
add_ratios <- function(fr){ fr$FaTa<-fr$log10Fa-fr$log10Ta; fr$FluT90<-fr$log10Fluence-fr$log10T90; fr$AlBe<-fr$Alpha-fr$Beta; fr }
frames4_r  <- lapply(frames,   add_ratios)
frames16_r <- lapply(frames16, add_ratios)

base <- default_cfg()
ratio_vars <- c(VARS, "FaTa","FluT90","AlBe")

cat("\n=== isolate draws: base on 4 vs 16 frames ===\n")
r_b4  <- run_config(base, train_all, frames,   "base_nf4")
r_b16 <- run_config(base, train_all, frames16, "base_nf16")
cat(sprintf("draws delta (16-4): d_pool=%+.4f  d_by_rep_mean=%+.4f\n",
            r_b16$r_pool-r_b4$r_pool, mean(r_b16$r_by_rep-r_b4$r_by_rep)))

cat("\n=== isolate ratios at 16 frames ===\n")
r_rat16 <- run_config(modifyList(base, list(vars=ratio_vars)), train_all, frames16_r, "ratios_nf16")
cat(sprintf("ratios delta @16 (rat-base): d_pool=%+.4f  d_by_rep_mean=%+.4f\n",
            r_rat16$r_pool-r_b16$r_pool, mean(r_rat16$r_by_rep-r_b16$r_by_rep)))

cat(sprintf("\nSUMMARY: base_nf4=%.4f  base_nf16=%.4f  ratios_nf16=%.4f  (combined lift vs base_nf4 = %+.4f)\n",
            r_b4$r_pool, r_b16$r_pool, r_rat16$r_pool, r_rat16$r_pool-r_b4$r_pool))
cat("ROUND4_DONE\n")
