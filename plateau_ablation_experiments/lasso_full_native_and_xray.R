#!/usr/bin/env Rscript
# LASSO ranking on the fully-populated dataset: both the standard X-ray-scale
# plateau params AND their optical-native-scale counterparts (real where
# measured, emcee-calibration-projected -- in whichever direction was needed
# -- everywhere else), for every GRB. Answers: does LASSO prefer the X-ray
# or the native-optical representation of Fa/Ta/Alpha/Beta?
suppressMessages(library(glmnet))

setwd("/Users/petelin/Desktop/test/GRB_app/plateau_ablation_experiments")
d <- read.csv("full_native_and_xray.csv", header = TRUE)

FEATURES <- c("log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
              "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux",
              "log10Fa_native", "log10Ta_native", "Alpha_native", "Beta_native")
X <- as.matrix(d[, FEATURES])
Y <- log10(d$Redshift_crosscheck + 1)

LASSO <- function(X, Y) cv.glmnet(X, Y, alpha = 1)

lasso_coef <- vector()
for (a in 1:100) {
  lasmod <- LASSO(X, Y)
  lasso_coef <- cbind(lasso_coef, lasmod$glmnet.fit$beta[, lasmod$glmnet.fit$lambda == lasmod$lambda.1se])
}
lasso_coef_avg <- rowMeans(lasso_coef)

cat("\n--- Mean LASSO coefficients (full native + X-ray dataset, N=", nrow(d), ") ---\n", sep = "")
print(lasso_coef_avg)

ranked <- names(sort(abs(lasso_coef_avg), decreasing = TRUE))
cat("\nFull ranking by |coefficient|:\n")
print(ranked)

cat("\nX-ray-scale vs native-optical-scale, side by side:\n")
pairs <- list(c("log10Fa","log10Fa_native"), c("log10Ta","log10Ta_native"),
              c("Alpha","Alpha_native"), c("Beta","Beta_native"))
for (p in pairs) {
  cat(sprintf("  %-16s coef=%9.5f  rank=%2d/14   |   %-16s coef=%9.5f  rank=%2d/14\n",
              p[1], lasso_coef_avg[p[1]], which(ranked == p[1]),
              p[2], lasso_coef_avg[p[2]], which(ranked == p[2])))
}

saveRDS(list(coef = lasso_coef_avg, ranked = ranked), "lasso_full_native_and_xray_result.rds")

# ---- Plot -------------------------------------------------------------------
png(filename = "LassoFeatures_full_native_and_xray.png", width = 950, height = 780)
par(mar = c(5, 13, 4, 1))
ord <- order(abs(lasso_coef_avg))
is_native <- grepl("_native$", names(lasso_coef_avg)[ord])
is_xray_plateau <- names(lasso_coef_avg)[ord] %in% c("log10Fa","log10Ta","Alpha","Beta")
cols <- ifelse(is_native, "#eb6834", ifelse(is_xray_plateau, "#eda100", "#2a78d6"))
barplot(abs(lasso_coef_avg)[ord],
        horiz = TRUE, las = 1, xlab = "Mean |LASSO coefficient|",
        main = "LASSO feature importance: X-ray-scale vs native-optical-scale plateau params\n(orange = native-optical, yellow = X-ray-scale plateau, blue = other)",
        col = cols,
        cex.names = 1.15, cex.axis = 1.1, cex.lab = 1.1, cex.main = 1.0,
        font.axis = 2, font.lab = 2)
axis(1, lwd = 2, cex.axis = 1.1)
dev.off()
cat("\nSaved LassoFeatures_full_native_and_xray.png\n")
