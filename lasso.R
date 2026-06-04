#' lasso.R — LASSO feature-importance ranking.
#'
#' Ranks the candidate predictors by how strongly LASSO retains them when
#' regressing log10(z+1) on the imputed feature set. Because a single
#' cross-validated LASSO fit is sensitive to the random fold assignment, we
#' refit `n_repeats` times and average the coefficients at lambda.1se, then
#' rank features by mean |coefficient|.
#'
#' Inputs (globals expected from the calling pipeline):
#'   GRBPred                 imputed feature dataframe (rows = GRBs)
#'   raw_xray_data           raw catalog, supplies Redshift_crosscheck
#'   features_for_mice_preds dataframe whose column names define the candidate
#'                           predictor set fed to LASSO
#'   PLOTaddr                output directory prefix for the diagnostic PNG
#'
#' Output:
#'   LassoFeatures.png  barplot of mean |LASSO coefficient| per feature
#'   lassovar           character vector of predictor names, ordered by
#'                      descending mean |coefficient| (consumed downstream)

library(glmnet)

# Response: log10(z+1), the same target the SuperLearner predicts.
Y <- log10(raw_xray_data$Redshift_crosscheck + 1)

#' Fit a cross-validated LASSO model.
#' @param X matrix or data.frame. Predictor matrix (rows = observations).
#' @param Y numeric. Response vector.
#' @return cv.glmnet object (alpha = 1, i.e. pure LASSO).
LASSO <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.vector(Y)
  lasso_model <- cv.glmnet(X, Y, alpha = 1)
  return(lasso_model)
}

# Refit repeatedly and collect the lambda.1se coefficients so the ranking is
# averaged over many random CV fold splits rather than one lucky/unlucky one.
n_repeats <- 100
lasso_coef <- vector()
for (a in 1:n_repeats) {
  lasmod <- LASSO(Y = Y, X = GRBPred[, colnames(features_for_mice_preds)])
  lasso_coef <- cbind(
    lasso_coef,
    lasmod$glmnet.fit$beta[, lasmod$glmnet.fit$lambda == lasmod$lambda.1se]
  )
}
lasso_coef_avg <- rowMeans(lasso_coef)

# Diagnostic: horizontal barplot of mean |coefficient| per feature.
png(filename = paste(PLOTaddr, "LassoFeatures.png", sep = ""))
par(mar = c(5, 11, 4, 1))
barplot(sort(abs(lasso_coef_avg)),
        horiz = T, las = 1, xlab = "Mean |LASSO coefficient|",
        main = "LASSO feature importance\n(mean |coef| over CV folds; larger = stronger predictor)",
        cex.names = 1.5,
        cex.axis = 1.5, cex.lab = 1.4, cex.main = 0.95,
        font.axis = 2, font.lab = 2)
axis(1, lwd = 3, cex.axis = 1.5)
dev.off()

# Predictor names ranked by descending importance.
lassovar <- names(lasso_coef_avg[order(abs(lasso_coef_avg), decreasing = TRUE)])
