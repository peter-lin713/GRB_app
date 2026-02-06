# Custom SuperLearner Wrappers for Best Caret Models
library(SuperLearner)
library(caret)
library(Cubist)
library(gbm)
library(randomForest)
library(ipred)
library(earth)
library(ranger)


# 1. CUBIST (Best performer)
SL.cubist <- function(Y, X, newX, family, obsWeights = NULL, ...) {
  if (is.null(obsWeights)) obsWeights <- rep(1, length(Y))
  
  fit <- cubist(x = X, y = Y, weights = obsWeights)
  pred <- predict(fit, newdata = newX)
  
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.cubist")
  return(out)
}

predict.SL.cubist <- function(object, newdata, ...) {
  predict(object, newdata = newdata)
}

# 2. GBM (Gradient Boosting)
SL.gbm.custom <- function(Y, X, newX, family, obsWeights = NULL, 
                          n.trees = 1000, interaction.depth = 2, 
                          shrinkage = 0.01, ...) {
  if (is.null(obsWeights)) obsWeights <- rep(1, length(Y))
  
  df <- data.frame(Y = Y, X)
  
  fit <- gbm(Y ~ ., data = df, distribution = "gaussian",
             n.trees = n.trees, interaction.depth = interaction.depth,
             shrinkage = shrinkage, weights = obsWeights, verbose = FALSE)
  
  best_trees <- gbm.perf(fit, method = "OOB", plot.it = FALSE)
  if (is.na(best_trees) || best_trees < 1) best_trees <- n.trees
  
  pred <- predict(fit, newdata = newX, n.trees = best_trees)
  
  out <- list(pred = pred, fit = list(model = fit, best_trees = best_trees))
  class(out$fit) <- c("SL.gbm.custom")
  return(out)
}

predict.SL.gbm.custom <- function(object, newdata, ...) {
  predict(object$model, newdata = newdata, n.trees = object$best_trees)
}


# 3. Random Forest
SL.rf.custom <- function(Y, X, newX, family, mtry = max(floor(ncol(X)/3), 1),
                         ntree = 500, ...) {
  fit <- randomForest(x = X, y = Y, mtry = mtry, ntree = ntree)
  pred <- predict(fit, newdata = newX)
  
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.rf.custom")
  return(out)
}

predict.SL.rf.custom <- function(object, newdata, ...) {
  predict(object, newdata = newdata)
}

# 4. Treebag (Bagged Trees)
SL.treebag <- function(Y, X, newX, family, nbagg = 100, ...) {
  df <- data.frame(Y = Y, X)
  
  fit <- bagging(Y ~ ., data = df, nbagg = nbagg, coob = TRUE)
  pred <- predict(fit, newdata = newX)
  
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.treebag")
  return(out)
}

predict.SL.treebag <- function(object, newdata, ...) {
  predict(object, newdata = newdata)
}

# 5. bagEarth (Bagged MARS)
SL.bagEarth <- function(Y, X, newX, family, degree = 2, nprune = NULL, 
                        nbagg = 50, ...) {
  fit <- bagEarth(x = X, y = Y, degree = degree, nprune = nprune, B = nbagg)
  pred <- predict(fit, newdata = newX)
  
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.bagEarth")
  return(out)
}

predict.SL.bagEarth <- function(object, newdata, ...) {
  predict(object, newdata = newdata)
}

# 6. Ranger (Fast Random Forest)
SL.ranger.custom <- function(Y, X, newX, family, num.trees = 500, 
                             mtry = max(floor(ncol(X)/3), 1), ...) {
  df <- data.frame(Y = Y, X)
  
  fit <- ranger(Y ~ ., data = df, num.trees = num.trees, mtry = mtry)
  pred <- predict(fit, data = newX)$predictions
  
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.ranger.custom")
  return(out)
}

predict.SL.ranger.custom <- function(object, newdata, ...) {
  predict(object, data = newdata)$predictions
}


custom_caret_libs <- c(
  "SL.cubist",
  "SL.gbm.custom", 
  "SL.rf.custom",
  "SL.treebag",
  "SL.bagEarth",
  "SL.ranger.custom"
)

cat("Custom caret wrappers loaded:\n")
print(custom_caret_libs)