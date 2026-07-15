# Round 5: Daume c-sweep (priority #2, incl c=0 = drop optical copies) + untried
# nonlinear learners (earth/MARS, xgboost) as a final model-side probe.
N_REP <- 4L; K <- 5L; V_INNER <- 5L; N_FRAMES <- 4L
source("../rounds/harness_prep.R")

SL.earth2 <<- function(Y,X,newX,family,obsWeights,...){
  df <- data.frame(Y=Y, X[, intersect(linear_vars, colnames(X))], check.names=FALSE)
  fit <- earth::earth(Y~., data=df, degree=2, penalty=3, nk=21)
  nd <- data.frame(newX[, intersect(linear_vars, colnames(newX))], check.names=FALSE)
  out <- list(object=fit); class(out) <- "SL.earth2"
  list(pred=as.numeric(predict(fit, newdata=nd)), fit=out)
}
predict.SL.earth2 <<- function(object,newdata,...)
  as.numeric(predict(object$object, newdata=data.frame(newdata, check.names=FALSE)))

SL.xgb2 <<- function(Y,X,newX,family,obsWeights,...){
  vv <- intersect(c(linear_vars, paste0(top7,"Sqr"),"is_optical"), colnames(X))
  dtr <- xgboost::xgb.DMatrix(as.matrix(X[,vv]), label=Y)
  fit <- xgboost::xgboost(data=dtr, nrounds=150, max_depth=3, eta=0.05,
                          subsample=0.8, colsample_bytree=0.8, verbose=0, nthread=1)
  out <- list(object=fit, vv=vv); class(out) <- "SL.xgb2"
  list(pred=as.numeric(predict(fit, as.matrix(newX[,vv]))), fit=out)
}
predict.SL.xgb2 <<- function(object,newdata,...)
  as.numeric(predict(object$object, as.matrix(newdata[,object$vv])))

base <- default_cfg()
cfgs <- list(
  base       = base,
  daume_c0   = modifyList(base, list(c_scale = 0.0)),
  daume_c0.5 = modifyList(base, list(c_scale = 0.5)),
  daume_c1.0 = modifyList(base, list(c_scale = 1.0)),
  add_earth  = modifyList(base, list(extra_learners = "SL.earth2")),
  add_xgb    = modifyList(base, list(extra_learners = "SL.xgb2")),
  add_both   = modifyList(base, list(extra_learners = c("SL.earth2","SL.xgb2")))
)
t0 <- Sys.time()
res <- run_set(cfgs)
cat("elapsed:", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
cat("ROUND5_DONE\n")
