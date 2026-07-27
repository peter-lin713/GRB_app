#source('Load_Imports.R')

##### LOADING UP PACKAGES #####

require(doParallel)
require(mice)
#require(VIM) # mice graph 1
require(ggplot2) # mice graph 2
require(lattice) # GAM predicted vs observed plot
require(stringr) # string remove
require(dplyr) # arrange (sorting for GAM)
require(MASS) # fit normal dist in GAM (super balanced sampling)
require(randomForest)
require(earth)
library(glmnet) # LASSO, ElasticNet
require(SuperLearner)
require(mgcv)
require(xgboost)
require(gbm)
require(caret)
#require(cforest)
require(biglasso)
require(MASS)
require(dplyr)
#library(VIM)
#library(biglasso)
library(arm)
#require(ck37r) # SL.mgcv

##### DEFINING THE CUSTOM GAM AND GLM FUNCTIONS #####
# THIS IS DONE SO THAT IT DOES NOT NEED TO READ IN
# THE FUNCTION FILES. HELPS REDUCE OCCURANCE OF ERRORS

####### GAM #####
SL.mgcv_gam <- function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, gam.model = NA, verbose = F, ...) {
  # using require instead of requireNamespace() to allow the formula to parse correctly with s(), gam::s() doesn't work, is not recognized as a special function
  #if(!require('gam')) {stop("SL.gam requires the gam package, but it isn't available")} 
  #if("mgcv" %in% loadedNamespaces()) warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  # create the formula for gam with a spline for each continuous variable
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if(is.atomic(gam.model)){
    if (sum(!cts.x) > 0) { 
      gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
    } else {
      gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
    }
    # fix for when all variables are binomial
    if (sum(!cts.x) == length(cts.x)) {
      gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
    }
  } else {
    
    if(verbose){print(gam.model)}
    fn <- strsplit(as.character(gam.model),split = '~')
    gam.model <- as.formula(paste('Y ~',fn[[3]][1])) # CREATING THE FORMULA WITH THE RESPONSE VARIABLE SET AS Y
  }
  
  if(verbose){print(gam.model)}
  
  
  fit.gam <- mgcv::gam(gam.model, data = X, family = family, control = mgcv::gam.control(maxit = 50), weights = obsWeights)
  
  #fit.gam <- mgcv::gam(test_formula, data = TrainingData[,-c(1,3)], family = gaussian(), control = mgcv::gam.control(maxit = 50))#, weights = obsWeights)
  
  if(verbose){print('fitted')}
  
  if(packageVersion('gam') >= 1.15) {
    pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv_gam")
  return(out)
}

predict.SL.mgcv_gam <- function(object, newdata, ...){
  #.SL.require('mgcv')
  if(packageVersion('gam') >= 1.15) {
    pred <- mgcv::predict.gam(object = object$object, newdata = newdata, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  
  return(pred)
}

####### GLM ######

require('arm')

SL.custom_glm <- function(Y, X, newX, family, obsWeights, glm.model = NA, ...){
  
  
  if(is.atomic(glm.model)){
    glm.model <- as.formula(paste("Y~", paste(colnames(X), collapse="+")))
  } else {
    print(glm.model)
    fn <- strsplit(as.character(glm.model),split = '~')
    glm.model <- as.formula(paste('Y ~',fn[[3]][1])) # CREATING THE FORMULA WITH THE RESPONSE VARIABLE SET 
    print(glm.model)
  }
  
  fit.glm <- glm(glm.model, data = X, family = family, weights = obsWeights)
  
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.custom_glm")
  return(out)
}

predict.SL.custom_glm <- function(object, newdata, ...){
  #.SL.require('arm')
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  return(pred)
}

####### READING IN THE DATA #####

run_locally = T

if(run_locally){
  raw_xray_data <- read.csv("superlearner_training_ot_v3_native_opt.csv", header = TRUE, row.names = 1)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- args[1]
  output_file <- args[2]
  raw_xray_data <- read.csv(input_file, header = TRUE, row.names = 1)
}

do_mice = T

###USERS INPUTS THEIR OWN TRAINING SET INSTEAD OF THIS ONE
#GRBPred <- read.csv(file = "OutputFiles/MEstimator/grb_xray_m_est.csv", header = TRUE, row.names = 1) 

# colnames(GRBPred)

##### NEEDS TO BE EDITED WHEN WE WORK ON MICE
# checking for NA values, if using MICE, should have none
#md.pattern(GRBPred, rotate.names = T)


SqrTermGen <- function(inputData) {
  indVar <- colnames(inputData)
  for (i in 1:length(indVar)) { # Loop over all variables in indVar
    for (j in i:(length(indVar))) { # Loop over all variables at index i and greater than index i
      # If i and j correspond to the same varaible call the variable varSqr
      if (indVar[i] == indVar[j]) {
        inputData[[paste(indVar[i], "Sqr", sep = "")]] <- inputData[, indVar[i]] * inputData[, indVar[j]]
      } # else{
      # inputData[[paste(indVar[i],indVar[j],sep="")]] <- inputData[,indVar[i]]*inputData[,indVar[j]]
      # }
    }
  }
  return(inputData)
}

#addr <-paste("Results_",Sys.Date(),"/Files/",sep='') #THIS HOLDS THE ADDRESS AT WHICH THE FILES ARE OUTPUT
PLOTaddr <-paste("Plot_Output") #THIS HOLDS THE ADDRESS AT WHICH THE PLOTS ARE OUTPUT

sz<-0.8
rez=120


## CREATE DIRECTORIES IF THEY DONT EXIST
# if(!dir.exists(PLOTaddr)){dir.create(PLOTaddr)}
# if(!dir.exists(addr)){dir.create(addr)}


# Pick out Long GRBs — handle both old format (T90 linear) and new format (log10T90 directly)
if ("T90" %in% colnames(raw_xray_data)) {
  raw_xray_data = raw_xray_data[raw_xray_data$T90 > 2,]
  raw_xray_data$log10T90 = log10(raw_xray_data$T90)
} else {
  raw_xray_data = raw_xray_data[raw_xray_data$log10T90 > log10(2),]
}


# creating new subset for only numeric features that are not tied to the response (Redshift_crosscheck)
# The 4 "_native" columns are the EXTRA candidates for this experiment: real
# native optical-band Alpha/Beta/Fa/Ta, present only for the 73 optically-
# sourced GRBs (NA -- and thus MICE-imputed -- for the other 223). They sit
# alongside, not in place of, the standard OT-projected Alpha/Beta/Fa/Ta,
# so LASSO's own ranking decides which version (if either) carries signal.
features_for_mice_preds = subset(raw_xray_data,select = c(log10T90,
                                                          log10Fa,
                                                          log10Ta,
                                                          Alpha,
                                                          Beta,
                                                          Gamma,
                                                          log10Fluence,
                                                          PhotonIndex,
                                                          log10NH,
                                                          log10PeakFlux,
                                                          log10Fa_native,
                                                          log10Ta_native,
                                                          Alpha_native,
                                                          Beta_native))

# ---- Non-negotiable outlier detection: single source of truth --------------
# Delegates to generate_nonnegotiable_outliers.R (repo root) so there is only
# ONE place that implements this logic, not one copy per formula-generation
# folder. Recomputes candidate_outliers_*.csv and unions into
# confirmed_outliers_to_drop.txt fresh every run -- see that file's header for
# the full rationale (err>value restricted to Alpha/Beta/log10Fa/log10Ta only,
# manual additions from 4D-plot review preserved across reruns).
source("../../generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = ".", method = "ot_v3_native_opt")
to_drop    <- .outliers$to_drop
reason_for <- setNames(.outliers$candidates$reasons, .outliers$candidates$GRB)

grb_ids      <- rownames(raw_xray_data)
dropped_mask <- grb_ids %in% to_drop
if (any(dropped_mask)) {
  dropped_log <- data.frame(
    GRB    = grb_ids[dropped_mask],
    reason = ifelse(grb_ids[dropped_mask] %in% names(reason_for),
                     reason_for[grb_ids[dropped_mask]], "manual/visual review")
  )
  write.csv(dropped_log, "removed_outliers_final.csv", row.names = FALSE)
  cat("Removed", sum(dropped_mask), "GRB(s) confirmed physically infeasible -- see removed_outliers_final.csv\n")
  raw_xray_data           <- raw_xray_data[!dropped_mask, ]
  features_for_mice_preds <- features_for_mice_preds[!dropped_mask, ]
} else {
  cat("confirmed_outliers_to_drop.txt is empty or matched nothing -- no GRBs removed.\n")
}

# Safety net only -- everything above should already be gone. Kept as NA+impute
# (not a second removal pass) purely to catch anything that slips through
# unreviewed; should find nothing to touch if the checkpoint above did its job.
features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20]       <- NA
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 2]              <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3]            <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3]            <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

# Fluence/PeakFlux errors may be stored linear (FluenceErr/PeakFluxErr) or
# already dex (log10FluenceErr/log10PeakFluxErr) depending on the source
# catalog -- detect which is present instead of hardcoding the linear names.
fluence_err_col  <- if ("FluenceErr"  %in% colnames(raw_xray_data)) "FluenceErr"  else "log10FluenceErr"
peakflux_err_col <- if ("PeakFluxErr" %in% colnames(raw_xray_data)) "PeakFluxErr" else "log10PeakFluxErr"

features_for_mice_errs = raw_xray_data[, c("T90Err", "log10FaErr", "log10TaErr",
                                            "AlphaErr", "BetaErr", fluence_err_col,
                                            "PhotonIndexErr", peakflux_err_col,
                                            "log10FaErr_native", "log10TaErr_native",
                                            "AlphaErr_native", "BetaErr_native")]

if(do_mice){
  set.seed(1)

  # features_for_mice_all <- cbind(features_for_mice_preds, features_for_mice_errs)
  # The 4 native columns are strongly linearly related to their non-native
  # counterparts by construction (that's literally how the emcee projection
  # was built) -- this trips mice's global find.collinear() screen, which
  # drops them to method="" (never imputed, left all-NA), which then crashes
  # the downstream glmnet LASSO call with "x has missing values". Fix:
  # exclude them from being a predictor for anything else (removes them from
  # that screen entirely, since it only considers columns used as a predictor
  # for >=1 other variable) and keep their own predictor row away from the
  # plateau cluster (avoids a near-singular fit inside midastouch's own
  # per-iteration regression). Same issue, same fix, for the *_native error
  # columns below.
  native_cols   <- c("log10Fa_native", "log10Ta_native", "Alpha_native", "Beta_native")
  plateau_cols  <- c("log10Fa", "log10Ta", "Alpha", "Beta")
  pred_mat_preds <- mice::make.predictorMatrix(features_for_mice_preds)
  pred_mat_preds[, native_cols] <- 0
  for (nc in native_cols) pred_mat_preds[nc, c(plateau_cols, native_cols)] <- 0

  mice_model_preds <- mice(data = features_for_mice_preds,
                           m = 20,
                           method = 'midastouch',
                           predictorMatrix = pred_mat_preds,
                           printFlag = F)
  features_for_mice_preds <- complete(mice_model_preds,20)

  err_native_cols  <- c("log10FaErr_native", "log10TaErr_native", "AlphaErr_native", "BetaErr_native")
  err_plateau_cols <- c("log10FaErr", "log10TaErr", "AlphaErr", "BetaErr")
  pred_mat_errs <- mice::make.predictorMatrix(features_for_mice_errs)
  pred_mat_errs[, err_native_cols] <- 0
  for (nc in err_native_cols) pred_mat_errs[nc, c(err_plateau_cols, err_native_cols)] <- 0

  mice_model_errs <- mice(data = features_for_mice_errs,
                          m = 20,
                          method = 'midastouch',
                          predictorMatrix = pred_mat_errs,
                          printFlag = F)
  features_for_mice_errs <- complete(mice_model_errs,20)

  GRBPred <- cbind(features_for_mice_preds, features_for_mice_errs)
  # GRBPred<-complete(mice_model,20)
  
} else {
  #NEEDS TO BE EXPANDED IF NO MICE
  GRBPred <- rs_data_preds
}

########### Adding log error columns ###########
# Convert linear-scale errors to dex only if the dex columns aren't already present
if (!"log10T90Err" %in% colnames(GRBPred)) {
  T90 <- 10^GRBPred$log10T90
  GRBPred$log10T90Err <- GRBPred$T90Err / (T90 * log(10))
}
if (!"log10FluenceErr" %in% colnames(GRBPred)) {
  Fluence <- 10^GRBPred$log10Fluence
  GRBPred$log10FluenceErr <- GRBPred$FluenceErr / (Fluence * log(10))
}
if (!"log10PeakFluxErr" %in% colnames(GRBPred)) {
  PeakFlux <- 10^GRBPred$log10PeakFlux
  GRBPred$log10PeakFluxErr <- GRBPred$PeakFluxErr / (PeakFlux * log(10))
}


######### LASSO FEATURE SELECTION ######

library(glmnet)

LASSO <- function(X,Y)
{
  X<-as.matrix(X) # THE TRAINING DATA
  Y<-as.vector(Y) # THE RESPONSE VECTOR
  lasso_model<-cv.glmnet(X,Y,alpha=1) # LASSO REGRESSION
  return(lasso_model)
}

lasso_coef=vector()

for (a in 1:100) {
  lasmod=LASSO(Y = log10(raw_xray_data$Redshift_crosscheck + 1), X = features_for_mice_preds)
  lasso_coef<-cbind(lasso_coef,lasmod$glmnet.fit$beta[,lasmod$glmnet.fit$lambda==lasmod$lambda.1se])
}
lasso_coef_avg = rowMeans(lasso_coef)

# png(filename = paste(PLOTaddr,'LassoFeatures.png',sep = ''))#,width = 600,height = 1000)
# 
# par(mar=c(5,9,1,1))
# barplot(sort(abs(lasso_coef_avg))
#         ,horiz = T,las=1,xlab="Coefficient"
#         ,cex.names = 1.5
#         ,cex.axis = 1.5,cex.lab=1.75
#         ,font.axis=2,font.lab=2
# )
# axis(1,lwd=3,cex.axis=1.5)
# dev.off()
print(lasso_coef_avg)
lassovar = names(lasso_coef_avg[order(abs(lasso_coef_avg), decreasing=TRUE)])

##### SELECT HOW MANY VARIABLES ####
lassovar=head(lassovar,9)  # experiment: paper keeps top-7, but here we test whether
# ranks 8-9 (weak/near-zero LASSO coefficients on the standard data) do anything once
# combined with the native-optical Alpha/Beta/Fa/Ta substitution -- see the 9-variable
# cost/signal discussion this folder was built to test. 9 vars -> 18 O2-columns ->
# 2^18-1 = 262,143 formulas (16x the paper's 16,383), hence supercomputer, not laptop.


# generating squared terms for future ML methods
Variables <- subset(GRBPred, select = c(lassovar))
Err <- subset(GRBPred, select = colnames(GRBPred) %in% paste(lassovar,"Err",sep=""))

#GRBPred <- SqrTermGen(Variables)
#GRBPred <- cbind(GRBPred, Err)

# adding z back into original dset
GRBPred$Redshift_crosscheck <- raw_xray_data$Redshift_crosscheck
# unscaled_GRBPred$Redshift_crosscheck <- GRBPred$Redshift_crosscheck

# adding log10z back into original dataset
GRBPred$log10z <- log10(raw_xray_data$Redshift_crosscheck + 1)
# unscaled_GRBPred$log10z <- GRBPred$log10z 

# GRBPred$invz <- 1/(raw_xray_data$Redshift_crosscheck + z_e)
# unscaled_GRBPred$invz <- GRBPred$invz


# writing file to output directory
if (do_mice){
  write.csv(GRBPred, "grb_xray_imputed.csv")
  # write.csv(GRB_Err, "OutputFiles/DataHandle/grb_xray_errors.csv")
} else {
  write.csv(GRBPred, "grb_xray.csv")
  # write.csv(GRB_Err, "OutputFiles/DataHandle/grb_xray_errors.csv")
}


#WE'LL DECIDE IF WE NEED IT
##source("m_estimator.R")

#PredictionData <- tail(GRBPred, n = 0.20 * nrow(GRBPred))

#TrainData <- GRBPred[!(rownames(GRBPred) %in% rownames(PredictionData)),]

#Response <- TrainData$log10z
#Predictors <- subset(TrainData
#                     ,select = -c(log10z, Redshift_crosscheck)) # EXCLUDING LOG10Z, INVZ AND Z

TrainData=GRBPred
#TrainData=read.csv("Total_useable_data_w_MICE_NH&Peak_imputed.csv",header = T,row.names = 1)
PredictionData <- tail(TrainData,n = nrow(TrainData)*(0.1))



### FORMULA GENERATION #####


#source('Formula_Generation/GAM_formula_generator.R')
# NEW GAM FORMULA GENERATION
# THIS USES THE combn FUNCTION
# WITH THE PARAMETER m LOOPED OVER FROM
# 1 TO 9
# THIS SHOULD GENERATE MANY COMBINATIONS

#Function to add second order data to the input data set and return this modified data set
SqrTermGen <- function(inputData){
  
  indVar = colnames(inputData)
  
  for (i in 1:length(indVar)){ #Loop over all variables in indVar
    for (j in i:(length(indVar))){ #Loop over all variables at index i and greater than index i
      #If i and j correspond to the same varaible call the variable varSqr
      if (indVar[i] == indVar[j]){
        inputData[[paste(indVar[i],"Sqr",sep="")]] <- inputData[,indVar[i]]*inputData[,indVar[j]]
      }#else{
      # inputData[[paste(indVar[i],indVar[j],sep="")]] <- inputData[,indVar[i]]*inputData[,indVar[j]]
      #}
    }
  }
  
  return(inputData) 
}

quadTermGen <- function(inputData){
  #indVar = c("LogFluence","LogT90","LogPeak","PhotonIndex","Alpha","LogTa","LogNH","LogFlux","Gamma")
  #indVar = c("LogFlux1.100m", "LogEnergy_Flux100", "Frac_Variability", "LogHighest_Energy", "PL_Index", "LogPivot_Energy", "LP_Index", "LP_beta","Gaia_G_Magnitude","w4","Lognufnu","Lognu")
  #indVar = NumOvar[c(-1,-2)]
  #indVar = NumOvar[c(-1,-2,-3,-4)]
  indVar = colnames(inputData)
  #print(indVar)
  #indVar
  for (i in 1:length(indVar)){ #Loop over all variables in indVar
    for (j in i:(length(indVar))){ #Loop over all variables at index i and greater than index i
      #If i and j correspond to the same varaible call the variable varSqr
      if (indVar[i] == indVar[j]){
        inputData[[paste(indVar[i],"Sqr",sep="")]] <- inputData[,indVar[i]]*inputData[,indVar[j]]
      }else{
        inputData[[paste(indVar[i],indVar[j],sep="")]] <- inputData[,indVar[i]]*inputData[,indVar[j]]
      }
    }
  }
  
  return(inputData)
}

#TrainData = read.csv("ServerTrainingData.csv",header = T,row.names = 1)
#TrainData = read.csv(paste("OutputFiles/MEstimator","/grb_xray_m_est.csv",sep = ""),header = T,row.names = 1)

Response = TrainData$log10z

Predictor = subset(TrainData,select = lassovar)

O1Predictor_names = colnames(Predictor)

SQRPredictors = SqrTermGen(Predictor)

O2Predictors = quadTermGen(Predictor)

{
  O2=T
  if(O2){
    Predictor_names = colnames(SQRPredictors)
  }else{
    Predictor_names = colnames(Predictor)
  }
  
  term_matrix=c();
  
  for (i in 1:length(Predictor_names)) {
    term_matrix[[i]] <- combn(Predictor_names, m = i)
  }
  
  dim(term_matrix[[3]])[2]
  
  formula_list=vector()
  index=1
  gam_summary=c()
  
  for (m in 1:length(Predictor_names)) {
    
    j = dim(term_matrix[[m]])[2]
    print(j)
    
    for(j in 1:dim(term_matrix[[m]])[2] ){
      
      test1=c()
      for(i in 1:dim(term_matrix[[m]])[1]){
        
        if(dim(term_matrix[[m]])[1] != 1){
          test1 = paste( term_matrix[[m]][i,j], test1 ,sep="+")
        }else{
          test1 = paste( term_matrix[[m]][i,j],":",term_matrix[[m]][i,j],"+",sep="")
        }
      }  
      #print(Predictor_names[!(Predictor_names %in% term_matrix[[3]][,1])])
      
      term_to_add = Predictor_names[!(Predictor_names %in% term_matrix[[m]][,j])]
      
      if(dim(term_matrix[[m]])[1] != 1){
        test1 = paste("(",substr(test1,1,str_length(test1)-1),")^2",sep="") # str_length(xx) - 1 TAKES CARE OF THE + SYMBOL ADDED AT THE END
      }else{
        test1=substr(test1,1,str_length(test1)-1)
      }  
      
      for(ii in term_to_add){ # THIS ADDS THE TERMS NOT INCLUDED
        test1 = paste(test1,ii,sep="+") # SO WE ALWAYS HAVE THE 9 PREDICTORS
      }
      
      #final_formula = as.formula(paste("Response ~ ",test1))
      final_formula = paste("Response ~ ",test1)
      
      formula_list[index] <- final_formula
      
      # THE LINE BELOW IS FOR TESTING IF THE FORMULA IS WORKING OR NOT
      #gam_summary[[index]]=summary(mgcv::gam(final_formula,data = cbind(Response,SQRPredictors),family = gaussian()))
      
      index=index+1
      # A+B+C == A^2 +B+C == A:A +B+C
    }
  }
  
  if(O2){
    write.csv(formula_list,paste0("O2_formula_list"))
  }else{
    write.csv(formula_list,paste0("O1_formula_list"))
  }
}

{
  O2=F
  if(O2){
    Predictor_names = colnames(SQRPredictors)
  }else{
    Predictor_names = colnames(Predictor)
  }
  
  term_matrix=c();
  
  for (i in 1:length(Predictor_names)) {
    term_matrix[[i]] <- combn(Predictor_names, m = i)
  }
  
  dim(term_matrix[[3]])[2]
  
  
  formula_list=vector()
  index=1
  gam_summary=c()
  
  for (m in 1:length(Predictor_names)) {
    
    j = dim(term_matrix[[m]])[2]
    print(j)
    
    for(j in 1:dim(term_matrix[[m]])[2] ){
      
      test1=c()
      for(i in 1:dim(term_matrix[[m]])[1]){
        
        if(dim(term_matrix[[m]])[1] != 1){
          test1 = paste( term_matrix[[m]][i,j], test1 ,sep="+")
        }else{
          test1 = paste( term_matrix[[m]][i,j],":",term_matrix[[m]][i,j],"+",sep="")
        }
      }  
      #print(Predictor_names[!(Predictor_names %in% term_matrix[[3]][,1])])
      
      term_to_add = Predictor_names[!(Predictor_names %in% term_matrix[[m]][,j])]
      
      if(dim(term_matrix[[m]])[1] != 1){
        test1 = paste("(",substr(test1,1,str_length(test1)-1),")^2",sep="") # str_length(xx) - 1 TAKES CARE OF THE + SYMBOL ADDED AT THE END
      }else{
        test1=substr(test1,1,str_length(test1)-1)
      }  
      
      for(ii in term_to_add){ # THIS ADDS THE TERMS NOT INCLUDED
        test1 = paste(test1,ii,sep="+") # SO WE ALWAYS HAVE THE 9 PREDICTORS
      }
      
      #final_formula = as.formula(paste("Response ~ ",test1))
      final_formula = paste("Response ~ ",test1)
      
      formula_list[index] <- final_formula
      
      # THE LINE BELOW IS FOR TESTING IF THE FORMULA IS WORKING OR NOT
      #gam_summary[[index]]=summary(mgcv::gam(final_formula,data = cbind(Response,SQRPredictors),family = gaussian()))
      
      
      index=index+1
      # A+B+C == A^2 +B+C == A:A +B+C  
    }
  }
  
  
  # sink("Gam_summary.txt")
  # print(gam_summary)
  # sink()
  # 
  # sink("All_formula.txt")
  # print(formula_list)
  # sink()
  
  if(O2){
    write.csv(formula_list,paste0("O2_formula_list"))
  }else{
    write.csv(formula_list,paste0("O1_formula_list"))
  }
}



############# FOR SMOOTH FUNCTION ####

Predictor_names = colnames(Predictor)

term_matrix=c();

for (i in 1:length(Predictor_names)) {
  term_matrix[[i]] <- combn(Predictor_names, m = i)
}

dim(term_matrix[[3]])[2]


formula_list=vector()
index=1
gam_summary=c()

for (m in 1:length(Predictor_names)) {
  
  j = dim(term_matrix[[m]])[2]
  print(j)
  
  for(j in 1:dim(term_matrix[[m]])[2] ){
    
    test1=c()
    test2=c()
    for(i in 1:dim(term_matrix[[m]])[1]){
      #print(term_matrix[[m]][i,j])
      
      if(dim(term_matrix[[m]])[1] != 1){
        test1 = paste( term_matrix[[m]][i,j], test1 ,sep=",")
        #test1 = paste( "s(",term_matrix[[m]][i,j],")+",sep="")
        test2 = paste( "s(",term_matrix[[m]][i,j],") +",test2,sep="")
        #print(test1)
        #print(test2)
      }else{
        test1 = paste( " s(",term_matrix[[m]][i,j],") +",sep="")
        #print(test1)
      }
    }  
    #print(Predictor_names[!(Predictor_names %in% term_matrix[[3]][,1])])
    
    term_to_add = Predictor_names[!(Predictor_names %in% term_matrix[[m]][,j])]
    
    if(dim(term_matrix[[m]])[1] != 1){
      test1 = paste("s(",substr(test1,1,str_length(test1)-1),") ",sep="") # str_length(xx) - 1 TAKES CARE OF THE + SYMBOL ADDED AT THE END
      #print(test1)
      
      test2 = paste(substr(test2,1,str_length(test2)-1),sep="") # str_length(xx) - 1 TAKES CARE OF THE + SYMBOL ADDED AT THE END
      #print(test2)
    }else{
      test1=substr(test1,1,str_length(test1)-1)
      #print(test1)
    }
    
    for(ii in term_to_add){ # THIS ADDS THE TERMS NOT INCLUDED
      test1 = paste(test1,ii,sep=" + ") # SO WE ALWAYS HAVE THE 9 PREDICTORS
      test2 = paste(test2,ii,sep=" + ")
    }
    
    #final_formula = as.formula(paste("Response ~ ",test1))
    #print(test1)
    #print(test2)
    
    
    if(m==1){ # IF M=1 THEN ONLY STORE THE FORMULA IN TEST1 NOT IN TEST2
      final_formula = paste("Response ~ ",test1)
      formula_list[index] <- final_formula
      index=index+1
      #print(final_formula)
    }
    if(m==2){ # IF M==2 THEN STORE BOTH TEST1 AND TEST2 FORMULAS
      final_formula = paste("Response ~ ",test1)
      formula_list[index] <- final_formula
      index=index+1
      #print(final_formula)
      
      final_formula = paste("Response ~ ",test2)
      formula_list[index] <- final_formula
      index=index+1
      #print(final_formula)
    }
    if(m>=3){ # IF M>2 THEN ONLY STORE TEST2 BECAUSE WE DONT WANT SMOOTH OF MORE THAN 2 VARIABLES AT A TIME
      # THAT IS WE DONT WANT s(A,B,C,...). WE ONLY WANT s(A,B) max
      final_formula = paste("Response ~ ",test2)
      formula_list[index] <- final_formula
      index=index+1
      #print(final_formula)
    }
  }
}

#formula_list
write.csv(formula_list,paste0("SmoothedO1_formula_list"))

# source("Find_Best_GAM.R")


###### PAPER-STYLE 100 RANDOM TRAIN/TEST SPLITS WITH PER-SPLIT CHECKPOINTS #####
# Follows Dainotti et al. 2025 Sec. 4.2.3: 100 randomized 80:20 train/test
# splits. Per split, every formula gets one 10fCV pass on the training set;
# the subset passing (r >= 99.9% quantile) AND (RMSE <= 2% quantile) is then
# scored on the held-out test set. Scoring is in LINEAR z-space throughout
# (predictions converted via z = 10^pred - 1 before computing r/RMSE), matching
# the paper's own Fig 6 exactly -- its axes are labeled "Linear scale RMSE" /
# "Linear scale correlation" and its stated cutoffs (r=0.565, RMSE=1.167) are
# only sensible in linear z, not log10(z+1). The GAM/GLM models still fit
# log10(z+1) internally; only the selection metric is converted back to z.
# Bias (Sec 4.4: <z_pred - z_obs>) is a GATE,
# not a vote -- candidates worse than the split's median |bias| are dropped,
# then the best-r and best-rmse formula among the survivors each cast one
# vote. (An earlier version gave bias its own equal-weight vote; diagnostics
# showed that criterion picks a near-different formula on ~90/100 splits --
# a single split's mean signed residual is too noisy a statistic to vote on
# directly -- so it was demoted to a pre-filter instead.) Aggregating votes
# across splits (formula_win_frequency.csv) replaces the paper's "appeared
# the maximum number of times" count. Results are checkpointed after EVERY
# split into checkpoints/, so a killed run keeps all completed splits, and a
# restarted run resumes where it left off instead of redoing finished splits.

Response  = TrainData$log10z
Predictor = subset(TrainData, select = lassovar)
SqrData   = SqrTermGen(Predictor)
SqrData$Response <- TrainData$log10z

O2_formula_list  = read.csv("O2_formula_list", header = T, row.names = 1)
all_formula_list = na.omit(O2_formula_list)
all_formula <- apply(as.matrix(all_formula_list), 1, as.formula)
n_formula   <- length(all_formula)
cat("Formulas:", n_formula, "| GRBs:", nrow(SqrData), "\n")

MyLaptop = F  # supercomputer job: 262,143 formulas (9 vars) x 100 splits is
              # ~16x the standard 7-var search -- not attempted locally.
if (MyLaptop) {
  slaves <- detectCores() - 1
  cl_onenode <- makeCluster(slaves)
  registerDoParallel(cl_onenode)
} else {
  slaves <- 72 - 1
  sink("/dev/null"); cl_onenode <- makeCluster(slaves, type = "MPI"); sink()
  registerDoParallel(cl_onenode)
}

N_SPLITS <- 100
CKPT_DIR <- "checkpoints"
if (!dir.exists(CKPT_DIR)) dir.create(CKPT_DIR)
summary_file <- file.path(CKPT_DIR, "best_formulas_per_split.csv")

# Resume support: skip splits that already have a checkpoint on disk.
done_ids <- as.integer(gsub("[^0-9]", "",
              list.files(CKPT_DIR, pattern = "^split_\\d+\\.rds$")))

for (s in 1:N_SPLITS) {
  if (s %in% done_ids) { cat("split", s, "already done, skipping\n"); next }
  tick <- proc.time()

  set.seed(s)
  test_idx <- sample(nrow(SqrData), size = round(0.2 * nrow(SqrData)))
  TrainSet <- SqrData[-test_idx, ]
  TestSet  <- SqrData[test_idx, ]

  # Everything below scores/selects in LINEAR z-space, not log10(z+1), matching
  # the paper exactly: Fig 6's own axis labels are "Linear scale RMSE" / "Linear
  # scale correlation", and its stated cutoffs (r=0.565, RMSE=1.167) only make
  # sense there -- an RMSE of 1.167 is meaningless for log10(z+1) values, which
  # only range ~0-1 for this sample, but is a perfectly ordinary linear-z RMSE.
  # The GAM/GLM models themselves still fit log10(z+1) (that's the whole point
  # of the log link -- linear residuals in z blow up at high z), but selection
  # converts predictions back to z = 10^pred - 1 before scoring.
  train_z_obs <- 10^TrainSet$Response - 1
  scores <- foreach(j = 1:n_formula, .combine = rbind,
                    .multicombine = TRUE, .maxcombine = 1000,
                    .packages = c("mgcv", "caret")) %dopar% {
    set.seed(s * 1e6 + j)
    fm   <- all_formula[[j]]
    pred <- rep(NA_real_, nrow(TrainSet))
    ok <- tryCatch({
      folds <- createFolds(y = TrainSet$Response, k = 10)
      for (i in seq_along(folds)) {
        g <- mgcv::gam(fm, data = TrainSet[-folds[[i]], ], family = gaussian())
        pred[folds[[i]]] <- predict(g, TrainSet[folds[[i]], ])
      }
      TRUE
    }, error = function(e) FALSE)
    if (!ok || anyNA(pred)) c(NA_real_, NA_real_)
    else {
      z_pred <- 10^pred - 1
      c(cor(z_pred, train_z_obs),
        sqrt(mean((z_pred - train_z_obs)^2)))
    }
  }
  colnames(scores) <- c("r", "rmse")

  # Candidate subset: same cutoffs as the paper (r above 99.9% quantile,
  # RMSE below 2% quantile), now computed on the linear-z scores above.
  # Fall back to top-10 by r if the AND is empty.
  r_cut    <- quantile(scores[, "r"],    0.999, na.rm = TRUE)
  rmse_cut <- quantile(scores[, "rmse"], 0.02,  na.rm = TRUE)
  cand <- which(scores[, "r"] >= r_cut & scores[, "rmse"] <= rmse_cut)
  if (length(cand) == 0) cand <- order(-scores[, "r"])[1:10]

  # bias, per the paper's own definition (Sec 4.4): <z_pred - z_obs>, the mean
  # SIGNED difference in linear z-space -- not an absolute-value metric, and
  # not the same as RMSE/MAD (which are magnitude-only). r/rmse here are also
  # linear-z now, consistent with the training-CV scores above.
  z_obs_test <- 10^TestSet$Response - 1
  test_eval <- t(sapply(cand, function(j) {
    tryCatch({
      g <- mgcv::gam(all_formula[[j]], data = TrainSet, family = gaussian())
      p <- predict(g, TestSet)
      z_pred <- 10^p - 1
      c(r    = cor(z_pred, z_obs_test),
        rmse = sqrt(mean((z_pred - z_obs_test)^2)),
        bias = mean(z_pred - z_obs_test))
    }, error = function(e) c(r = NA_real_, rmse = NA_real_, bias = NA_real_))
  }))

  # Bias is a GATE, not a vote: a single split's mean signed residual is a
  # noisy statistic (57ish held-out GRBs; +/- errors cancel by luck), and
  # letting it cast an equal-weight best-|bias| vote alongside r/rmse just
  # injects that noise into the tally -- diagnosed empirically by checking how
  # many DISTINCT formulas ever won each per-split criterion across the 100
  # splits: r and rmse concentrate on ~60-70 consistently-good formulas, while
  # best-|bias| scattered across ~90/100 splits (i.e. a different, often much
  # worse-fitting, formula "won" bias almost every single time just by chance
  # cancellation). So: drop the worse-|bias| half of this split's candidates
  # first, THEN vote best-r / best-rmse only among the survivors -- this still
  # screens out systematically-offset formulas without rewarding lucky bias.
  bias_cut <- median(abs(test_eval[, "bias"]), na.rm = TRUE)
  surv     <- cand[abs(test_eval[, "bias"]) <= bias_cut]
  surv_idx <- match(surv, cand)  # positions of survivors within test_eval's row order

  best_r    <- surv[which.max(test_eval[surv_idx, "r"])]
  best_rmse <- surv[which.min(test_eval[surv_idx, "rmse"])]
  best_bias <- cand[which.min(abs(test_eval[, "bias"]))]  # diagnostic only, not a vote

  saveRDS(list(split = s, test_idx = test_idx, scores = scores,
               candidates = cand, test_eval = test_eval, bias_gate_cut = bias_cut,
               best = c(r = best_r, rmse = best_rmse),
               best_bias_diag = best_bias),
          file.path(CKPT_DIR, sprintf("split_%03d.rds", s)))

  row <- data.frame(
    split          = s,
    best_r_idx     = best_r,
    best_r_val     = round(max(test_eval[surv_idx, "r"],    na.rm = TRUE), 4),
    best_rmse_idx  = best_rmse,
    best_rmse_val  = round(min(test_eval[surv_idx, "rmse"], na.rm = TRUE), 4),
    best_bias_idx  = best_bias,
    best_bias_val  = round(test_eval[cand == best_bias, "bias"], 4),
    n_candidates   = length(cand),
    n_survivors    = length(surv),
    minutes        = round((proc.time() - tick)[3] / 60, 1),
    best_r_formula = as.character(all_formula_list[best_r, 1]))
  write.table(row, summary_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(summary_file),
              append = file.exists(summary_file))

  cat(sprintf("split %d/%d done in %.1f min | best test r=%.3f\n",
              s, N_SPLITS, (proc.time() - tick)[3] / 60,
              max(test_eval[, "r"], na.rm = TRUE)))
}

# Final tally over every completed split (works on partial runs too --
# rerun just this block to summarize whatever checkpoints exist).
files <- list.files(CKPT_DIR, pattern = "^split_\\d+\\.rds$", full.names = TRUE)
wins  <- unlist(lapply(files, function(f) readRDS(f)$best))
tab   <- sort(table(wins), decreasing = TRUE)
freq  <- data.frame(formula_idx = as.integer(names(tab)),
                    wins        = as.integer(tab),
                    formula     = as.character(all_formula_list[as.integer(names(tab)), 1]))
write.csv(freq, "formula_win_frequency.csv", row.names = FALSE)
cat("Splits completed:", length(files),
    "| top formulas written to formula_win_frequency.csv\n")
print(head(freq, 10))

stopCluster(cl_onenode)


############# M-ESTIMATOR OUTLIER CUT USING THE MOST-FREQUENT WINNER #############
# Dainotti et al. 2025 Sec. 4.3: "we selected the formula [...] which obtained
# the highest counts" and used THAT formula (not a generic regression) to fit
# the M-estimator on the full sample, dropping the bottom 5% by weight. Only
# run this once every split above has completed -- formula_win_frequency.csv
# is only final at that point.
if (length(files) == N_SPLITS) {
  Formula_for_outlier <- all_formula[[freq$formula_idx[1]]]
  cat("Most-frequent winning formula (", freq$wins[1], "/", N_SPLITS, "splits):\n")
  print(Formula_for_outlier)
  writeLines(deparse(Formula_for_outlier), "Formula_for_outlier.txt")

  require(MASS)
  M_est <- MASS::rlm(Formula_for_outlier, data = SqrData, method = "M", maxit = 50)
  weights <- M_est$w
  weight_threshold <- quantile(weights, 0.05)

  kept_rows    <- rownames(SqrData)[weights > weight_threshold]
  removed_rows <- rownames(SqrData)[weights <= weight_threshold]

  cat("M-estimator outlier cut:", length(removed_rows), "of", nrow(SqrData),
      "GRBs removed (weight <=", round(weight_threshold, 4), ")\n")
  writeLines(removed_rows, "removed_outliers.txt")

  final_cut <- raw_xray_data[rownames(raw_xray_data) %in% kept_rows, ]
  write.csv(final_cut, "final_outliers_removed.csv")
  cat("Final outlier-removed data saved:", nrow(final_cut), "GRBs -> final_outliers_removed.csv\n")
} else {
  cat("Only", length(files), "/", N_SPLITS, "splits completed -- skipping the",
      "M-estimator cut until the full 100-split search finishes",
      "(rerun this script; completed splits are checkpointed and skipped).\n")
}
