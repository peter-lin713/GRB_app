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
  #raw_xray_data <- read.csv("SORTED_FINAL_X-Ray_DATA.csv", header = TRUE, row.names = 1)
  raw_xray_data <- read.csv("combined_data_with_redshift.csv", header = TRUE, row.names = 1)
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
features_for_mice_preds = subset(raw_xray_data,select = c(log10T90,
                                                          log10Fa,
                                                          log10Ta,
                                                          Alpha,
                                                          Beta,
                                                          Gamma,
                                                          log10Fluence,
                                                          PhotonIndex,
                                                          log10NH,
                                                          log10PeakFlux))

features_for_mice_errs = subset(raw_xray_data,select = c(T90Err,
                                                         log10FaErr,
                                                         log10TaErr,
                                                         AlphaErr,
                                                         BetaErr,
                                                         FluenceErr,
                                                         PhotonIndexErr,
                                                         PeakFluxErr))

# replacing inf in log10PeakFlux feature with NAs
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

# removing all log10NH values lower than 20
features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20] <- NA

features_for_mice_preds$Beta[features_for_mice_preds$Beta > 3] <- NA

features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3] <- NA

features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3] <- NA

features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA


if(do_mice){
  set.seed(1)
  
  # features_for_mice_all <- cbind(features_for_mice_preds, features_for_mice_errs)
  mice_model_preds <- mice(data = features_for_mice_preds,
                           m = 20,
                           method = 'midastouch',
                           printFlag = F)
  features_for_mice_preds <- complete(mice_model_preds,20)
  
  mice_model_errs <- mice(data = features_for_mice_errs,
                          m = 20,
                          method = 'midastouch',
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
lassovar=head(lassovar,7)  # matches the paper: "the top seven features... are picked"


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


###### TRAIN TEST SPLIT #####


Response = TrainData$log10z
Predictor = subset(TrainData,select = lassovar)

O2Predictor = quadTermGen(Predictor)
SQRPredictor = SqrTermGen(Predictor)


SqrTrainData = SqrTermGen(Predictor)

SqrTrainData$Response <- TrainData$log10z

GamTrainData <- subset(SqrTrainData, select = !colnames(SqrTrainData) %in% c("log10z"))

GamTestData <- head(GamTrainData,n = (0.05*nrow(GamTrainData)))
dim(GamTestData)

GamValidationData <- tail(GamTrainData,n = 0.05*nrow(GamTrainData) )
dim(GamValidationData)

GamTrainData <- GamTrainData[!(rownames(GamTrainData)%in%c(rownames(GamTestData),rownames(GamValidationData))),]

intersect(rownames(GamTestData),rownames(GamTrainData))
intersect(rownames(GamTestData),rownames(GamValidationData))

######### SETTING UP GAM FUNCTION ######


bestGAM <- as.formula("Response ~ (logPeak + Beta)^2 + (PhotonIndex + logNH)^2 + (PhotonIndex + Beta)^2 + (logNH + logFa)^2 + (logT_a + Alpha)^2")

###### SNIPPET TO READ IN THE FORMULA FILE #####
O1_formula_list = read.csv("O1_formula_list",header = T,row.names = 1)
#O1_formula_list = read.csv(paste0("O1_formula_list_",ncol(Predictor),"features.csv"),row.names = 1)

#O1_gam_formulas = apply(as.matrix(O1_formula_list), 1, as.formula)

O2_formula_list = read.csv("O2_formula_list",header = T,row.names = 1)
#O2_gam_formulas = apply(as.matrix(O2_formula_list), 1, as.formula)

SO1_formula_list = read.csv("SmoothedO1_formula_list",header = T,row.names = 1)
#SO1_gam_formulas = apply(as.matrix(SO1_formula_list), 1, as.formula)
###############################################

# ########## TESTING SUPERLEARNER #############
# system.time({gam_sl = SuperLearner(Y = Response,
#                        X = SQRPredictor,
#                        family = gaussian(),
#                        SL.library = c(learners2$names),
#                        control = list(saveFitLibrary=T),
#                        verbose = F
#                        ,obsWeights = c(1:nrow(Predictors))
# )
# })

################ TESTING MGCV::GAM ###########

all_formula_list <- rbind(O1_formula_list
                          #,SO1_formula_list
                          ,O2_formula_list
)

#O1andSO1=dim(O1_formula_list)[1] + dim(SO1_formula_list)[1]

#best_O1_SO1 = as.formula("Response ~ (logPeak + logFluence + PhotonIndex + Alpha)^2 + logFa + logT_a + Beta + logT90 + logNH")

#best_O2_run1 = as.formula("Response ~ (logNHSqr + BetaSqr + logPeak + logNH + logT90 + Beta)^2 + logFa + logT_a + Alpha + PhotonIndex + logFluence + logFaSqr + logT_aSqr + AlphaSqr + PhotonIndexSqr + logT90Sqr + logFluenceSqr + logPeakSqr")

#best_O2_run2 = as.formula("Response ~ (logPeakSqr + logFluenceSqr + logT90Sqr + BetaSqr + logFaSqr + logPeak + logT90)^2 + logFa + logT_a + Alpha + Beta + PhotonIndex + logNH + logFluence + logT_aSqr + AlphaSqr + PhotonIndexSqr + logNHSqr")

# WE ARE KEEPING FORMULAS TILL 262450 BECAUSE AFTER THIS ERRORS OCCUR DUE TO TOO MANY VARIABLES
#all_formula_list <- all_formula_list[1:max_formula,]
all_formula_list = na.omit(all_formula_list)
#all_formula <- c(bestGAM,O1_gam_formulas,SO1_gam_formulas,O2_gam_formulas[18:10000])

all_formula <- apply(as.matrix(all_formula_list), 1, as.formula)



# formula_analysis_addr = paste0("Formula_Generation/GAM_",Sys.Date(),"/")
# if(!dir.exists(formula_analysis_addr)){dir.create(formula_analysis_addr)}


GamValidationData = rbind(GamValidationData,GamTestData)
#write.csv(GamTestData,file = "MGCV/Gam_test_set.csv",row.names = rownames(GamTestData))
write.csv(GamValidationData,file = paste0("Gam_validation_set.csv"),row.names = rownames(GamValidationData))
write.csv(GamTrainData,file = paste0("Gam_train_set.csv"),row.names = rownames(GamTrainData))

GamTrainData_unscaled = GamTrainData
GamValidationData_unscaled = GamValidationData

#### STORE THE SD AND MEAN OF THE GAM TRAIN DATA
# THEN SCALE THE VALIDATION SET WITH THAT SD AND MEAN
GamTrainData_sd = apply(GamTrainData,2,sd)
GamTrainData_mean = apply(GamTrainData,2,mean)


#print(max_formula)
#print(max_parallel)
#print(num_cores)

# print(paste(
#   "Number of formulas=",length(all_formula)
#   ,"| Dimension of training set=",dim(GamTrainData)
# )
#   #,"Scaling the data = ",Scaling)
# )


MyLaptop = F

if(MyLaptop){
  
  slaves <- detectCores()-2
  { #sink("/dev/null");
    cl_onenode <- makeCluster(slaves);
    #sink();
  } # number of MPI tasks to use
  registerDoParallel(cl_onenode)
  
}else{  
  
  slaves <- 72 - 1#detectCores() - 1
  { 
    sink("/dev/null"); 
    cl_onenode <- makeCluster(slaves, type="MPI"); 
    sink(); 
  } # number of MPI tasks to use
  registerDoParallel(cl_onenode)
}



#sink("MGCV/progress.txt",append=T)
tick <- proc.time()
manual_cv <- foreach(j = 1:length(all_formula)
                     #,.export = c(all_formula)
                     ,.packages=c("SuperLearner","mgcv", "caret" ,"xgboost", "randomForest", "gbm", "lattice", "Matrix", "glmnet", "biglasso","e1071",'earth','party')
) %dopar% {
  
  InnerLoop = 50
  
  test_preds <- data.frame(Predicted= numeric(nrow(GamTrainData)),Observed= numeric(nrow(GamTrainData)))
  test_preds_loop= matrix(nrow = nrow(GamTrainData),ncol = InnerLoop)
  
  set.seed(j)
  
  gam_formula = all_formula[[j]]
  
  for (k in 1:InnerLoop) {
    
    folds <- createFolds(y = GamTrainData$Response, k = 10)
    
    ###### THE 10FCV SECTION ##########
    for (i in 1:length(folds)) {
      train_set <- GamTrainData[-c(folds[[i]]),]
      test_set <- GamTrainData[c(folds[[i]]),]
      
      #print(rownames(head(train_set)))
      #print(rownames(test_set))
      
      #gam_model <- MASS::rlm(gam_formula,train_set,method = 'M')
      gam_model <- mgcv::gam(formula = gam_formula
                             ,data = train_set
                             ,family = gaussian())
      
      
      ###### TEST SET PREDICTION ###########
      #predict(gam_model,test_set)
      test_preds$Predicted[c(folds[[i]])] <- predict(gam_model,test_set)
      test_preds$Observed[c(folds[[i]])] <- test_set$Response
    }
    
    test_preds_loop[,k] = test_preds$Predicted
    
  }
  
  test_preds$Predicted = rowMeans(test_preds_loop)
  
  # THE ACTUAL RETURN STATEMENT
  gam_model <- mgcv::gam(formula = gam_formula
                         ,data = GamTrainData
                         ,family = gaussian())
  #gam_model <- MASS::rlm(gam_formula,train_set,method = 'M')
  return(list(test_preds$Predicted,predict(gam_model,GamValidationData)))
}

tock <- proc.time() - tick

# NEED TO TEST THIS TIMING
cat(tock)

#sink();
# CREATE DIRECTORIES IF THEY DONT EXIST
saveRDS(manual_cv,file = paste0("SuperLearner_complete_data.rds"))
saveRDS(all_formula,file=paste0("Formulas_used.rds"))


############# SELECT WINNING FORMULA + M-ESTIMATOR OUTLIER CUT #############
# Following Dainotti et al. 2025 (Sec. 4.2.3-4.3): the formula search above ran
# on the FULL, uncut sample. Now we pick the single best-performing formula
# (same selection logic as "Formula_for_outlier" in Formula_Search_Aditya_v2.R:
# CV correlation/RMSE quantile cutoffs, then best validation-set correlation
# among survivors) and use THAT formula -- not a generic regression -- to fit
# the M-estimator on the full sample. Only then do we drop the bottom 5% by
# weight. This must run after the full formula search, since the winning
# formula is not known in advance.

CV_Prediction_matrix         <- matrix(nrow = length(manual_cv[[1]][[1]]), ncol = length(all_formula))
Validation_Prediction_matrix <- matrix(nrow = length(manual_cv[[1]][[2]]), ncol = length(all_formula))

CV_correlation          <- vector(length = length(all_formula))
CV_RMSE                 <- vector(length = length(all_formula))
Validation_correlations <- vector(length = length(all_formula))

for (j in seq_along(all_formula)) {
  CV_Prediction_matrix[, j] <- manual_cv[[j]][[1]]
  CV_correlation[j] <- cor(CV_Prediction_matrix[, j], GamTrainData$Response)
  CV_RMSE[j]        <- sqrt(mean((CV_Prediction_matrix[, j] - GamTrainData$Response)^2))

  Validation_Prediction_matrix[, j] <- manual_cv[[j]][[2]]
  Validation_correlations[j] <- cor(Validation_Prediction_matrix[, j], GamValidationData$Response)
}

correlation_cutoff <- 0.999   # same cutoffs used elsewhere in this pipeline
RMSE_cutoff        <- 0.02    # (see runs/formula_search_*_paper_cutoffs)

CV_correlation_cutoff <- quantile(CV_correlation, correlation_cutoff)
CV_RMSE_cutoff        <- quantile(CV_RMSE, RMSE_cutoff)

Correlation_formula <- which(CV_RMSE < CV_RMSE_cutoff & CV_correlation > CV_correlation_cutoff)
if (length(Correlation_formula) == 0) {
  # Cutoffs too strict to pass anything (can happen with a small/capped
  # formula list) -- fall back to the single best-CV-correlation formula.
  Correlation_formula <- which.max(CV_correlation)
}

best_idx <- Correlation_formula[which.max(Validation_correlations[Correlation_formula])]
Formula_for_outlier <- all_formula[[best_idx]]

cat("Winning formula selected for outlier removal:\n")
print(Formula_for_outlier)
writeLines(deparse(Formula_for_outlier), "Formula_for_outlier.txt")

# ---- M-estimator outlier cut using the winning formula, on the FULL sample ----
# GamTrainData (90%) + GamValidationData (already rbind'd with GamTestData,
# 10%) together reconstruct the full sample the formula search started from.
FullData <- rbind(GamTrainData, GamValidationData)

require(MASS)
M_est <- MASS::rlm(Formula_for_outlier, data = FullData, method = "M", maxit = 50)
weights <- M_est$w
weight_threshold <- quantile(weights, 0.05)

kept_rows    <- rownames(FullData)[weights > weight_threshold]
removed_rows <- rownames(FullData)[weights <= weight_threshold]

cat("M-estimator outlier cut:", length(removed_rows), "of", nrow(FullData),
    "GRBs removed (weight <=", round(weight_threshold, 4), ")\n")

writeLines(removed_rows, "removed_outliers.txt")

final_cut <- raw_xray_data[rownames(raw_xray_data) %in% kept_rows, ]
write.csv(final_cut, "final_outliers_removed.csv")

cat("Final outlier-removed data saved:", nrow(final_cut), "GRBs -> final_outliers_removed.csv\n")


