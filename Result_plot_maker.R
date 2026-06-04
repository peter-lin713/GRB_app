#' Result_plot_maker.R — cross-validation result statistics and plots.
#'
#' Defines result_plotter(), the reporting stage of the pipeline. Given the
#' cross-validated predictions and their observed truths (both on the
#' log10(z+1) scale), it computes the standard redshift-estimation diagnostics
#' (bias, sigma, RMS, MAD/NMAD, the 1- and 2-sigma "cones") and writes a fixed
#' set of CSVs and PNGs, then returns the per-GRB results table.
#'
#' Scale conventions:
#'   InvZspec / InvZphot   observed / predicted log10(z+1)  ("inverse"-named for
#'                         historical reasons; these are the log-scale values)
#'   Zspec / Zphot         observed / predicted linear redshift z = 10^log - 1
#'   Dlogz, Dz, normDz     residuals on the log, linear, and normalized scales
#'
#' Globals expected from the calling pipeline:
#'   addr       output directory prefix for CSV / .Rdata files
#'   PLOTaddr   output directory prefix for PNGs
#'   plotnames  suffix tag appended to every output filename
#'   sz, rez    plot size multiplier and resolution (also set locally below)
#'
#' Outputs (filenames suffixed by `plotnames`):
#'   10fCVResults*.csv, Results_wo_catout*.csv          result tables
#'   BoxPlotOfRedshift.png, BoxPlotOfNormalizedRedshift.png
#'   DeltaZSpreadBLL.png, NormDeltaZSpreadBLL.png         residual histograms
#'   z_pred_v_obs_log*.png, z_pred_v_obs_linear*.png      pred-vs-obs scatters
#'   For_proposal_z_pred_v_obs_{log,linear}*.png          publication variants
#'   Workspace_outputs*.Rdata                             full workspace dump
#'
#' @param names character. GRB identifiers; become the result rownames.
#' @param p numeric. Predicted log10(z+1) per GRB.
#' @param o numeric. Observed log10(z+1) per GRB.
#' @param pred.max numeric. Per-GRB upper error-bar bound, log scale (default 0).
#' @param pred.min numeric. Per-GRB lower error-bar bound, log scale (default 0).
#' @param pred.sd numeric. Per-GRB prediction SD, log scale (default 0).
#' @param linpred.sd numeric. Per-GRB prediction SD, linear z scale (default 0).
#' @return data.frame. Per-GRB results table (also written to CSV).

###### STATISTICAL PARAMETERS ####
result_plotter <- function(names,p,o,pred.max=0,pred.min=0,pred.sd = 0, linpred.sd = 0){

# Colour-blind-friendly palette used for the scatter plots.
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Assemble the results table: predicted/observed log10(z+1), then derive the
# linear-z columns and the residuals on each scale.
results<-data.frame(InvZphot=round(p,4)
                    ,InvZspec=round(o,4)
                    )

rownames(results) <- names

results$Dlogz <- (results$InvZspec - results$InvZphot)
results$Zphot  <- 10^results$InvZphot - 1
results$Zspec  <- 10^results$InvZspec - 1
results$Dz     <- (results$Zspec - results$Zphot) # DIFFERENCE BETWEEN OBS AND PRED Z
results$normDz <- results$Dz/(1+results$Zspec) # NORMALIZED RESIDUALS
results$SD <- pred.sd # STANDARD DEVIATIONS (log scale)
results$z_SD <- linpred.sd # STANDARD DEVIATION IN Z (linear scale)
results$pred_max <- pred.max
results$pred_min <- pred.min

results$linpred_max <- 10^(pred.max) - 1
results$linpred_min <- 10^(pred.min) - 1

write.csv(results,file = paste(addr,'10fCVResults',plotnames,'.csv',sep = ''))

#### CALCULATING SIGMA IN LOG SCALE ###

L_Sigma <- sd(results$Dlogz) # STANDARD DEVIATION IN LOGSCALE
L_2Sigma <-2*L_Sigma

normRMS<-sqrt(mean((results$normDz)^2))
InvRMS<-sqrt(mean((results$Dlogz)^2))

InvBias<-mean(results$Dlogz)
Bias<-mean(results$Dz) #BIAS CALCULATED FROM DELTA Z
normBias<-mean(results$normDz) # NORMALIZED BIAS CALCULATED FROM NORMALIZED DELTA Z

plot(results$Zspec,results$Dz,pch='*',main='Redshift vs Bias',xlab='z',ylab='Bias')
abline(h=Bias)

plot(results$Zspec,results$normDz,main='Redshift vs Normalized Bias',xlab='z',ylab='Normalized Bias')
abline(h=normBias)

Sigma <- sqrt(sum((results$Dz-Bias)^2)/(nrow(results))) # SIGMA = RMS OF (DELTA Z - BIAS)

normSigma<-sqrt(mean((results$normDz-normBias)^2 )) # NORM SIGMA = RMS OF (NORM DELTA Z - NORM BIAS)

MAD<-median(abs(results$Dz))
NMAD<-median(abs(results$normDz))

CVboxplot<-boxplot(results$InvZphot-results$InvZspec) # BOXPLOT IN Log(z+1)


{ ##### BOX PLOT: residuals in linear z ####
  png(filename = paste(PLOTaddr,'BoxPlotOfRedshift', '.png', sep = ''),width = (750*sz),height = (750*sz),res=rez)
  CVboxplot<-boxplot(results$Dz) # BOXPLOT IN z
  boxplot(results$Dz, ylab = "Dz = observed z - predicted z",
          cex.main = 0.9,
          main = paste0('Redshift residuals  (n=', CVboxplot$n, ')\n',
                        'Dz = observed z - predicted z\n',
                        'median=', signif(CVboxplot$stats[3,], 2),
                        ', points beyond whiskers=', length(CVboxplot$out))
  )
  dev.off()
}

{ ##### BOX PLOT: normalized residuals ####
  png(filename = paste(PLOTaddr,'BoxPlotOfNormalizedRedshift','.png', sep = ''),width = (750*sz),height = (750*sz),res=rez)
  CVboxplot<-boxplot(results$normDz) # BOXPLOT IN z
  boxplot(results$normDz, ylab = "normalized Dz = (obs z - pred z)/(1 + obs z)",
          cex.main = 0.9,
          main = paste0('Normalized redshift residuals  (n=', CVboxplot$n, ')\n',
                        'normalized Dz = (obs z - pred z)/(1 + obs z)\n',
                        'median=', signif(CVboxplot$stats[3,], 2),
                        ', points beyond whiskers=', length(CVboxplot$out))
  )
  dev.off()
}

{ ##### RESIDUAL HISTOGRAM: linear z ####
  png(filename = paste(PLOTaddr,'DeltaZSpreadBLL', '.png',sep = ''),width = (1000*sz),height = (750*sz),res=rez)
  hist(results$Dz,breaks = 50,cex.main = 0.9,main = paste0('Redshift residual distribution  (Dz = obs z - pred z)\n','Sigma=',signif(Sigma,3),', Bias=',signif(Bias,3),'\n(blue = +/-Sigma, red = Bias)')
       ,xlab='Dz = observed z - predicted z'
  )
  abline(v=c(-Sigma,Sigma),col='blue')
  abline(v=Bias,col='red')
  dev.off()
  hist(results$Dz,breaks = 50,cex.main = 0.9,main = paste0('Redshift residual distribution  (Dz = obs z - pred z)\n','Sigma=',signif(Sigma,3),', Bias=',signif(Bias,3),'\n(blue = +/-Sigma, red = Bias)')
       ,xlab='Dz = observed z - predicted z'
  )
  abline(v=c(-Sigma,Sigma),col='blue')
  abline(v=Bias,col='red')

}

{ ##### RESIDUAL HISTOGRAM: normalized ####
  png(filename = paste(PLOTaddr,'NormDeltaZSpreadBLL', '.png',sep = ''),width = (1000*sz),height = (750*sz),res=rez)
  hist(results$normDz,breaks = 25,cex.main = 0.9,main = paste0('Normalized residual distribution  (normalized Dz = (obs z - pred z)/(1+obs z))\n','Sigma=',signif(normSigma,3),', Bias=',signif(normBias,3),'\n(blue = +/-Sigma, red = Bias)')
       ,xlab='normalized Dz = (obs z - pred z)/(1 + obs z)'
  )
  abline(v=c(-normSigma,normSigma),col='blue')
  abline(v=normBias,col='red')
  dev.off()
  hist(results$normDz,breaks = 25,cex.main = 0.9,main = paste0('Normalized residual distribution  (normalized Dz = (obs z - pred z)/(1+obs z))\n','Sigma=',signif(normSigma,3),', Bias=',signif(normBias,3),'\n(blue = +/-Sigma, red = Bias)')
       ,xlab='normalized Dz = (obs z - pred z)/(1 + obs z)'
  )
  abline(v=c(-normSigma,normSigma),col='blue')
  abline(v=normBias,col='red')

}


##### SIGMA CONE CALCULATIONS #####
# A point is "inside the 2-sigma cone" if its log-scale residual stays within
# +/-2*L_Sigma of the identity line; likewise for the 1-sigma cone.

belowUpperLine<-results[results$InvZphot-(1)*results$InvZspec < 2*L_Sigma,]
aboveLowerLine<-results[results$InvZphot-(1)*results$InvZspec > -2*L_Sigma,]

InsideTheCone<-results[intersect(rownames(belowUpperLine),rownames(aboveLowerLine)),]

print(rownames(InsideTheCone))
write.csv(InsideTheCone,file = paste(addr,'Results_wo_catout',plotnames,'.csv',sep = ''))

sigma1_belowUpperLine<-results[results$InvZphot-(1)*results$InvZspec < L_Sigma,]
sigma1_aboveLowerLine<-results[results$InvZphot-(1)*results$InvZspec > -L_Sigma,]

sigma_1InsideTheCone<-results[intersect(rownames(sigma1_belowUpperLine),rownames(sigma1_aboveLowerLine)),]


######### Cross-validation correlation plot: log10(z+1) ####

sz<-0.8
rez=120

uplim<-max(results$InvZphot,results$InvZspec)
lowlim<-min(results$InvZphot,results$InvZspec)

print(plotnames)

{
  png(filename = paste(PLOTaddr,"z_pred_v_obs_log",plotnames,".png",sep = ""),width = 1000*sz,height = 1000*sz,res=rez)
  plt1<-ggplot(results, aes(x=InvZspec, y=InvZphot,col='blue'))+ xlim(lowlim,uplim)+ylim(lowlim,uplim)+
    geom_point(shape=1)+
    geom_point(aes(x=InvZspec,y=InvZphot,col='red'),data = InsideTheCone, shape=1)+
    geom_abline(slope=1,intercept = 0,col="red")+
    geom_abline(slope=1,intercept = 2*L_Sigma,col="blue")+
    geom_abline(slope=1,intercept = -2*L_Sigma,col="blue")+
    geom_errorbar(ymin=results$pred_min,ymax=results$pred_max)+
    labs(color='Points',
         title = paste('Cross-validated predicted vs. observed log10(z+1)\nSample = ',nrow(results),' GRBs  |  within 2-sigma cone = ',nrow(InsideTheCone),' (',signif(100*(nrow(InsideTheCone)/nrow(results)),2),'%)',sep=""),
         subtitle = paste("r = ", signif(cor(results$InvZphot,results$InvZspec),4),
                          '| Sigma = ',signif(L_Sigma,3),
                          "| RMS = ", signif(sqrt(mean((results$Dlogz)^2)),3), # RMS VALUE
                          "| Bias = ", signif(InvBias,2), # BIAS VALUE
                          "| NMAD =",signif(1.48*median(abs(results$Dlogz)),3)
         )) +
    scale_color_manual(values=cbPalette, labels=c('All GRBs','Within 2-sigma cone'))+theme_bw()+theme(plot.background = element_rect(color = 'black',size=1))+
    xlab("Observed log10(z+1)")+ylab("Predicted log10(z+1)  (CV mean; bars = min-max over folds)")

  print(plt1)

  dev.off()
  print(plt1)
}


uplim<-max(results$Zphot,results$Zspec)
lowlim<-min(results$Zphot,results$Zspec)

Yuplim<-max(results$linpred_max)

{
  png(filename = paste(PLOTaddr,"z_pred_v_obs_linear",plotnames,".png",sep = ""),width = 1000*sz,height = 1000*sz,res=rez)
  plt1<-ggplot(results, aes(x=Zspec, y=Zphot),col='black')+xlim(0,8.5)+ylim(0,11)+

    geom_point(shape=16)+

    geom_point(aes(x=Zspec,y=Zphot),col='black',data = InsideTheCone, shape=16)+

    geom_abline(slope=1,intercept = 0,col="red")+

    geom_errorbar(ymin=results$linpred_min,ymax=results$linpred_max,width=0.2,col='black')+

    geom_abline(slope=10^(2*L_Sigma),intercept = (10^(2*L_Sigma) - 1),col="blue")+
    geom_abline(slope=10^(-2*L_Sigma),intercept =(10^(-2*L_Sigma) - 1),col="blue")+
    geom_abline(slope=10^(L_Sigma),intercept = (10^(L_Sigma) - 1),col="green")+
    geom_abline(slope=10^(-L_Sigma),intercept =(10^(-L_Sigma) - 1),col="green")+

labs(color = '',
     title = paste('Cross-validated predicted vs. observed redshift z\nSample = ',nrow(results)
                   ,' GRBs  | within 2-sigma cone = ',nrow(InsideTheCone),' (',signif(100*(nrow(InsideTheCone)/nrow(results)),2),'%)'
                   ,'\nwithin 1-sigma cone = ',nrow(sigma_1InsideTheCone),' (',signif(100*(nrow(sigma_1InsideTheCone)/nrow(results)),2),'%)'
                   ,sep=""),
     subtitle = paste(" r = ", signif(cor(results$Zphot,results$Zspec),3),
                      " | Sigma = ",signif(Sigma,3),
                      " | RMS = ", signif(sqrt(mean((results$Dz)^2)),2), # RMS VALUE
                      " | Bias = ", signif(Bias,2), # BIAS VALUE
                     " | NMAD = ", signif((1.48*mad(results$Dz)),3) # NMAD VALUES FOR normalized residuals
     )) +
    scale_color_manual(values=c('black','black'))+
    theme_bw()+
    scale_x_continuous(limits = c(lowlim,uplim),breaks = seq(round(lowlim),round(uplim),1))+
    scale_y_continuous(limits = c(lowlim,max(results$linpred_max)),breaks = seq(round(lowlim),round(max(results$linpred_max)),1))+
    theme(plot.background = element_rect(color = 'white',size=0,fill=NA)
          ,panel.grid = element_line(colour = 'white')
          ,panel.border = element_rect(colour = 'black',fill=NA,size=0.5)
          ,axis.line = element_line(colour = 'black',size=1)
          ,axis.text.y = element_text( face ="bold",size=10,colour = 'black')
          ,axis.text.x = element_text( face ="bold",size=10,colour = 'black')
          ,axis.ticks = element_line(colour = 'black',size=1))+
    xlab("Observed redshift z")+ylab("Predicted redshift z  (bars = min-max over CV folds)")

  print(plt1)
  dev.off()
  print(plt1)
}


######### Cross-validation correlation plot: publication variants ####

sz<-0.8
rez=200

uplim<-max(results$InvZphot,results$InvZspec)
lowlim<-min(results$InvZphot,results$InvZspec)

print(plotnames)

{
  png(filename = paste(PLOTaddr,"For_proposal_z_pred_v_obs_log",plotnames,".png",sep = ""),width = 1000*sz,height = 1000*sz,res=rez)
  plt1<-ggplot(results, aes(x=InvZspec, y=InvZphot,col='blue'))+ xlim(lowlim,uplim)+ylim(lowlim,uplim)+
    geom_point(shape=1)+
    geom_point(aes(x=InvZspec,y=InvZphot,col='red'),data = InsideTheCone, shape=1)+
    geom_abline(slope=1,intercept = 0,col="red")+
    geom_abline(slope=1,intercept = 2*L_Sigma,col="blue")+
    geom_abline(slope=1,intercept = -2*L_Sigma,col="blue")+
    geom_errorbar(ymin=results$pred_min,ymax=results$pred_max)+
    scale_color_manual(values=cbPalette)+theme_bw()+theme(plot.background = element_rect(color = 'black',size=1))+
    ggtitle("Predicted vs. observed log10(z+1)\n(cross-validated)")+
    xlab("Observed log10(z+1)")+ylab("Predicted log10(z+1)")

  print(plt1)

  dev.off()
  print(plt1)
}


uplim<-max(results$Zphot,results$Zspec)
lowlim<-min(results$Zphot,results$Zspec)

Yuplim<-max(results$linpred_max)

{
  png(filename = paste(PLOTaddr,"For_proposal_z_pred_v_obs_linear",plotnames,".png",sep = ""),width = 1000*sz,height = 1000*sz,res=rez)
  plt1<-ggplot(results, aes(x=Zspec, y=Zphot),col='black')+xlim(0,8.5)+ylim(0,11)+

    geom_point(shape=16)+

    geom_point(aes(x=Zspec,y=Zphot),col='black',data = InsideTheCone, shape=16)+

  geom_abline(slope=1,intercept = 0,col="red")+

    geom_errorbar(ymin=results$linpred_min,ymax=results$linpred_max,width=0.2,col='black')+

    geom_abline(slope=10^(2*L_Sigma),intercept = (10^(2*L_Sigma) - 1),col="blue")+
    geom_abline(slope=10^(-2*L_Sigma),intercept =(10^(-2*L_Sigma) - 1),col="blue")+
    geom_abline(slope=10^(L_Sigma),intercept = (10^(L_Sigma) - 1),col="green")+
    geom_abline(slope=10^(-L_Sigma),intercept =(10^(-L_Sigma) - 1),col="green")+

    scale_color_manual(values=c('black','black'))+
    theme_bw()+
    scale_x_continuous(limits = c(lowlim,uplim),breaks = seq(round(lowlim),round(uplim),1))+
    scale_y_continuous(limits = c(lowlim,max(results$linpred_max)),breaks = seq(round(lowlim),round(max(results$linpred_max)),1))+
    theme(plot.background = element_rect(color = 'white',size=0,fill=NA)
          ,panel.grid = element_line(colour = 'white')
          ,panel.border = element_rect(colour = 'black',fill=NA,size=0.5)
          ,axis.line = element_line(colour = 'black',size=1)
          ,axis.text.y = element_text( face ="bold",size=10,colour = 'black')
          ,axis.text.x = element_text( face ="bold",size=10,colour = 'black')
          ,axis.ticks = element_line(colour = 'black',size=1))+
    ggtitle("Predicted vs. observed redshift z\n(cross-validated)")+
    xlab("Observed redshift z")+ylab("Predicted redshift z  (bars = min-max over CV folds)")

  print(plt1)
  dev.off()
  print(plt1)
}


save.image(file = paste(addr,"Workspace_outputs",plotnames,".Rdata",sep = ""))

return(results)
}
