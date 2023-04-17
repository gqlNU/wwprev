###  source('/Volumes/WorkSpace/OnGitHub/wwprev/inst/scripts/summarising.R')


rm(list=ls())
library(wwprev)  #  run install_wwprev.R to install this package

workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/results/nimble_keep/forecast/'
setwd(workdir)

imodels <- 2
fit_ts <- c(1,20)
horizon <- 20
ndraws <- 50
use_pmean <- FALSE

###  get LTLA estimates over the disaggregation period
alldata_file <- '../../../alldata.rds'
alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
t_ids <- (1:horizon) + fit_ts[2]
obs <- alldata$p_mn[,t_ids]


for (imodel in imodels) {
  file <- paste0('summary_nimble_multiscale_m',imodel,'_t',fit_ts[1],'-t',fit_ts[2],'h',horizon)
  if (use_pmean | imodel==1) {
    file <- paste0(file,'_pmean.RData')
  } else {
    file <- paste0(file,'_',ndraws,'MIs.RData')
  }
  load(file)
  pred <- forecast$across_allMIs[,,t_ids]
  print(file)
  print(forecast_summary(pred_mean=pred[2,,],pred_ci=pred[c(1,3),,],obs))
}
