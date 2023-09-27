rm(list=ls())

library(wwprev) 
library(publicWW)

####################
#  user input
####################
setting <- list()
setting$imodel <- 2
setting$fit_ts <- c(1,20)
setting$horizon <- 20
setting$use_pmean <- !TRUE
setting$ndraws <- 50

setting$region <- 'separate'  #  separate = nowcast is done for each region separately
#setting$region <- NULL        #  NULL = nowcast is done for the whole nation
#setting$region <- 'London'    #  e.g. London = only summarising nowcast for that region with the model only fitted to that region

setting$reduced_coarse <- !FALSE  # TRUE  = national prevalence only available for a subset of weeks during nowcast
                                  # FALSE = national prevalence available for ALL weeks during nowcast
setting$workdir_root <- '/Volumes/WorkSpace/OnGitHub/wwprev_final/results/'
####################
#  end user input
####################

#  some preparation steps for summary
out <- prepare_for_summary(setting)
setting <- out$setting
alldata <- out$alldata

#  summarising nowcast
(out <- gather_nowcast_summary(setting,alldata))

#  saving nowcast regional + overall summary
saveRDS(out,file=setting$save_file)

