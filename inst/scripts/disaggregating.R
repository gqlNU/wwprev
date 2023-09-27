###  source('/Volumes/WorkSpace/OnGitHub/wwprev_final/inst/scripts/disaggregating.R')

rm(list=ls())
library(wwprev)

########################
###  user inputs
########################
setting <- list()
setting$imodel <- 2
setting$fit_ts <- c(1,20)
setting$horizon <- 20
setting$use_pmean <- FALSE
setting$reduced_coarse <- !FALSE  # TRUE  = national prevalence only available for a subset of weeks during nowcast
                                  # FALSE = national prevalence available for ALL weeks during nowcast
setting$which_region <- 0  #  0=all regions and 1-9 picks one of the
                           #  9 English regions (see the object regions in wwprev::prepare_to_disaggregate for definition)
                           
setting$workdir_root <- '/Volumes/WorkSpace/OnGitHub/wwprev_final/results/'

########################
###  end user inputs
########################

###  preparation for disaggregation
out <- prepare_to_disaggregate(setting)
alldata <- out$alldata
setting <- out$setting
setwd(setting$workdir)

###  get disaggregated prevalence from each fit
ds <- get_dis_prev(setting,alldata)

###  summarise disaggregated prevalence
my_print('summarising disaggregated prevalence')
forecast <- summarise_dis_prev(ds,setting)

###  saving disaggregated prevalence
my_print('saving disaggregated prevalence')
save(file=forecast$save_file,forecast)
my_print('done')
