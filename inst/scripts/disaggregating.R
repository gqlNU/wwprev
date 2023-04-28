###  source('/Volumes/WorkSpace/OnGitHub/wwprev/inst/scripts/disaggregating.R')


rm(list=ls())
library(wwprev)  #  run install_wwprev.R to install this package

workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/results/nimble_keep/'
setwd(workdir)
if (!dir.exists('params')) dir.create('params')
if (!dir.exists('forecast')) dir.create('forecast')

###  user inputs
setting <- list()
setting$imodel <- 3
setting$fit_ts <- c(1,20)
setting$horizon <- 5
setting$use_pmean <- FALSE
if (setting$imodel==1) setting$use_pmean <- TRUE
setting$recollect <- TRUE

alldata_file <- '../../alldata.rds'
alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)

###  get the list of files containing the fits
setting$draw_files <- get_nimble_files(setting)
setting$ndraws <- length(setting$draw_files)

###  get disaggregated prevalence from each fit
ds <- get_dis_prev(setting,alldata)

###  summarise disaggregated prevalence
my_print('summarising disaggregated prevalence')
forecast <- summarise_dis_prev(ds,setting)

###  saving disaggregated prevalence
my_print('saving disaggregated prevalence')
save(file=forecast$save_file,forecast)
my_print('done')
