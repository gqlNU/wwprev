rm(list=ls())

library(devtools)
library(nimble)
library(publicWW)
library(dplyr)
library(car)
library(wwprev)

####################
#  user input
####################
#   directory to save the results in
wd <- '/Volumes/WorkSpace/OnGitHub/wwprev_final/results/'

imodel <- 1       #  which model to fit
fit_ts <- c(1,20) #  define the period in which fine-scale prevalence estimates are available
                  #  e.g. fit_ts <- c(1,20) means fine scale prevalence available for the first 20 weeks
horizon <- 20     #  define the nowcast time interval
                  #  e.g. horizon <- 20 means nowcast for the subsequent 20 weeks
idraw <- 1        #  which set of wastewater values sampled from the predictive distribution to use
iters <- c(50000,10000,80)  #  setting for the MCMC run, e.g. iters <- c(50000,10000,80)
                            #  means a total of 50k iterations with first 10k discarded with a thin of 80
                            #  (saving every 80th iteration for posterior summary)
                            #  By default, two MCMC chains are run in Stan

coarse_avail_times <- NULL  #  By default, national prevalence are available for every week over the nowcast period
coarse_avail_times <- c(25,30,35,40)   #  national prevalence are only available for some selected weeks
                                        #  e.g. coarse_avail_times <- c(25,30,35,40) available every 5 weeks

####################
#  end user input
####################

reduced_coarse <- FALSE
if (!is.null(coarse_avail_times)) reduced_coarse <- TRUE

###  setting working directory for saving results
workdir <- paste0(wd,'full/')
if (reduced_coarse) workdir <- paste0(wd,'design/')
if (!dir.exists(workdir)) {
	dir.create(workdir,recursive=TRUE)
}
setwd(workdir)

###  create the master dataset
alldata_file <- '../alldata.rds'
if (!file.exists(alldata_file)) {
    alldata <- get_all_ww_prev_data(save_file=alldata_file,load_file=NULL)
} else {
    alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
}

###  format the master data for model fitting
fitdata <- extract_from_alldata(imodel,fit_ts,horizon,idraw,alldata)

###  format the fitdata if national prevalence is only available 
###  for some selected weeks during the nowcast period
if (reduced_coarse) {
	fitdata <- fitdata_with_reduced_coarse(fitdata,fit_ts,coarse_avail_times)
}

###  model fitting
output <- fit_nimble(imodel,iters,fitdata,reduced_coarse)

###  save model fit
save_nimble_fit(imodel,output,workdir,fitdata,iters)

