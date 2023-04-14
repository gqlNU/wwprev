###  source('/Volumes/WorkSpace/OnGitHub/wwprev/inst/scripts/fitting.R')


rm(list=ls())
library(devtools)
library(nimble)
library(publicWW)
library(dplyr)
library(wwprev)  #  run install_wwprev.R to install this package

###  setting working directory for saving results
workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/results/nimble_keep/'
setwd(workdir)

alldata_file <- '../../alldata.rds'
if (!file.exists(alldata_file)) {
    alldata <- get_all_ww_prev_data(save_file=alldata_file,load_file=NULL)
} else {
    alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
}

imodel <- as.numeric(readline(prompt='imodel =  '))
horizon <- as.numeric(readline(prompt='horizon = '))
start_draw <- as.numeric(readline(prompt='start draw = '))
end_draw <- as.numeric(readline(prompt='end draw = '))


fit_ts <- c(1,20)
iters <- c(50000,10000,80)
draws_ww <- c(floor(seq(1,1000,length.out=50)))[start_draw:end_draw]
if (imodel==1) draws_ww <- 0
ndraws_ww <- length(draws_ww)

idw <- 0
for (idraw in draws_ww) {
    idw <- idw + 1
    fitdata <- extract_from_alldata(imodel,fit_ts,horizon,idraw,alldata)
    output <- fit_nimble(imodel,iters,fitdata)
    save_nimble_fit(imodel,output,workdir,fitdata,iters)
    msg <- paste0('done ',idw,' out of ',ndraws_ww)
    my_print(msg)
}
