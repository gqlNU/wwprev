rm(list=ls())
library(wwprev)

workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev_final/results/full/forecast/'
setwd(workdir)

#  where to save the pdf to
output_file <- '/Volumes/WorkSpace/WBE/wwprev/doc/figs/ww_only_data_integration.pdf'

fit_ts <- c(1,20)
horizon <- 20
setting <- list(fit_ts=fit_ts,horizon=horizon)

###  observed
alldata_file <- '../../alldata.rds'
alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
allobs <- alldata$p_mn
fitdata <- extract_from_alldata(3,c(1,40),0,1,alldata)
fitdata <- add_region(fitdata)

###  forecasting
f <- system.file("extdata", "summary_nimble_multiscale_m2_t1-t20h20_50MIs.RData", package = "wwprev")
load(f)  #  wastewater + national prevalence
ft1 <- forecast
f <- system.file("extdata", "summary_nimble_multiscale_m7_t1-t20h20_50MIs.RData", package = "wwprev")
load(f)  #  wastewater only
ft2 <- forecast



blue <- c(0,71/255,171/255)
red <- c(150/255,0,0)
LTLA_cds <- rownames(alldata$p_mn)
setting$legtext <- c('Data integration','With WW only')
rgns <- unique(fitdata$rgn_nm)
set.seed(1)
chosen_areas <- sapply(rgns,function(x) {
	#ii <- which(fitdata$rgn_nm==x)[1]
	ii <- which(fitdata$rgn_nm==x)
	out <- sample(ii,1)
	return(out)
})

pdf(file=output_file,height=9,width=9)
par(mfrow=c(3,3))
par(mar=c(2,4.2,4,1))
ic <- 0
for (iarea in chosen_areas) {
	ic <- ic + 1
	f1 <- list(f=t(ft1$across_allMIs[,iarea,]),cols=red)
	f2 <- list(f=t(ft2$across_allMIs[,iarea,]),cols=blue)
	allfs <- list(f1,f2)
	setting$main <- names(chosen_areas)[ic]
	setting$this_ltla_cd <- LTLA_cds[iarea]
	if (ic!=1) setting$legtext <- NULL
	setting$ylab <- ''
	if (any(ic==c(1,4,7))) setting$ylab <- 'logit COVID prevalence'
	compare_forecasts(allfs,allobs[iarea,],setting)
}
dev.off()