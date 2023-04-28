###  source('/Volumes/WorkSpace/OnGitHub/wwprev/inst/scripts/summarising.R')


rm(list=ls())
library(wwprev)  #  run install_wwprev.R to install this package
library(ComplexHeatmap)
library(circlize)

workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/results/nimble_keep/forecast/'
setwd(workdir)

save_pdf <- TRUE  #  save figures to pdf?
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

###  overall summary
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

###  define color scheme for figures
diff <- pred[2,,] - obs
rownames(diff) <- NULL
rg <- round(range(c(diff)),digits=1)
m <- max(abs(rg))
#rg <- c(-m, 0, m)
rg <- c(rg[1],0,rg[2])
#rg <- c(-1.5,0,1.5)
col_fun <- colorRamp2(rg, c("brown", "lightgreen", "darkblue"))

###  disaggregated vs observed: scatterplot
l <- range(pred[2,,],obs)
#tiff(file='../../../figures/scatter.tif',height=500,width=500)
tiff(file='../../../figures/scatter.tif',res=200,width=6,height=6,unit='in')
plot(pred[2,,],obs,pch=19,cex=0.6
    ,xlab='Disaggregated logit prevalence'
    ,ylab='Estimated logit prevalence',xlim=l,ylim=l
    ,col=col_fun(diff))
points(pred[2,,],obs,lwd=0.02,cex=0.65)
abline(0,1,col=2,lwd=2)
dev.off()

###  heatmap of difference
tiff(file='../../../figures/heatmap.tif',res=200,width=6,height=6,unit='in')
Heatmap(diff,cluster_rows=FALSE,cluster_columns=FALSE
       ,col=col_fun
       ,heatmap_legend_param=list(title='',legend_height=unit(10,'cm')
                                 ,labels_gp = gpar(fontsize = 15))
       ,column_title='Difference in logit prevalence: Disaggregated - Observed'
       ,row_title='Lower Tier Local Authorities')
dev.off()

###  choropleth of LTLA RMSE averaged across weeks
md <- apply(diff,1,RMSE)
