rm(list=ls())
library(publicWW)
library(wwprev)
library(maptools)
library(spdep)
library(RColorBrewer)

###  setting working directory for saving results
workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev_final/results/full/'
setwd(workdir)

##   where to create the pdf of the map
output_file <- '/Volumes/WorkSpace/WBE/wwprev/doc/figs/wwprev_trends_ltla.pdf'


alldata_file <- '../alldata.rds'
alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
fitdata <- extract_from_alldata(3,c(1,43),0,1,alldata)
fitdata <- add_region(fitdata)

###   Kendall's Tau by LTLA trends
tau <- t(sapply(1:fitdata$nareas,function(x){unlist(cor.test(fitdata$w[x,],fitdata$p_mn[x,],method='kendall')[c('p.value','estimate')])}))
print(paste0(length(which(tau[,2]>0)),': number of LTLAs with tau >0 '))
print(paste0(length(which(tau[,2]>0 & tau[,1]<0.05)),': number of LTLAs with tau >0 with p.value < 5%'))
print(paste0(length(which(tau[,2]<0 & tau[,1]<0.05)),': number of LTLAs with tau <0 with p.value < 5%'))


###   load shapefile for mapping
shpfile <- system.file("extdata", "Local_Authority_Districts_(May_2021)_UK_BFE_V3/LAD_MAY_2021_UK_BFE_V2.shp", package = "wwprev")
ltla_shp <- readShapePoly(shpfile)
d <- ltla_shp@data
d$LAD21CD <- as.character(d$LAD21CD)
###  add City of London and Isle of Scilly back to the map
ltla_cds <- c(rownames(fitdata$w),'E06000053','E09000001')
ids <- sapply(ltla_cds,function(x){which(d$LAD21CD==x)})
poly <- (ltla_shp@polygons)[ids]
for (i in 1:length(poly)) poly[[i]]@ID <- ltla_cds[i]
pp <- SpatialPolygons(as.list(poly))

###   regional shapefile
shpfile_rgn <- system.file("extdata", "Regions_December_2022_EN_BFE_5024249782331847868/RGN_DEC_2022_EN_BFE.shp", package = "wwprev")
rgn_shp <- readShapePoly(shpfile_rgn)


###  create map
ncolors <- 6
cpts <- quantile(tau[,2],seq(0,1,length.out=(ncolors+1)))
ltla_l <- cut(tau[,2],breaks=cpts,labels=1:ncolors,include.lowest=TRUE)
shadings <- brewer.pal(ncolors-1,'BuGn')
ltla_shadings <- c(rgb(101,67,33,alpha=50,maxColorValue=255),shadings)
pdf(file=output_file,height=12,width=12)
par(mar=rep(0,4))
lcol <- c(ltla_shadings[ltla_l],rep('black',2))
plot(pp,col=lcol,border='transparent')
plot(rgn_shp,add=TRUE,col='transparent',lwd=0.5,border=1)
add_legend(16,48,40000,ltla_shadings,cpts,dt=18000,cex=2)
add_legend_main(24,52,'Correlation between \nWW and prevalence')
dev.off()



if (FALSE) {
	#  concordance between prevalence and wastewater measurements
	#  both are scaled to have mean 0 for ease of comparison
	ic <- 0
	par(mfrow=c(4,4))
	par(mar=rep(2,4))
	ids <- which(tau[,2]>0.4 & tau[,1]<0.05)
	for (i in ids) {
		matplot(cbind(scale(fitdata$w[i,]),scale(fitdata$p_mn[i,])),type='b',pch=19)
		if (ic==0) legend('bottomright',c('WW','prev'),col=c(1,2),pch=19,bty='n')
		ic <- ic + 1
	}
}
