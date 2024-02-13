#  rm(list=ls())
library(prettymapr)  #  scale bar and north arrow
library(wwprev)
library(sf)
library(RColorBrewer)

###  load parameter summary
f <- system.file("extdata", "all_para_summary.rds", package = "wwprev")
all_summary <- readRDS(f)

pdf(file='supp_fig3.pdf',height=18,width=24)
l <- c(rep(1,10),
       rep(1,10),
       rep(1,10),
       rep(2:3,each=5),
       rep(2:3,each=5),
       rep(2:3,each=5),
       rep(2:3,each=5),
       rep(2:3,each=5),
       rep(2:3,each=5),
       rep(2:3,each=5))
layout(matrix(l,nrow=10,byrow=TRUE))

###  national trend of prevalence
b <- all_summary$b
nweeks <- nrow(b)
matplot(b,col='transparent',ylab='',main='(A) Posterior estimate of the national trend of prevalence on the logit scale',xlab='',xaxt='n',cex.axis=4,cex.main=4)
axis(1, at = c(1, 15, 28, 41), labels = c("Jun", "Sep", "Dec", "Mar"),cex.axis=3)
abline(v=20.3,lwd=4,lty=2)
text(20-10,-6,'training period',cex=4)
text(25,-6,'nowcast period',cex=4)
polygon(c(1:nweeks,nweeks:1),c(b[,1],b[nweeks:1,3]),border='lightgrey',col='lightgrey')
points(b[,2],pch=19,cex=3)
lines(b[,2],lwd=3)

###   load shapefile for mapping
f <- system.file("extdata", "tau.rds", package = "wwprev")
tau <- readRDS(f)
ltla_cds <- c(rownames(tau),'E06000053','E09000001')

###   load shapefile for mapping
shpfile <- system.file("extdata", "Local_Authority_Districts_(May_2021)_UK_BFE_V3/LAD_MAY_2021_UK_BFE_V2.shp", package = "wwprev")
ltla_shp <- st_read(shpfile)
d <- as.character(ltla_shp$LAD21CD)
ids <- sapply(ltla_cds,function(x){which(d==x)})
ltla <- ltla_shp[ids,]

###   regional shapefile
shpfile_rgn <- system.file("extdata", "Regions_December_2022_EN_BFE_5024249782331847868/RGN_DEC_2022_EN_BFE.shp", package = "wwprev")
rgn_shp <- st_read(shpfile_rgn)

ncolors <- 9

###  map of U[i]
par(mar=rep(2,4))
m <- all_summary$U
shadings <- brewer.pal(ncolors,'RdBu')
ltla_shadings <- shadings[length(shadings):1]
cpts <- quantile(m,seq(0,1,length.out=(ncolors+1)))
ltla_l <- cut(m,breaks=cpts,labels=1:ncolors,include.lowest=TRUE)

lcol <- c(ltla_shadings[ltla_l],rep('black',2))
plot(st_geometry(ltla),col=lcol,border='lightgrey',lwd=0.1)
add_legend(18,55,40000,ltla_shadings,cpts,dt=18000,cex=3)
add_legend_main(32,60,'(B) posterior mean of Ui',cex=6)

###  map of d+M[i]
par(mar=rep(2,4))
m <- all_summary$dM
shadings <- brewer.pal(ncolors,'BuGn')
ltla_shadings <- shadings
cpts <- quantile(m,seq(0,1,length.out=(ncolors+1)))
ltla_l <- cut(m,breaks=cpts,labels=1:ncolors,include.lowest=TRUE)

lcol <- c(ltla_shadings[ltla_l],rep('black',2))
plot(st_geometry(ltla),col=lcol,border='lightgrey',lwd=0.1)
add_legend(18,55,40000,ltla_shadings,cpts,dt=18000,cex=3)
add_legend_main(35,60,'(C) posterior mean of d+Mi',cex=6)
addnortharrow(padin=c(1.5,1),scale=2)
addscalebar(plotunit='m',pos='topright',padin=c(1,0.5),htin=0.2, widthhint=0.28,label.cex=3)
dev.off()
