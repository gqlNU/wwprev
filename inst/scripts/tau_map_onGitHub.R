library(prettymapr)  #  scale bar and north arrow
library(wwprev)
library(sf)
library(RColorBrewer)

f <- system.file("extdata", "tau.rds", package = "wwprev")
tau <- readRDS(f)
###  add City of London and Isle of Scilly back to the map
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

ncolors <- 6
cpts <- quantile(tau[,2],seq(0,1,length.out=(ncolors+1)))
ltla_l <- cut(tau[,2],breaks=cpts,labels=1:ncolors,include.lowest=TRUE)
shadings <- brewer.pal(ncolors-1,'BuGn')
ltla_shadings <- c(rgb(101,67,33,alpha=50,maxColorValue=255),shadings)
lcol <- c(ltla_shadings[ltla_l],rep('black',2))


pdf(file='fig2.pdf',height=12,width=12)
par(mar=rep(0,4))
plot(st_geometry(ltla),col=lcol,border='transparent')
plot(st_geometry(rgn_shp),add=TRUE,col='transparent',lwd=0.5,border=1)
add_legend(16,48,40000,ltla_shadings,cpts,dt=18000,cex=2)
add_legend_main(26,52,'Correlation between \nWW and prevalence')
addnortharrow(padin=c(1.5,1),scale=1.6)
addscalebar(plotunit='m',pos='topright',padin=c(1,0.5),htin=0.1, widthhint=0.2,label.cex=2)
dev.off()

