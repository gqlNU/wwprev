###  source('/Volumes/WorkSpace/OnGitHub/wwprev/inst/scripts/summarising.R')


rm(list=ls())
library(wwprev)  #  run install_wwprev.R to install this package
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(publicWW)
library(sf)
library(cowplot)


workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/results/nimble_keep/forecast/'
setwd(workdir)

save_pdf <- TRUE  #  save figures to pdf?
imodel <- 2
fit_ts <- c(1,20)
horizon <- 15
ndraws <- 50
use_pmean <- FALSE

###  get LTLA estimates over the disaggregation period
alldata_file <- '../../../alldata.rds'
alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)
t_ids <- (1:horizon) + fit_ts[2]
obs <- alldata$p_mn[,t_ids]

###  overall summary
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

###  LAD2Region lookup
lad2rgn_lookup <- read.csv('/Volumes/WorkSpace/OnGitHub/wwprev/inst/extdata/Local_Authority_District_to_Region_(April_2021)_Lookup_in_England.csv',header=TRUE)

###  define color scheme for figures
diff <- pred[2,,] - obs
rownames(diff) <- NULL
rg <- round(range(c(diff)),digits=1)
m <- max(abs(rg))
#rg <- c(-m, 0, m)
rg <- c(rg[1],0,rg[2])
#rg <- c(-1.5,0,1.5)
col_fun <- colorRamp2(rg, c("brown", "lightgreen", "darkblue"))
col_fun <- colorRamp2(rg, c("brown", "white", "darkblue"))
col_fun <- colorRamp2(rg, c("brown", "beige", "darkblue"))
###  disaggregated vs observed: scatterplot
l <- range(pred[2,,],obs)
#tiff(file='../../../figures/scatter.tif',height=500,width=500)
f <- paste0('../../../figures/scatter_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.tif')
tiff(file=f,res=200,width=6,height=6,unit='in')
plot(obs,pred[2,,],pch=19,cex=0.6
    ,ylab='Disaggregated logit prevalence'
    ,xlab='Estimated logit prevalence',xlim=l,ylim=l
    ,col=col_fun(diff))
points(obs,pred[2,,],lwd=0.02,cex=0.65)
abline(0,1,col=2,lwd=2)
dev.off()

rgn21nm <- unique(lad2rgn_lookup$RGN21NM)
nregions <- length(rgn21nm)
rgn_rmse <- NULL
for (i in 1:length(rgn21nm)) {
	l <- lad2rgn_lookup$LAD21CD[which(lad2rgn_lookup$RGN21NM==rgn21nm[i])]
	xx <- unlist(sapply(l,function(x){which(rownames(obs)==x)}))
	d <- c(diff[xx,])
	rgn_rmse <- c(rgn_rmse,sqrt(mean(d^2)))
}

orgn21nm <- rgn21nm[order(rgn_rmse)]
rgn_ids <- rgn <- NULL
for (i in 1:nregions) {
	rg <- orgn21nm[i]
	l <- lad2rgn_lookup$LAD21CD[which(lad2rgn_lookup$RGN21NM==rg)]
	xx <- unlist(sapply(l,function(x){which(rownames(obs)==x)}))
	rgn_ids <- c(rgn_ids,xx)
	rgn <- c(rgn,rep(i,length(xx)))
}
txts <- as.list(orgn21nm)
ha <- rowAnnotation(foo=anno_empty(border=FALSE,width=max_text_width(unlist(txts)) + unit(0,'mm')))

###  heatmap of difference
f <- paste0('../../../figures/heatmap_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.tif')
tiff(file=f,res=200,width=8,height=6,unit='in')
rgn_diff <- diff[rgn_ids,]
Heatmap(rgn_diff,cluster_rows=FALSE,cluster_columns=FALSE
       ,col=col_fun
       ,heatmap_legend_param=list(title='',legend_height=unit(10,'cm')
                                 ,labels_gp = gpar(fontsize = 10))
#       ,column_title='Difference in logit prevalence: Disaggregated - Observed'
       ,column_title=''
       ,row_title='Lower Tier Local Authorities'
       ,row_split=rgn
       ,right_annotation=ha)

for(i in 1:nregions) {
    decorate_annotation("foo", slice = i, {
        grid.rect(x = 0, width = unit(0.7, "mm"), gp = gpar(fill = 'darkgrey', col = NA), just = "left")
        grid.text(paste(txts[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left")
    })
}
dev.off()

###  choropleth of LTLA RMSE averaged across weeks
md <- apply(diff,1,RMSE)
names(md) <- rownames(obs)


ltla_map <- function(ltla.shp,values,inset.region='London',colour.scheme=NULL) {
  shp <- merge(ltla.shp,values,by='LAD21CD',all.x=FALSE,all.y=TRUE)
  #  main map
  if (is.null(colour.scheme)) {
    main.map <- shp %>%
    ggplot() +
      geom_sf(aes(fill = mean),lwd = 0.001,colour = "light grey") +
#      scale_fill_gradient(low = "grey", high = "brown")
	  scale_fill_distiller(palette = 1,name='RMSE',direction=1) + 
	  theme_void() +
	  theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18),
        legend.key.size = unit(.7, 'cm')
      ) + 
      coord_sf(expand = FALSE)  # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
  } else {
    main.map <- shp %>%
      ggplot() +
	  geom_sf(aes(fill = cut(mean,breaks= colour.scheme$breaks,labels= colour.scheme$labels)),
	                 lwd = 0.001,colour = "light grey") +
	  scale_fill_manual(values=colour.scheme$cols,limits=colour.scheme$labels) +
	  theme_void() +
	  theme(legend.position='none') +
	  coord_sf(expand = FALSE)  # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
  }
  #  putting the main map and the inset together
  if (!is.null(inset.region)) {
    coords <- get_inset_coords(inset.region,ltla.shp)
    zx <- coords$zx
    zy <- coords$zy
    #  add a rectangle on the main map to indicate inset
    main.map <- main.map + geom_rect(xmin = zx[1],ymin = zy[1],xmax = zx[2],ymax = zy[2],
                                                                fill = NA,colour = "black",size = 0.6)
    final <- ggdraw(main.map) +
                draw_plot({main.map + coord_sf(xlim = zx,ylim = zy,expand = FALSE) +
                                   theme(legend.position = "none")
                                  },
                                  # The distance along a (0,1) x-axis to draw the left edge of the plot
                                  x = 0.02,
                                  # The distance along a (0,1) y-axis to draw the bottom edge of the plot
                                  y = 0.2,
                                  # The width and height of the plot expressed as proportion of the entire ggdraw object
                                  width = 0.35,height = 0.35)
  } else {#  no inset
  	final <- map.map
  }
  return(final)
}


###  choropleth of LTLA-RMSE
shpfile <- '/Volumes/WorkSpace/publicWW_dev/inst/extdata/Local_Authority_Districts_(May_2021)_UK_BFE_V3/LAD_MAY_2021_UK_BFE_V2.shp'
ltla_shp <- sf::st_read(shpfile)

values <- data.frame(mean=apply(diff,1,RMSE),LAD21CD=rownames(obs))
p2 <- ltla_map(ltla_shp,values,inset.region='London')
f <- paste0('../../../figures/RMSE_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.tif')
tiff(file=f,res=200,width=6,height=6,unit='in')
print(p2)
dev.off()






###  observed vs disaggregated time series plots
obs_pred_ts_plot <- function(id,obs,disaggregated,fit_ts,ltla_names=NULL,rgn_names=NULL,ylim=NULL) {
	a <- disaggregated[,id,]
	b <- obs[id,]
	tt <- length(b)
	nm <- rownames(obs)[id]
	if (!is.null(ltla_names)) nm <- paste0(nm,' ',ltla_names[id])
	t_ids <- (fit_ts[2]+1):length(b)
	pm <- matrix(a[2,t_ids],nrow=1)
	pci <- array(0,c(2,1,length(t_ids)))
	pci[,1,] <- a[c(1,3),t_ids]
	o <- matrix(b[t_ids],nrow=1)
	s <- forecast_summary(pm,pci,o)
	nm <- paste0(nm,'\n RMSE=',round(s['RMSE'],digits=3),'   95%cov=',round(s['coverage']*100,digits=2),'%')
	if (is.null(ylim)) ylim <- range(c(a,b))
	matplot(t(a),type='l',ylim=ylim,col='transparent'
	       ,xlab='week ID',ylab='logit prevalence',main=nm)
	polygon(c(1:tt,tt:1),c(a[1,],a[3,tt:1]),col='lightgrey',border='lightgrey')
	lines(a[2,])
	points(b,pch=19)
	abline(v=fit_ts[2]+0.4,lty=2)
	legend('bottomright',legend=rgn_names[id],bty='n',lty=-1)
}



obs <- alldata$p_mn[,1:(fit_ts[2]+horizon)]
pred <- forecast$across_allMIs
rg <- range(c(obs,pred))
ltla_names <- ltla_shp$LAD21NM[sapply(rownames(obs),function(x) {which(ltla_shp$LAD21CD==x)})]
rgn_names <- lad2rgn_lookup$RGN21NM[sapply(rownames(obs),function(x) {which(lad2rgn_lookup$LAD21CD==x)})]


f <- paste0('../../../figures/ts_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.pdf')
pdf(file=f,width=5*3,height=3*3)
#ids <- order(values$mean,decreasing=TRUE)[1:9]
#ids <- which(rgn_names=='North East')[1:9]
for (ir in 1:nregions) {
	ids <- which(rgn_names==orgn21nm[ir])
	par(mfrow=c(3,5))
	for (id in ids) obs_pred_ts_plot(id,obs,pred,fit_ts,ltla_names,rgn_names,ylim=rg)
}
dev.off()


f <- paste0('../../../figures/rmse_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.RData')
values$rgn_names <- rgn_names
save(file=f,values)

if (FALSE) {
###  regional variation in RMSE
boxplot(values$mean~rgn_names,cex.axis=0.6)
summary(lm(values$mean~rgn_names))

f <- paste0('../../../figures/rmse_m',imodel,'_t',fit_ts[1],'to',fit_ts[2],'_h',horizon,'.RData')
values$rgn_names <- rgn_names
save(file=f,values)


l <- pred[c(1,3),,(fit_ts[2]+1):ncol(obs)]
dl <- l[2,,]-l[1,,]
summary(c(dl))
}

f1 <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/figures/rmse_m1_t1to20_h20.RData'
f2 <- '/Volumes/WorkSpace/OnGitHub/wwprev/forTesting/figures/rmse_m2_t1to20_h20.RData'
if (file.exists(f1) & file.exists(f2)) {
	library(ggplot2)
	load(f1)
	m1 <- values
	m1$model <- 1
	load(f2)
	m2 <- values
	m2$model <- 2
	m <- rbind(m1,m2)
	m$r <- m$mean
	m$model <- as.character(m$model)
	m$rgn_names[which(m$rgn_names=='Yorkshire and The Humber')] <- 'Yorkshire and \nThe Humber'
	p <- m %>%
  ggplot(aes(x=rgn_names,y=r, color=model)) +
  labs(y= "RMSE", x = "Regions") +
#  geom_boxplot(width=.5) + 
  geom_boxplot(outlier.shape=NA)#+
#  geom_point(position=position_jitterdodge())

#  geom_jitter(width=0.15)

f <- paste0('../../../figures/compareRMSE.pdf')
pdf(file=f,width=10,height=3)
print(p)
dev.off()
}