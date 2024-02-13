#' used for producing the wastewater vs wastewater+national-prevalence nowcast
#'
#'
#' @return 
#' @export
compare_forecasts <- function(fs,obs,setting) {
	tds <- 1:nrow(fs[[1]]$f)
	nfs <- length(fs)
	yl <- NULL
	for (i in 1:nfs) yl <- range(c(yl,fs[[i]]$f))
	yl <- c(-7.5,-2)
	m <- paste0(setting$this_ltla_cd,'\n(',setting$main,')')
	matplot(fs[[1]]$f,type='l',ylim=yl,col='transparent',
	        main=m,xaxt='n',ylab=setting$ylab,cex.lab=1.2)
	axis(1,at=c(1,15,28,41),labels=c('Jun','Sep','Dec','Mar'))
	abline(v=setting$fit_ts[2]+0.3,col='grey',lwd=1.3)
	##   going through the different forecasts
	cols <- NULL
	for (i in 1:nfs) {
		ff <- fs[[i]]
		add_one_forecast(ff$f,cols=ff$cols)
		cols <- c(cols,rgb(ff$cols[1],ff$cols[2],ff$cols[3],1))
	}
	points(obs[tds],pch=19,cex=1.4)
	##   legends
	if (!is.null(setting$legtext)) {
		legend('topleft',legend=setting$legtext,bty='n',col=cols,lwd=2)
		text(1,-7.3,'Local prevalence \navailable',cex=1.1,pos=4)
		text(23,-7.3,'Local prevalence \nunavailable',cex=1.1,pos=4)
	}
}


#' used for producing the wastewater vs wastewater+national-prevalence nowcast
#'
#'
#' @return 
#' @export
add_one_forecast <- function(f,cols) {
	tds <- 1:nrow(f)
	##   pmean line
	lines(tds,f[,2],col=rgb(cols[1],cols[2],cols[3],1),lwd=2)
	##   uncertainty
	rtds <- rev(tds)
	polygon(c(tds,rtds),c(f[,1],f[rtds,3]),
	        border=NA,col=rgb(cols[1],cols[2],cols[3],0.2))
}


#' Add a legend to the correlation plot
#'
#'
#' @return 
#' @export
add_legend <- function(x,y,d,cols,txt,dt,cex=2) {
	txt <- format(round(txt,digits=2))
	txt <- sapply(txt,function(tt){sub(' ','',tt)})
	x <- x*10000
	y <- y*10000
	ncs <- length(cols)
	ic <- 0
	xx <- c(x,x+d)
	for (i in ncs:1) {
		ic <- ic + 1
		yy <- c(y-d*ic,y-(ic-1)*d)
		rect(xx[1],yy[1],xx[2],yy[2],col=cols[i],border=cols[i])
		#  add text for each colour
		if (i!=1) tt <- paste0('(',txt[i],', ',txt[i+1],']')
		if (i==1) tt <- paste0('[',txt[i],', ',txt[i+1],']')
		text(x-dt,mean(yy),tt,pos=2,cex=cex)
	}
	ic <- ic + 1
	yy <- c(y-d*ic,y-(ic-1)*d)
	rect(xx[1],yy[1],xx[2],yy[2],col=1,border=1)
	#  add text for each colour
	text(x-dt,mean(yy),'Not included',pos=2,cex=cex)
	
}

#' Add a title to the legend of the correlation plot
#'
#'
#' @return 
#' @export
add_legend_main <- function(x,y,main, cex=2.3) {
	text(x*10000,y*10000,main,pos=2,cex=2.3)
}
