#' Root mean square error
#'
#'
#'
#' @param
#' @return
#' @export
RMSE <- function(x) {
    o <- sqrt(mean(c(x)^2))
    return(o)
}

#' mean bias
#'
#'
#' @param
#' @return
#' @export
MB <- function(x) {
    o <- mean(c(x))
    return(o)
}

#' mean absolute bias
#'
#'
#'
#' @param
#' @return
#' @export
MAB <- function(x) {
    o <- mean(abs(c(x)))
    return(o)
}

#' coverage rate
#'
#'
#'
#' @param
#' @return
#' @export
coverage <- function(obs,pred_ci) {
	nareas <- nrow(obs)
	ntimes <- ncol(obs)
	cr <- array(0,c(nareas,ntimes))
	for (i in 1:nareas) {
	    for (j in 1:ntimes) {
			cr[i,j] <- point_in_interval(obs[i,j],pred_ci[,i,j])
    	}
	}
  	out <- mean(cr)
  	return(out)
}


#' summarise quality of disaggregated prevalence
#'
#'
#'
#' @param
#' @return
#' @export
forecast_summary <- function(pred_mean,pred_ci,obs) {
  #  summary of bias
  x <- pred_mean - obs # bias
  out <- c(MB(x),MAB(x),RMSE(x))
  #  coverage
  cover <- coverage(obs,pred_ci)
  out <- c(out,cover)
  names(out) <- c('MB','MAB','RMSE','coverage')
  return(out)
}

#' to check if a point (pt) is inside a given interval (itn)
#'
#' Is itn[1]<= pt <= itn[2] true?
#'
#' @param
#' @return
#' @export
point_in_interval <- function(pt,itn) {
    out <- 0
    if (pt<=itn[2] & pt>=itn[1]) out <- 1
    return(out)
}


#' read prediction
#'
#'
#'
#' @param
#' @return
#' @export
read_pred <- function(setting,alldata) {
	region <- setting$region
	reduced_coarse <- setting$reduced_coarse
	imodel <- setting$imodel
	fit_ts <- setting$fit_ts
	horizon <- setting$horizon
	use_pmean <- setting$use_pmean
	ndraws <- setting$ndraws
	if (is.null(region)) {
		rgn <- ''
	} else {
		rgn <- paste0('_',region)
	}
	
	###  overall summary
	file <- paste0('summary_nimble_multiscale_m',imodel,rgn,'_t',fit_ts[1],'-t',fit_ts[2],'h',horizon)
	if (use_pmean | imodel==1) {
	  file <- paste0(file,'_pmean.RData')
	} else {
	  file <- paste0(file,'_',ndraws,'MIs.RData')
	}
	load(file)
	
	t_ids <- (1:horizon) + fit_ts[2]
	if (!is.null(region)) {
		ltla <- dimnames(forecast$across_allMIs)[[2]]
		ids <- sapply(ltla,function(x){which(rownames(alldata$p_mn)==x)})
		obs <- alldata$p_mn[ids,t_ids]
	} else {
		obs <- alldata$p_mn[,t_ids]
	}
	pred <- forecast$across_allMIs[,,t_ids]
	out <- list(obs=obs,pred=pred)
	return(out)
}




#' preparing for nowcast summary (setting up working directory, file to save etc)
#'
#'
#'
#' @param
#' @return
#' @export
prepare_for_summary <- function(setting) {
	f_root <- setting$workdir_root
	###  get LTLA estimates over the disaggregation period
	alldata_file <- paste0(setting$workdir_root,'alldata.rds')
	alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)

	###  file path to save nowcast summary
	f <- paste0(f_root,'overall_summary/m',setting$imodel,'.rds')
	if (setting$reduced_coarse) {
		if (is.null(setting$region)) {
			f <- paste0(f_root,'overall_summary/m',setting$imodel,'_national_reduced.rds')
		} else {
			f <- paste0(f_root,'overall_summary/m',setting$imodel,'_pooled_regions_reduced.rds')
		}
	} else {
		if (is.null(setting$region)) {
			f <- paste0(f_root,'overall_summary/m',setting$imodel,'_national.rds')
		} else {
			f <- paste0(f_root,'overall_summary/m',setting$imodel,'_pooled_regions.rds')
		}
	}
	setting$save_file <- f
	
	setting$regions <- c('London','North East','East Midlands','East of England','South East',
	                     'South West','West Midlands','Yorkshire and The Humber','North West')
	setting$nregions <- length(setting$regions)
	
	workdir_root <- setting$workdir_root
	workdir <- paste0(workdir_root,'full/forecast/')
	if (setting$reduced_coarse | !is.null(setting$region)) {
		workdir <- paste0(workdir_root,'design/forecast/')
		if (setting$reduced_coarse) workdir <- paste0(workdir_root,'design/subset/forecast/')
	}
	setting$workdir <- workdir
	setwd(workdir)
	
	if (setting$imodel==1) setting$use_pmean <- TRUE
	
	out <- list(alldata=alldata,setting=setting)
	return(out)
}

#' summarising nowcast from a model
#'
#'
#'
#' @param
#' @return
#' @export
gather_nowcast_summary <- function(setting,alldata) {
	regions <- setting$regions
	nregions <- length(regions)
	s <- NULL
	if (is.null(setting$region)) {
		out <- read_pred(setting,alldata)
		all_obs <- out$obs
		pred <- out$pred
		###  by region
		alldata$w <- alldata$w_mn
		d <- add_region(alldata)
		for (ir in 1:nregions) {
			rgn <- regions[ir]
			ltla <- rownames(d$w_mn)[which(d$rgn_nm==rgn)]
			ids <- sapply(ltla,function(x){which(rownames(d$p_mn)==x)})
			#print(paste0('===', rgn))
			#rs <- forecast_summary(pred_mean=pred[2,ids,],pred_ci=pred[c(1,3),ids,],all_obs[ids,])
			#print(rs)
			s <- rbind(s,forecast_summary(pred_mean=pred[2,ids,],pred_ci=pred[c(1,3),ids,],all_obs[ids,]))
		}
		rownames(s) <- regions
	} else {
		if (setting$region=='separate') {
			nrs <- nregions
		} else {
			nrs <- 1
			regions <- setting$region
		}
		all_obs <- NULL
		pred_mn <- pred_ci1 <- pred_ci2 <- NULL
		for (ir in 1:nrs) {
			setting$region <- regions[ir]
			out <- read_pred(setting,alldata)
			obs <- out$obs
			pred <- out$pred
			all_obs <- rbind(all_obs,obs)
			pred_mn <- rbind(pred_mn,out$pred[2,,])
			pred_ci1 <- rbind(pred_ci1,out$pred[1,,])
			pred_ci2 <- rbind(pred_ci2,out$pred[3,,])
			#print(paste0('===', setting$region))
			rs <- forecast_summary(pred_mean=pred[2,,],pred_ci=pred[c(1,3),,],obs)
			#print(rs)
			s <- rbind(s,rs)
		}
		rownames(s) <- regions
		pred <- array(0,c(3,nrow(pred_mn),setting$horizon))
		pred[1,,] <- pred_ci1
		pred[2,,] <- pred_mn
		pred[3,,] <- pred_ci2
	}
	rs <- forecast_summary(pred_mean=pred[2,,],pred_ci=pred[c(1,3),,],all_obs)
	out <- list(region=s,overall=rs)
	return(out)
}
