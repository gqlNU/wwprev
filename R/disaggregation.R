#' disaggregating national prevalence to LTLA
#'
#'
#'
#' @param
#' @return
#' @export
disaggregate <- function(sims_list_file,imodel,alldata,region=NULL) {
    load(sims_list_file)
    nareas <- fitdata$nareas
    ntimes <- fitdata$ntimes
    ntot_times <- fitdata$ntot_times
    d <- dim(sims.list$mu_pv)
    nsims <- d[1]
    nchains <- d[4]
    d[3] <- ntot_times
    all <- array(0,d)
    idraw <- fitdata$idraw
    ltla <- 1:alldata$nareas
    if (fitdata$nareas < alldata$nareas) {
    	ltla <- sapply(rownames(fitdata$w),function(x){which(dimnames(alldata$w_sims)[[2]]==x)})
    }
    w <- alldata$w_sims[idraw,ltla,] - 8
    eps <- sims_noise(sims.list$sd_pv)
    for (i in 1:nareas) {
        for (t in (ntimes+1):ntot_times) {
        	if (imodel==1) {
                f <- sims.list$u[,i,] + sims.list$b[,t,]
            }
        	if (imodel==2 | imodel==4) {
                f <- sims.list$u[,i,] + sims.list$b[,t,] + (sims.list$gamma[,i,]+sims.list$mu_a)*w[i,t]
            }
            if (imodel==3) {
                f <- sims.list$u[,i,] + sims.list$b[,t,] + (sims.list$gamma[,i,]+sims.list$mu_a+sims.list$a[,t,])*w[i,t]
            }
            if (imodel==5) {
                f <- sims.list$u[,i,] + sims.list$b[,t,] + (sims.list$gamma[,i,]+sims.list$mu_a)*w[i,t] + 
                     sims.list$beta.imd * fitdata$std_imd[i]
            }
            if (imodel==6) {
                f <- sims.list$u[,i,] + sims.list$b[,t,] + (sims.list$gamma[,i,]+sims.list$mu_a)*w[i,t] + 
                     sims.list$beta.imd * fitdata$std_imd[i] + 
                     sims.list$beta.bame * fitdata$std_bame[i]
            }
            if (imodel==7) {
                f <- sims.list$u[,i,] + (sims.list$gamma[,i,]+sims.list$mu_a)*w[i,t]
            }
            if (imodel==8) {
                f <- sims.list$u[,i,] + sims.list$b[,t,] + (sims.list$gamma[,i,]+sims.list$mu_a)*w[i,t] + 
                     sims.list$beta.bame * fitdata$std_bame[i]
            }
            all[,i,t,] <- f + eps
        }
    }
    all[,,1:ntimes,] <- sims.list$mu_pv
    out <- list(all=all,nsims=nsims,nchains=nchains)
    if (fitdata$nareas < alldata$nareas) {
    	out$ltla <- rownames(fitdata$w)
    }
    return(out)
}

#' simulate a set of values from a normal distribution with mean 0 and given variance
#'
#'
#'
#' @param
#' @return
#' @export
sims_noise <- function(s) {
	if (is.null(dim(s))) {
		out <- sapply(s,function(x){rnorm(1,0,x)})
	} else {
		out <- apply(s,c(1,2),function(x){rnorm(1,0,x)})
	}
	return(out)
}

#' extract MCMC draws from a Nimble fit
#'
#'
#'
#' @param
#' @return
#' @export
get_nimble_MCMC_draws <- function(imodel,draw_file,return_type,save_sims_list=TRUE) {
    load(draw_file)
    nareas <- fitdata$nareas
    ntimes <- fitdata$ntimes
    ntot_times <- fitdata$ntot_times
    nchains <- length(output)
    nsims <- nrow(output[[1]])
    ###  a list of parameters to be extracted from the Nimble results
    all_params <- list(alpha=0,sd_pv=0,sd_P=0
                      ,mu_pv=c(nareas,ntimes)
                      ,sd_b=0,b=ntot_times
                      ,sd_u=0,u=nareas)
    addp <- NULL
    if (imodel>1) {
        addp <- list(mu_a=0,sd_gamma=0,gamma=nareas)
        if (imodel==3) {
            addp <- c(addp,list(sd_a=0,a=ntot_times))
        }
        if (imodel==4) {
	        addp <- c(addp,list(sd_s=0,s=nareas))
        }
        if (imodel==5 | imodel==6) {
	        addp <- c(addp,list(beta.imd=0))
	        if (imodel==6) {
	        		addp <- c(addp,list(beta.bame=0))
	        }
        }
        if (imodel==8) {
        	addp <- c(addp,list(beta.bame=0))
        }
    }
    all_params <- c(all_params,addp)
    my_print('extracting draws')
    sims.list <- extract_all_parameters(params=all_params,output=output,start=1)
    save_file <- paste0('params/',draw_file)
    if (save_sims_list) {
        save(file=save_file,sims.list,fitdata)
        my_print('saving sims.list')
    }
    if (return_type=='MCMC_draws') {
        return(sims.list)
    }
    if (return_type=='file_path') {
        return(save_file)
    }
}


#' get the file names for all fits associated with a given model and horizon
#'
#'
#'
#' @param
#' @return
#' @export
get_nimble_files <- function(setting) {
    imodel <- setting$imodel
    fit_ts <- setting$fit_ts
    horizon <- setting$horizon
    use_pmean <- setting$use_pmean
    recollect <- setting$recollect
    region <- rgn <- setting$region
    if (region=='') {
#    if (is.null(region)) {
    	rgn <- ''
    } else {
    	rgn <- paste0('_',rgn)
    }
    p1 <- paste0('nimble_multiscale_m',imodel)
    p2 <- paste0(rgn,'_t',fit_ts[1],'-t',fit_ts[2],'h',horizon,'idraw')
    if (use_pmean) {
        p2 <- paste0(rgn,'_t',fit_ts[1],'-t',fit_ts[2],'h',horizon,'_pmean')
        recollect <- TRUE
    }
    p3 <- paste0(p1,p2)
    f <- list.files()
#    f1 <- f[grep(p1,f)]
#    out <- f1[grep(p2,f1)]
    out <- f[grep(p3,f)]
    if (!recollect) {
    ###  see how many of the sims.list of those draws have been collected
        cl <- list.files('params')
        cl1 <- cl[grep(p1,cl)]
        collected <- cl1[grep(p2,cl1)]
        ids <- sapply(collected,function(x){which(out==x)})
        out <- out[-c(ids)]
    }
    return(out)
}




#' extract a given parameter from nimble fit
#'
#'
#'
#' @param
#' @return
#' @export
extract_params <- function(param,dms,output,start,end=NULL,nms) {
    pm <- colnames(output[[1]])
    if (is.null(end)) end <- nrow(output[[1]])
    iters <- start:end
    niters <- length(iters)
    nchains <- length(output)
    if (length(dms)==1) {
        if (dms==0) {
            ###  a scale parameter
            pm <- param
            id <- which(nms==pm)
            out <- cbind(output[[1]][iters,id],output[[2]][iters,id])
        } else {
            ###  a vector of parameters
            out <- array(0,c(niters,dms,nchains))
            for (i in 1:dms) {
                pm <- paste0(param,'[',i,']')
                id <- which(nms==pm)
                out[,i,] <- cbind(output[[1]][iters,id],output[[2]][iters,id])
            }
        }
    }
    if (length(dms)==2) {
        ###  a matrix of parameters
        out <- array(0,c(niters,dms[1],dms[2],nchains))
        for (i in 1:dms[1]) {
            for (j in 1:dms[2]) {
                pm <- paste0(param,'[',i,',',j,']')
                id <- which(nms==pm)
                out[,i,j,] <- cbind(output[[1]][iters,id],output[[2]][iters,id])
            }
        }
    }
    if (length(dms)>2) stop('Error in extract_params: cannot deal with a parameter array with 3 or more dimensions yet')
    return(out)
}


#' a wrapper function of extract_params to extract a collection of parameters
#'
#'
#'
#' @param
#' @return
#' @export
extract_all_parameters <- function(params,output,start,end=NULL,ndraws=NULL) {
#	if (is.null(ndraws)) output <- list(output=ouput)
    nms <- sub(' ','',colnames(output[[1]]))
    dms <- params
    nparams <- length(params)
    out <- params
    for (ip in 1:nparams) {
        if (is.null(ndraws)) {
            out[[ip]] <- extract_params(names(params)[ip],dms[[ip]],output,start=start,end=end,nms=nms)
        } else {
            tmp <- as.list(1:ndraws)
            for (idraw in 1:ndraws) {
                tmp[[idraw]] <- extract_params(names(params)[ip],dms[[ip]],output[[idraw]],start=start,end=end,nms=nms)
            }
            out[[ip]] <- list2array(tmp)
        }

    }
    names(out) <- names(params)
    return(out)
}

#' getting disaggregated prevalence from each fit
#'
#'
#'
#' @param
#' @return
#' @export
get_dis_prev <- function(setting,alldata) {
    ic <- 0
    for (draw_file in setting$draw_files) {
        ic <- ic + 1
        save_file <- get_nimble_MCMC_draws(setting$imodel,draw_file,return_type='file_path')
        tmp <- disaggregate(save_file,setting$imodel,alldata,region=setting$region)
        if (ic==1) {
            ###  initiate storage
            nsims <- tmp$nsims
            nchains <- tmp$nchains
            nareas <- dim(tmp$all)[2]
            ds <- array(0,c(setting$ndraws,nsims,nareas#alldata$nareas
                           ,setting$fit_ts[2]+setting$horizon,nchains))
            if (!is.null(setting$region)) dimnames(ds)[[3]] <- tmp$ltla
        }
        ###  store for this MI fit
        ds[ic,,,,] <- tmp$all
    }  #  next MI fit
    return(ds)
}


#' summarising disaggregated prevalence from multiple fits
#'
#'
#'
#' @param
#' @return
#' @export
summarise_dis_prev <- function(ds,setting) {
	rgn <- setting$region
	if (setting$region=='') {
	    file <- paste0('forecast/summary_nimble_multiscale_m'
                       ,setting$imodel,'_t',setting$fit_ts[1]
                       ,'-t',setting$fit_ts[2],'h',setting$horizon)
    } else {
    	file <- paste0('forecast/summary_nimble_multiscale_m'
                       ,setting$imodel,'_',rgn,'_t',setting$fit_ts[1]
                       ,'-t',setting$fit_ts[2],'h',setting$horizon)
    }
    if (setting$use_pmean) {
        file <- paste0(file,'_pmean.RData')
    } else {
        file <- paste0(file,'_',setting$ndraws,'MIs.RData')
    }
###  by each multiple imputation
    by_MI <- apply(ds,c(1,3,4),posterior_summary)

###  over all multiple imputations
    across_allMIs <- apply(ds,c(3,4),posterior_summary)
    forecast <- list(by_MI=by_MI,across_allMIs=across_allMIs,save_file=file)
    return(forecast)
}



#' carry out the preparation steps (specifying working directory, creating folders etc) for disaggregation
#'
#'
#'
#' @param
#' @return
#' @export
prepare_to_disaggregate <- function(setting) {
	#  region or national
	regions <- c('London','North East','East Midlands','East of England','South East',
                 'South West','West Midlands','Yorkshire and The Humber','North West')
	if (setting$which_region==0) {
		region <- NULL
		setting$region <- ''
	} else {
		region <- setting$region <- regions[setting$which_region]
	}
	if (setting$imodel==1) setting$use_pmean <- TRUE
	setting$recollect <- TRUE
	
	
	alldata_file <- paste0(setting$workdir_root,'alldata.rds')
	#  get all prevalence and wastewater data
	alldata <- get_all_ww_prev_data(save_file=NULL,load_file=alldata_file)


	#  define working directory (where to find the fitted files)
	if (setting$reduced_coarse | setting$region!='') {
		workdir <- paste0(setting$workdir_root,'design/')
		if (setting$reduced_coarse) {
			workdir <- paste0(setting$workdir_root,'design/subset/')
		}
	}
	setting$workdir <- workdir
	setwd(workdir)
	
	#  create folders (if not already existed) for storing results 
	if (!dir.exists('params')) dir.create('params')
	if (!dir.exists('forecast')) dir.create('forecast')

	#  get the list of the fits associated with this model (50 fits to account for uncertainty in wastewater)
	setting$draw_files <- get_nimble_files(setting)
	setting$ndraws <- length(setting$draw_files)

	#  get regional indicator if required
	if (setting$region!='') {
		alldata <- add_region(alldata)
	}
	
	out <- list(alldata=alldata,setting=setting)
	return(out)
}
