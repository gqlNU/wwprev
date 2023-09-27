#' Fit a model in Nimble
#'
#'
#'
#' @param
#' @return
#' @export
fit_nimble <- function(imodel,iters,fitdata,reduced_coarse=FALSE,region=NULL) {
    fit_ts <- fitdata$fit_ts
    horizon <- fitdata$horizon
    if (!reduced_coarse) {
    	#  coarse-level prevalence available for all weeks during nowcast
	    ###  model specification
   		model_spec <- nimble_models(imodel,fitdata)
    	###  specify constants in the model
	    fitdata <- get_final_data(imodel,fitdata)
	} else {
		#  coarse-level prevalence available ONLY for some weeks during nowcast
	    ###  model specification
   		model_spec <- nimble_models_subset(imodel,fitdata)
    	###  specify constants in the model
	    fitdata <- get_final_data_subset(imodel,fitdata)		
	}
    ###  nimble fitting
    msg <- paste0(' [*] fitting m',imodel,' t',fit_ts[1],'-t',fit_ts[2],'_horizon',horizon)
    my_print(msg)
    output <- nimbleMCMC(code=model_spec$code
                        ,data=fitdata$fitdata
                        ,constants=fitdata$constants
                        ,inits=model_spec$inits
                        ,monitors=model_spec$params
                        ,niter=iters[1],nburnin=iters[2],thin=iters[3]
                        ,nchains=length(model_spec$inits))
    return(output)
}

#' Save the nimble fit to file
#'
#'
#'
#' @param
#' @return
#' @export
save_nimble_fit <- function(imodel,output,save_dir,fitdata,iters,region=NULL) {
    fit_ts <- fitdata$fit_ts
    horizon <- fitdata$horizon
    idraw <- fitdata$idraw
    msg <- paste0(' [*] saving m',imodel,' t',fit_ts[1],'-t',fit_ts[2],'_horizon',horizon)
    my_print(msg)
    file_end <- paste0('idraw_ww',idraw,'.RData')
    if (idraw==0 | imodel==1) {
        file_end <- '_pmean.RData'
    }
    save_file <- paste0(save_dir,'nimble_multiscale_m',imodel,
                        '_t',fit_ts[1],'-t',fit_ts[2],'h',horizon,file_end)
    if (!is.null(region)) {
	    save_file <- paste0(save_dir,'nimble_multiscale_m',imodel,'_',region,
    	                    '_t',fit_ts[1],'-t',fit_ts[2],'h',horizon,file_end)
    }
    save(file=save_file,output,fitdata,iters)
    return()
}
