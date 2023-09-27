#' Specification of models in Nimble with spatially-coarse prevalence data only available for a subset of weeks during nowcast
#'
#'
#'
#' @param
#' @return
#' @export
nimble_models_subset <- function(imodel,fitdata) {
    code <- as.list(1:20)
######################################################################
#   using national prevalence only (no wastewater)
######################################################################
    im <- 1
    code[[im]] <- nimbleCode({
      #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t]
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
        }
        sd_u ~ dunif(0,10)
        sd_pv ~ dunif(0,10)

		#  model for the national data (during training period)
        for (t in 1:ntimes) {
      		P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
        }
        #  model for the national data (during nowcast period)
        for (t in (ntimes+1):(ntimes+nd_times)) {
	        P_mn[t] ~ dnorm(mu[coarse_avail_times[t-ntimes]],sd=P_sd[t])  #  P_mn only available for a subset of weeks
     	}
     	
        for (t in 1:ntot_times) {  #  defining the underlying prevalence process for all weeks 
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t]
        }
        alpha ~ dflat()
        sd_P ~ dunif(0,10)
        #  overall trend for ALL weeks
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })


######################################################################
#   a data integration model (the simplified model in the paper)
#   nowcasting using wastewater and national prevalence
#   with a spatial-varying ww-prevalence relationship
######################################################################
    im <- 2
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a)*(w[i,t]-8)
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)

        #  model for the national data (during training period)
        for (t in 1:ntimes) {
      		P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
        }
        #  model for the national data (during nowcast period)
        for (t in (ntimes+1):(ntimes+nd_times)) {
	        P_mn[t] ~ dnorm(mu[coarse_avail_times[t-ntimes]],sd=P_sd[t])  #  P_mn only available for a subset of weeks
     	}
        for (t in 1:ntot_times) {  #  defining the underlying prevalence process for all weeks 
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend for ALL weeks
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

################################
###  end model spec
################################

    ###  specification of the chosen model
    spec <- list(code=code[[imodel]])
    ###  get initial values and parameters from the corresponding model
    tmp <- get_inits_and_params_subset(imodel,fitdata)
    spec <- c(spec,tmp)
    return(spec)
}

#' get constants for nimble fit
#'
#'
#'
#' @param
#' @return
#' @export
get_final_data_subset <- function(imodel,fitdata) {
  
    ###  constants
    const <- c('nareas','ntot_times','ntimes','nd_times','coarse_avail_times')
    ids <- sapply(const,function(x){which(names(fitdata)==x)})
    constants <- fitdata[c(ids)]
    fitdata <- fitdata[-c(ids)]
    ###  remove items that are not used in the model
    not_use <- c('ww_type','fit_ts','idraw','horizon')
    ids <- sapply(not_use,function(x){which(names(fitdata)==x)})
    fitdata <- fitdata[-c(ids)]
    out <- list(fitdata=fitdata,constants=constants)
    return(out)
}


#' Specify initial values and parameters to save
#'
#'
#'
#' @param
#' @return
#' @export
get_inits_and_params_subset <- function(imodel,fitdata) {
	inits1 <- list(mu=rep(-6,fitdata$ntot_times)
                  ,sd_P=0.2
                  ,alpha=-6
                  ,b=c(NA,rep(0.1,fitdata$ntot_times-1))
                  ,sd_b=0.1
                  ,mu_pv=array(-6,c(fitdata$nareas,fitdata$ntimes))
                  ,sd_pv=0.1
                  ,u=rep(0.1,fitdata$nareas)
                  ,sd_u=0.5
                   )
    inits2 <- list(mu=rep(-7,fitdata$ntot_times)
                  ,sd_P=0.5
                  ,alpha=-7
                  ,b=c(NA,rep(0.2,fitdata$ntot_times-1))
                  ,sd_b=0.2
                  ,mu_pv=array(-5,c(fitdata$nareas,fitdata$ntimes))
                  ,sd_pv=0.2
                  ,u=rep(0.2,fitdata$nareas)
                  ,sd_u=0.2
                   )
    inits1add <- inits2add <- NULL
    if (imodel>1) {
        inits1add <- list(mu_a=0.1
                         ,gamma=rep(0.1,fitdata$nareas)
                         ,sd_gamma=0.2)
        inits2add <- list(mu_a=0.2
                         ,gamma=rep(0.2,fitdata$nareas)
                         ,sd_gamma=0.1)
	}    
    inits1 <- c(inits1,inits1add)
    inits2 <- c(inits2,inits2add)
    inits <- list(inits1,inits2)
    ###  get parameters to save
    params <- names(inits[[1]])

    out <- list(inits=inits,params=params)
    return(out)
}
