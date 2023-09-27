#' Specification of three models in Nimble
#'
#'
#'
#' @param
#' @return
#' @export
nimble_models <- function(imodel,fitdata) {
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

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t]
        }
        alpha ~ dflat()
        sd_P ~ dunif(0,10)
        #  overall trend
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

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

######################################################################
#   a data integration model (the full model in the paper)
#   nowcasting using wastewater and national prevalence
#   with a space-time-varying ww-prevalence relationship
######################################################################
    im <- 3
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a+a[t])*(w[i,t]-8)
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + (mu_a+a[t])*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        #  RW2 as ICAR with sum2zero constraint
        a[1:ntot_times] ~ dcar_normal(adj_tm_rw2[1:K_rw2],weights_tm_rw2[1:K_rw2],num_tm_rw2[1:ntot_times],
                                      prec_a,zero_mean=1)
        prec_a <- 1/(sd_a*sd_a)
        sd_a ~ dunif(0,10)
        sd_P ~ dunif(0,10)
        #   overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })
    
######################################################################
#   a data integration model (the simplified model with BYM on U; see paper)
#   nowcasting using wastewater and national prevalence
#   with a spatial-varying ww-prevalence relationship
######################################################################
    im <- 4
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a)*(w[i,t]-8)
            }
            u[i] ~ dnorm(mu_u[i],sd=sd_u)
            mu_u[i] <- alpha + s[i]
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        ###  ICAR on the spatial random effects
        s[1:nareas] ~ dcar_normal(adj_sp[1:K_sp],weights_sp[1:K_sp],num_sp[1:nareas],prec_s,zero_mean=1)
        prec_s <- 1/(sd_s*sd_s)
        sd_s ~ dunif(0,10)

        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

######################################################################
#   the simplified model + IMD
#   nowcasting using wastewater and national prevalence
#   with a spatial-varying ww-prevalence relationship
######################################################################
    im <- 5
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a)*(w[i,t]-8) + 
                          beta.imd*std_imd[i]
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)
        
        beta.imd ~ dnorm(0,sd=1000)

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

######################################################################
#   the simplified model + IMD + BAME
#   nowcasting using wastewater and national prevalence
#   with a spatial-varying ww-prevalence relationship
######################################################################
    im <- 6
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a)*(w[i,t]-8) + 
                          beta.imd*std_imd[i] +
                          beta.bame*std_bame[i]
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)
        
        beta.imd ~ dnorm(0,sd=1000)
        beta.bame ~ dnorm(0,sd=1000)

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })


######################################################################
#   Nowcast using wastewater only, no national prevalence
######################################################################
    im <- 7
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + (gamma[i]+mu_a)*(w[i,t]-8)
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)

        #  model for the national data (not used here)
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- b[t]
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

######################################################################
#   the simplified model + BAME
#   nowcasting using wastewater and national prevalence
#   with a spatial-varying ww-prevalence relationship
######################################################################
    im <- 8
    code[[im]] <- nimbleCode({
        #  model for the LTLA data
        for (i in 1:nareas) {
            for (t in 1:ntimes) {
                p_mn[i,t] ~ dnorm(mu_pv[i,t],sd=p_sd[i,t]) # dealing with estimate uncertainty
                mu_pv[i,t] ~ dnorm(m[i,t],sd=sd_pv)  #  Type I space-time interaction
                m[i,t] <- u[i] + b[t] + (gamma[i]+mu_a)*(w[i,t]-8) + 
                          beta.bame*std_bame[i]
            }
            u[i] ~ dnorm(alpha,sd=sd_u)
            gamma[i] ~ dnorm(0,sd=sd_gamma)
        }
        sd_u ~ dunif(0,10)
        sd_gamma ~ dunif(0,10)
        sd_pv ~ dunif(0,10)
        
        beta.bame ~ dnorm(0,sd=1000)

        #  model for the national data
        for (t in 1:ntot_times) {
            P_mn[t] ~ dnorm(mu[t],sd=P_sd[t])
            mu[t] ~ dnorm(mu_P[t],sd=sd_P)
            mu_P[t] <- alpha + b[t] + mu_a*(W[t]-8)
        }
        alpha ~ dflat()
        mu_a ~ dnorm(0,sd=1000)
        sd_P ~ dunif(0,10)
        #  overall trend
        b[1] <- 0
        for (t in 2:ntot_times) {
            b[t] ~ dnorm(b[t-1],sd=sd_b)
        }
        sd_b ~ dunif(0,10)
    })

##############################
###  end model spec
##############################

    ###  specification of the chosen model
    spec <- list(code=code[[imodel]])
    ###  get initial values and parameters from the corresponding model
    tmp <- get_inits_and_params(imodel,fitdata)
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
get_final_data <- function(imodel,fitdata) {
    ###  constants
    const <- c('nareas','ntot_times','ntimes')
    if (imodel==3) {
        const <- c(const,'K_rw2','adj_tm_rw2','num_tm_rw2','weights_tm_rw2')
    }
    if (imodel==4) {
	    	const <- c(const,'K_sp','adj_sp','num_sp','weights_sp')
    }
    if (imodel==5) {
    		const <- c(const,'std_imd')
    }
    if (imodel==6) {
    		const <- c(const,'std_imd','std_bame')
    }
    if (imodel==8) {
    		const <- c(const,'std_bame')
    }
    ids <- sapply(const,function(x){which(names(fitdata)==x)})
    constants <- fitdata[c(ids)]
    fitdata <- fitdata[-c(ids)]
    
    ###  remove items that are not used in the model
    not_use <- c('ww_type','fit_ts','idraw','horizon')
    if (imodel==1) not_use <- c(not_use,'W','w')
    if (imodel==5) not_use <- c(not_use,'imd')
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
get_inits_and_params <- function(imodel,fitdata) {
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
		if (imodel==3) {
            inits1add$a <- rep(0.1,fitdata$ntot_times)
            inits1add$sd_a <- 0.1
            inits2add$a <- rep(0.2,fitdata$ntot_times)
            inits2add$sd_a <- 0.2
        }
        if (imodel==4) {
        	inits1add$s <- rep(0.1,fitdata$nareas)
        	inits1add$sd_s <- 0.1
        	inits2add$s <- rep(0.2,fitdata$nareas)
        	inits2add$sd_s <- 0.2
        }
        if (imodel==5) {
        	inits1add$beta.imd <- 0.1
        	inits2add$beta.imd <- -0.1	
        }
        if (imodel==6) {
        	inits1add$beta.imd <- 0.1
        	inits2add$beta.imd <- -0.1
        	inits1add$beta.bame <- 0.1
        	inits2add$beta.bame <- -0.1
        }
        if (imodel==8) {
        	inits1add$beta.bame <- 0.1
        	inits2add$beta.bame <- -0.1
        }
    }
    inits1 <- c(inits1,inits1add)
    inits2 <- c(inits2,inits2add)
    inits <- list(inits1,inits2)
    ###  get parameters to save
    params <- names(inits[[1]])

    out <- list(inits=inits,params=params)
    return(out)
}
