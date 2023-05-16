
#' Extract prevalence and wastewater data from the entire set for a specific model fit
#'
#'
#'
#' @param
#' @return
#' @export
extract_from_alldata <- function(imodel,fit_ts,horizon,idraw,alldata) {
    t_ids_short <- fit_ts[1]:fit_ts[2]     #  the length of time over which ltla prevalence is available
    t_ids <- fit_ts[1]:(fit_ts[2]+horizon) #  the entire set of time points
    out <- list(ntot_times=length(t_ids),nareas=dim(alldata$w_sims)[2])
    out$fit_ts <- fit_ts
    out$idraw <- idraw
    out$horizon <- horizon
    out$ntimes <- length(t_ids_short)
    #=============
    #  prevalence (posterior mean + SD)
    #
    #  national prevalence
    out$P_mn <- alldata$P_mn[t_ids]
    out$P_sd <- alldata$P_sd[t_ids]
    #  LTLA prevalence
    out$p_mn <- alldata$p_mn[,t_ids_short]
    out$p_sd <- alldata$p_sd[,t_ids_short]
    #==============
    #  wastewater: from a specific draw (idraw>0) or posterior means (idraw=0)
    #
    out$ww_type <- paste0('from draw ',idraw)
    out$W <- alldata$W_sims[idraw,t_ids]
    out$w <- alldata$w_sims[idraw,,t_ids]
    if (idraw==0) {
        #  use the posterior means for W and w
        out$W <- alldata$W_mn[t_ids]
        out$w <- alldata$w_mn[,t_ids]
        out$ww_type <- 'using posterior mean'
    }
    #  add RW neighbourhood structure
    if (imodel==3) {
        out <- add_tmajd_to_fitdata(out,nts=out$ntot_times,RW.order=2)
    }
    
    #  add spatial neighbourhood structure
    if (imodel==4) {
        out <- add_spajd_to_fitdata(out)
    }
    #  add IMD and BAME
    if (imodel==5) {
        out <- add_IMD(out)
        out <- add_ethnicity(out)
    }
    return(out)
}

#' Get all LTLA/national weekly prevalence and WW data across all 43 weeks
#'
#'
#'
#' @param
#' @return
#' @export
get_all_ww_prev_data <- function(save_file=NULL,load_file=NULL) {
    ww_data_dir <- '/Volumes/WorkSpace/wwprev/inst/extdata/' #system.file("extdata", "ww_model_fit-20220801-195150_CTRYpred_pw.RData", package = "wwprev")
    if (!is.null(load_file)) {
        if (!file.exists(load_file)) stop(' Error in get_all_ww_prev_data: cannot find the alldata file, construct it first!')
        my_print('loading alldata')
###  all the data have been constructed so just load the data
        alldata <- readRDS(load_file)
        my_print('alldata loaded')
    } else {
        my_print('start constructing alldata')
###  construct all the data from scratch
###  get LTLA-weekly prevalence
        tmp <- get_ltla_weekly_debiased_prevalence()
        nts_prev <- tmp$nts_prev
        ltla_prev <- tmp$ltla_prev

###  get LTLA-weekly WW estimates
        ww_ltla_weekly_sims <- get_ltla_weekly_ww(nts_prev,ww_data_dir)
        nsims <- dim(ww_ltla_weekly_sims)[1]
        ltla_ids_ww <- dimnames(ww_ltla_weekly_sims)[[2]]
        weekly_ids_ww <- dimnames(ww_ltla_weekly_sims)[[3]]

###  get national-weekly WW estimates
        tmp <- get_national_WW(nts_prev,ww_data_dir)
        ww_national_summary <- tmp$ww_national_summary
        ww_national_sims <- tmp$ww_national

###  format LTLA-weekly prevalence
        prev_wide <- format_ltla_prevalence_wide_sims(ltla_prev,nsims,nts_prev,ltla_ids_ww,weekly_ids_ww,get_sims=TRUE)

		#  get population weighted national-weekly prevalence
		prev_national <- get_national_prevalence(prev_wide)

		#  putting everything together
		nareas <- dim(ww_ltla_weekly_sims)[2]
		ntot_times <- dim(ww_ltla_weekly_sims)[3]
		alldata <- list(ntot_times=ntot_times,nareas=nareas,nsims=nsims)
		#  wastewater data
		alldata$W_sims <- ww_national_sims
    alldata$W_mn <- apply(alldata$W_sims,2,mean)
    alldata$w_sims <- ww_ltla_weekly_sims
        alldata$w_mn <- apply(alldata$w_sims,c(2,3),mean)

    #  prevalence data
		alldata$P_sims <- prev_national$sims
                alldata$P_mn <- prev_national$mean
                alldata$P_sd <- sqrt(prev_national$var)
		alldata$p_sims <- prev_wide$sims
                alldata$p_mn <- prev_wide$mean
                alldata$p_sd <- sqrt(prev_wide$var)

		if (!is.null(save_file)) {
                    saveRDS(alldata,file=save_file)
                    my_print('alldata saved')
		}
	}
  return(alldata)
}


#' Import LTLA-weekly debiased prevalence
#'
#' This function loads the LTLA-weekly debiased prevalence
#'
#' @param file Path to the data file
#' @return debiased prevalence
#' @export
get_ltla_weekly_debiased_prevalence <- function(file='/Volumes/WorkSpace/wwprev/inst/extdata/logit_moments.rds') {
    prev <- readRDS(file=file)
    nts_prev <- max(prev$time_index_ww)
    out <- list(ltla_prev=prev,nts_prev=nts_prev)
}


#' Import LTLA-weekly wastewater estimates
#'
#' This function loads the LTLA-weekly wastewater estimates
#'
#' @param nts Number of weeks
#' @return
#' @export
get_ltla_weekly_ww <- function(nts,ww_data_dir
                              ,ww_model='ww_model_fit-20220801-195150' #  regional random effect
                              ,remove_3LTLAs=TRUE) {
    cw <- getwd() #  current working directory
    setwd(ww_data_dir)
    ###  load WW
    aggregate.method <- 'pw'
    pred <- publicWW::load_pred(ww_model,'LAD',aggregate.method)
    setwd(cw)  # change work directory to the current one
    lad21cd <- dimnames(pred)[[2]]
    nlads <- length(lad21cd)

    ###  remove City of London and Isle of Scilly and Carlisle
    if (remove_3LTLAs) {
        lad.ids <- (1:nlads)[-sapply(c('E06000053','E09000001','E07000028'),function(x){which(lad21cd==x)})]
        final_pred <- pred[,lad.ids,]
    } else {
        final_pred <- pred
    }

    ###  match the time span with prevalence
    final_pred <- final_pred[,,1:nts]

    ###  assign weeks
    week.start <- read.csv(system.file("extdata", "time_index.csv",package = "publicWW"))
    dimnames(final_pred)[[3]] <- week.start$week_start[1:dim(final_pred)[3]]
    return(final_pred)
}


#' Import weekly wastewater estimates for the whole of England
#'
#' This function loads the weekly wastewater estimates for the whole of England. These are population weighted average of the
#' LSOA-level weekly wastewater estimates
#'
#' @param nts Number of weeks
#' @return
#' @export
get_national_WW <- function(nts,ww_data_dir
                           ,ww_model='ww_model_fit-20220801-195150') {
    cw <- getwd()  #  current working directory
    setwd(ww_data_dir)
    ###  load national WW
    aggregate.method <- 'pw'
    pred_national <- publicWW::load_pred(ww_model,'CTRY',aggregate.method)
    setwd(cw) #  change back to the current work directory

    ww_national <- pred_national[,1,1:nts]
    ww_national_summary <- as.data.frame(t(apply(ww_national,2,function(x){c(mean(x),var(x))})))
    names(ww_national_summary) <- c('mean','var')

    week.start <- read.csv(system.file("extdata", "time_index.csv",package = "publicWW"))
    ww_national_summary$time_index <- week.start$week_start[1:nrow(ww_national_summary)]

	out <- list(ww_national_summary=ww_national_summary
	           ,ww_national=ww_national)
    return(out)
}

#' Import weekly wastewater estimates for the whole of England
#'
#' This function loads the weekly wastewater estimates for the whole of England. These are population weighted average of the
#' LSOA-level weekly wastewater estimates
#'
#' @param nts Number of weeks
#' @return
#' @export
format_ltla_prevalence_wide_sims <- function(ltla_prev,nsims,nts_prev,ltla_ids_ww,weekly_ids_ww,get_sims=TRUE) {
    nareas <- length(ltla_ids_ww)
###  format debiased prev
    prev_mn <- prev_var <- array(0,c(nareas,nts_prev))
    prev_sims <- NULL
    if (get_sims) {  #  simulate ltla prevalence based on the uncertainty
        prev_sims <- array(0,c(nsims,nareas,nts_prev))
    }
    for (i in 1:nareas) {
        this_lad <- ltla_ids_ww[i]
        for (tt in 1:nts_prev) {
            id <- which(ltla_prev$LAD21CD==this_lad & ltla_prev$time_index_ww==tt)
            if (length(id)>0) {
            	m <- ltla_prev[['mean']][id]
            	s <- ltla_prev[['sd']][id]
                prev_mn[i,tt] <- m
                prev_var[i,tt] <- s^2
                if (get_sims) prev_sims[,i,tt] <- rnorm(nsims,m,s)  #  draws for prevalence
            }
        }
    }
    rownames(prev_mn) <-
        rownames(prev_var) <- ltla_ids_ww
    colnames(prev_mn) <-
      colnames(prev_var) <- weekly_ids_ww
    if (get_sims) {
        dimnames(prev_sims)[[2]] <- ltla_ids_ww
        dimnames(prev_sims)[[3]] <- weekly_ids_ww
    }
    out <- list(mean=prev_mn,var=prev_var,sims=prev_sims)
    return(out)
}

#' Read in LTLA level population
#'
#' data source: https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
#'
#' @param
#' @return
#' @export
get_LTLA_population <- function() {
	d <- readxl::read_excel('/Volumes/WorkSpace/wwprev/inst/extdata/ukpopestimatesmid2021on2021geographyfinal.xls',sheet='MYE2 - Persons',skip=7)[c('Code','Name','All ages')]
	return(d)
}

#' Obtain weekly national prevalence by taking population-weighted average of the LTLA-weekly prevalence
#'
#'
#'
#' @param
#' @return
#' @export
get_national_prevalence <- function(prev_wide) {
	#  get and format LTLA population
    pop <- get_LTLA_population()

	#  match the ltla in pop with that in prev_wide
    ltla <- rownames(prev_wide$mean)
    ids <- sapply(ltla,function(x){which(pop$Code==x)})
    spop <- pop[ids,]
    totpop <- sum(spop[['All ages']])

    sims <- prev_wide$sims
    prev_sims <- apply(sims,c(1,3),function(x) {
        cases <- expit(x)*spop[['All ages']]
        logit_pv <- car::logit(sum(cases)/totpop)
    })
    mn <- apply(prev_sims,2,mean)
    v <- apply(prev_sims,2,var)
    out <- list(mean=mn,var=v,sims=prev_sims)
    return(out)
}


#' Add RW neighbourhood structure to fitdata
#'
#'
#'
#' @param
#' @return
#' @export
add_tmajd_to_fitdata <- function(fitdata,nts,RW.order=1) {
    tmadj <- tm_adjacency(nts,RW.order=RW.order)
    if (RW.order==1) {
        fitdata$adj_tm_rw1 <- tmadj$adj
        fitdata$num_tm_rw1 <- tmadj$num
        fitdata$weights_tm_rw1 <- tmadj$weights
        fitdata$K_rw1 <- length(fitdata$adj_tm_rw1)
    }
    if (RW.order==2) {
        fitdata$adj_tm_rw2 <- tmadj$adj
        fitdata$num_tm_rw2 <- tmadj$num
        fitdata$weights_tm_rw2 <- tmadj$weights
    	fitdata$K_rw2 <- length(fitdata$adj_tm_rw2)
    }
    return(fitdata)
}


#' Construct RW1/2 neighbourhood structure
#'
#'
#'
#' @param
#' @return
#' @export
tm_adjacency <- function(T,RW.order) {
################################################################################
#    create the temporal adjacency matrix for the random walk prior
#    the codes were adapted from the GeoBUGS manual (V1.2 Example 9)
#
#    input:
#=====================
#       T 			= number of time points
#       RW.order	= 1 or 2 (order of the random walk prior)
#    output:
#=====================
#       a list containing three components: adj, num and weights
#
#       num:	 an array of size T with Entry i showing number of temporal neighbours
#            	 for time point i;
#       adj: 	 an array showing the neighbours of each time point;
#       weights: an array specifying the weights of the temporal
#                neighbours.
################################################################################
    if (RW.order==1) {
####   for random walk order 1
        num <- numeric(T)
        weights <- numeric((T-2)*2+2)
        adj <- numeric((T-2)*2+2)
        for (t in 1:1) {
            weights[t] <- 1
            adj[t] <- t+1
            num[t] <- 1
        }
        for (t in 2:(T-1)) {
            weights[2+(t-2)*2] <- 1
            adj[2+(t-2)*2] <- t-1
            weights[3+(t-2)*2] <- 1
            adj[3+(t-2)*2] <- t+1
            num[t] <- 2
        }
        for (t in T:T) {
            weights[(T-2)*2+2] <- 1
            adj[(T-2)*2+2] <- t-1
            num[T] <- 1
        }
    }
    if (RW.order==2) {
  ####   for random walk order 2
        num <- numeric(T)
        weights <- numeric((T-4)*4+3*2+2*2)
        adj <- numeric((T-4)*4+3*2+2*2)
        for (t in 1:1) {
            weights[t] <- 2
            adj[t] <- t+1
            weights[t+1] <- -1
            adj[t+1] <- t+2
            num[t] <- 2
        }
        for (t in 2:2) {
            weights[t+1] <- 2
            adj[t+1] <- t-1
            weights[t+2] <- 4
            adj[t+2] <- t+1
            weights[t+3] <- -1
            adj[t+3] <- t+2
            num[t] <- 3
        }
        for (t in 3:(T-2)) {
            weights[6+(t-3)*4] <- -1
            adj[6+(t-3)*4] <- t-2
            weights[7+(t-3)*4] <- 4
            adj[7+(t-3)*4] <- t-1
            weights[8+(t-3)*4] <- 4
            adj[8+(t-3)*4] <- t+1
            weights[9+(t-3)*4] <- -1
            adj[9+(t-3)*4] <- t+2
            num[t] <- 4
        }
        for (t in (T-1):(T-1)) {
            weights[(T-4)*4+6] <- 2
            adj[(T-4)*4+6] <- t+1
            weights[(T-4)*4+7] <- 4
            adj[(T-4)*4+7] <- t-1
            weights[(T-4)*4+8] <- -1
            adj[(T-4)*4+8] <- t-2
            num[t] <- 3
        }
        for (t in T:T) {
            weights[(T-4)*4+9] <- 2
            adj[(T-4)*4+9] <- t-1
            weights[(T-4)*4+10] <- -1
            adj[(T-4)*4+10] <- t-2
            num[t] <- 2
        }
    }
    if (T==2) {
        adj <- c(2,1)
        weights <- c(1,1)
        num <- c(1,1)
    }
    adj.tm <- list()
    adj.tm$adj <- adj
    adj.tm$num <- num
    adj.tm$weights <- weights
    return(adj.tm)
}


#' Construct the W matrix for LTLAs in England
#'
#'
#'
#' @param
#' @return
#' @export
add_spajd_to_fitdata <- function(dat) {
	data_file <- system.file("extdata", "Local_Authority_Districts_December_2021_GB_BUC_2022_7123482767279458595/LAD_DEC_2021_GB_BUC.shp", package = "wwprev")
	ltla <- sf::st_read(data_file)
	ltla_cd <- rownames(dat$w)
	ids <- sapply(ltla_cd,function(x){which(ltla$LAD21CD==x)})
	sltla <- ltla[ids,]
	m <- spdep::poly2nb(sltla)
	###  linking Isle of Wight to Southampton
	id_wight <- which(sltla$LAD21NM=='Isle of Wight')
	id_southampton <- which(sltla$LAD21NM=='Southampton')
	m[[id_wight]] <- as.integer(id_southampton)
	m[[id_southampton]] <- c(m[[id_southampton]],as.integer(id_wight))
#	w <- spdep::nb2mat(m, style = "B")
	spadj <- spdep::nb2WB(m)
	names(spadj$num) <- sltla$LAD21CD
	dat$adj_sp <- spadj$adj
	dat$num_sp <- spadj$num
	dat$weights_sp <- spadj$weights
	dat$K_sp <- length(spadj$adj)
	return(dat)
}


#' add ethnicity
#'
#'
#'
#' @param
#' @return
#' @export
add_ethnicity <- function(dat) {
	data_file <- system.file("extdata", "ethnic2021.xlsx", package = "wwprev")
	d <- as.data.frame(readxl::read_excel(data_file,sheet='Figure 3',skip=4))
	###  get the number columns
	all <- d[,grep('number',colnames(d))]
	n <- apply(all,1,sum)
	###  get the white columns
	white <- all[,grep('White:',colnames(all))]
	nw <- apply(white,1,sum)
	###  non-white percentage
	bame <- data.frame(bame=(n-nw)/n*100,LAD21CD=d[['Area code']],LAD21NM=d[['Area name']])
	###  select and reorder the LTLAs
	ids <- sapply(rownames(dat$w),function(x){which(bame$LAD21CD==x)})
	dat$bame <- bame$bame[ids]
	dat$std_bame <- scale(bame$bame[ids])[,1]
	return(dat)
}



#' add IMD (boundary change https://www.gov.uk/government/statistics/2011-rural-urban-classification-lookup-tables-for-all-geographies: Rural Urban Classification 2011 lookup tables for local authority areas)
#'
#'
#'
#' @param
#' @return
#' @export
add_IMD <- function(dat) {
	data_file <- system.file("extdata", "File_10_-_IoD2019_Local_Authority_District_Summaries__lower-tier__.xlsx", package = "wwprev")
	d <- as.data.frame(readxl::read_excel(data_file,sheet='IMD'))
	score <- d[['IMD - Average score']]
	
	new_northampton <- c('E06000060','E06000061','E06000062') # boundary change from 2019 to 2021
	###  select and reorder the LTLAs
	imd <- NULL
	lad21cd <- rownames(dat$w)
	n <- length(lad21cd)
	for (ilad in 1:n) {
		lad <- lad21cd[ilad]
		if (!any(new_northampton==lad)) {
			ii <- which(d['Local Authority District code (2019)']==lad)
			imd <- c(imd,score[ii])
		} else {
			if (lad=='E06000060') x <- c('E07000004','E07000005','E07000006','E07000007')
			if (lad=='E06000061') x <- c('E07000150','E07000152','E07000153','E07000156')
			if (lad=='E06000062') x <- c('E07000151','E07000154','E07000155')
			ii <- sapply(x,function(xx){which(d['Local Authority District code (2019)']==xx)})
			imd <- c(imd,mean(score[ii]))
		}
	}
	dat$imd <- imd
	dat$std_imd <- (scale(imd))[,1]
	return(dat)
}
