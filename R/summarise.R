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
  nareas <- nrow(pred_mean)
  ntimes <- ncol(pred_mean)
  cr <- array(0,c(nareas,ntimes))
  for (i in 1:nareas) {
    for (j in 1:ntimes) {
      cr[i,j] <- point_in_interval(obs[i,j],pred_ci[,i,j])
    }
  }
  cover <- mean(cr)
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
