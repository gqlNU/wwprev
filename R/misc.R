#' expit transform (the inverse of logit)
#'
#'
#' @param a a logit value to be back-transformed
#' @return expit(a)=exp(a)/(exp(a)+1)
#' @export
expit <- function(a) exp(a)/(exp(a)+1)

#' Print a message on screen with a time stamp
#'
#'
#' @param m a message to print on screen
#' @return
#' @export
my_print <- function(m) {
	m <- paste0(' ',m,' @ ',Sys.time())
	print(m)
}


#' posterior summary
#'
#'
#' @param
#' @return
#' @export
posterior_summary <- function(x) {
    out <- c(quantile(x,0.025),mean(x),quantile(x,0.975))
    return(out)
}
