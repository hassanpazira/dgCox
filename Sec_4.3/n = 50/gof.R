###############################################################################################################
# Data:	8/10/2013
# Author: Luigi Augugliaro
#
# Description: The 'gof' function is used to compute the following measure-of-goodness of fit.
#
# -2 logLikelihood + k * pen
#
# where 'pen' is the quantity used to measure the complesity of the fitted model
# and 'k' is the associated weight.
#
# Arguments
#
# object = the output of the dgcoxphR function;
# k = weight of the completity term
# complexity = complexity used in the measure of goodness-of-fit. If 'complexity = df' (default)
# the complefity of the fitted model is measured as the number of non-zero estimates otherwise
# ('complexity = df') the GIC criterion is used.

gof <- function(object, k = c("BIC", "AIC"), complexity = c("df", "gic")){
	nd <- sum(object$status == 1)
	if(is.character(k)){
		k <- match.arg(k)
		k <- ifelse(k == "BIC", log(nd), 2)
	}
	if(is.numeric(k) & k < 0) stop("k shoud be a positive value")
	complexity <- match.arg(complexity)
	if(complexity == "df") pen <- object$df
	if(complexity == "gic") pen <- object$gic_compl
	gof <- -2 * object$logplik + k * pen
	gof
}
