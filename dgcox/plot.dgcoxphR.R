##############################################################################################################
# Data:	8/10/2013
# Author: Luigi Augugliaro
#
# note: plot method for an object with class "dgcoxphR"
#
# Arguments
# x = the output of the dgcoxphR function;
# k = parameter used to specify the measure-of-goodness of fit. By default (k = "BIC") the BIC criterion is used; the AIC
#       is specified setting 'k == "AIC"'. The user can also use a non negative value.
# ... = arguments passed to the functions 'points' and 'matpoints'

plot.dgcoxphR <- function(x, k = c("BIC", "AIC"), ...){
	n <- dim(x$X)[1]
	if(is.numeric(k)){
		if(k <= 0) stop("k must be greater than zero")
		knm <- "GoF"
	}
	else{
		knm <- match.arg(k)
		k <- ifelse(knm == "BIC", log(n), 2)
	}
	beta <- x$beta
	logplik <- x$logplik
	g <- x$g	
	np <- ncol(beta)
	vn <- rownames(beta)
	id <- abs(beta) > 0
	mdl.list <- apply(id, 2, function(z) vn[z])
	action <- rep("", length = np)
	for(j in 2:np) action[j - 1] <- setDiff(mdl.list[[j]], mdl.list[[j - 1]])
	g.action <- g[action != ""]
    df <- x$df
    gof <- -2 * logplik + k * df
    g.gof <- g[which.min(gof)]
    plot(g, gof, xlab = expression(gamma), ylab = knm, type = "n", main = "Model Selection Criterion")
    abline(v = g.action, lty = 2, col = 8)
    abline(v = g.gof, lty = 2, col = 2)
    points(g, gof, xlab = expression(gamma), ylab = knm, type = "o", pch = 20, lty = 2, ...)
    op <- par(ask = dev.interactive())
	matplot(g, t(beta[-1, ]), col = 1, type = "n", xlab = expression(gamma), ylab = "Regression Coefficients", main = "Coefficients Path")
	abline(v = g.action, lty = 2, col = 8)
	abline(v = g.gof, lty = 2, col = 2)
	matpoints(g, t(beta[-1, ]), col = 1, type = "l", lty = 1, ...)
	axis(3, g.gof, knm, padj = 1)
	if(!is.null(g.gof)) op <- par(ask = dev.interactive())
	ru <- x$ru
	matplot(g, t(ru), col = 1, type = "n", xlab = expression(gamma), ylab = "Rao Score Statistics" ,main = "Rao Score Path")
	abline(v = g.action, lty = 2, col = 8)
	abline(v = g.gof, lty = 2, col = 2)
	matpoints(g, t(ru), col = 1, type = "l", lty = 1, ...)
	if(length(logplik)!=1) axis(3, g.gof, knm, padj = 1)
	par(op)
}
