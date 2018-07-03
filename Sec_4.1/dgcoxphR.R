###################################################################################################################
# Version: 2.3
# Data:	8/10/2013
# Author: Luigi Augugliaro - University of Palermo (Italy)
# e-mail: luigi.augugliaro@unipa.it
#
# Arguments
# time = The follow up time.
# status = The status indicator. 0 means alive and 1 means dead.
# X = Design matrix of dimension n x p.
# g0 = The smallest value for the tuning parameter gamma. Default is g0 = ifelse(p < n, 1.0e-04, 0.05)
# maxrsq_ratio = Argument used to define a second stop criterion. The algorithm is stoped when the ratio
#				 between the Cox-Snell pseudo-R square and its theoretical maximum value is below maxrsq_ratio.
# np = Number of point of the solution curve.
# nstep = Number of Newton-Raphson iterations.
# eps = Agument used for the convergence of the Newton-Raphson algorithm.
# verbose = Flag used to print out information about Newton-Raphson algorithm (only for debug).

dgcoxphR <- function(time, status, X, g0 = NULL, maxrsq_ratio = 0.9, np = 500, nstep = 100, eps = 1e-13, verbose = FALSE){
	index <- order(time, decreasing = TRUE)
	Xs <- X[index, ]
	Xs2 <- Xs^2
	nx <- dim(Xs)[1]
	p <- dim(Xs)[2]
	if(is.null(colnames(Xs))) colnames(Xs) <- paste("x", 1:p, sep = "")
	if(is.null(g0)) g0 <- ifelse(p < nx, 1.0e-04, 0.05)
	if(g0==0) g0 <- 1e-6
	D <- which(status[index] == 1)
	nd <- length(D)
	ru <- beta <- matrix(0, p, np)
	rownames(ru) <- rownames(beta) <- colnames(Xs)
	g <- vector(length =np, mode = "numeric")
	logplik <- vector(length = np, mode = "numeric")
	dev <- vector(length = np, mode = "numeric")
	dev_ratio <- vector(length = np, mode = "numeric")
	rsq <- vector(length = np, mode = "numeric")
	steps <- vector(length = np, mode = "numeric")
	df <- vector(length = np, mode = "numeric")
	gic_compl <- vector(length = np, mode = "numeric")
	id <- 1:p
	Y <- apply(Xs[D, ], 2, sum)
	exp_eta <- rep(1, nx)
	wgh <- 1:nx
	Z <- apply(Xs, 2, cumsum)
	Z <- Z / wgh
	Cmm <- apply(Xs2, 2, cumsum)
	Cmm <- Cmm / wgh
	dl <- Y - apply(Z[D, ], 2, sum)
	s <- sign(dl)
	sqrtIm <- sqrt(apply(Cmm[D, ] - Z[D, ]^2, 2, sum))
	k <- 1
	ru[, k] <- dl / sqrtIm
	g[k] <- max(abs(ru[, k]))
	A <- which(abs(ru[, k]) == g[k])
	nav <- length(A)
	logplik[k] <- - sum(log(wgh[D]))
	maxrsq <- 1 - exp(2 * logplik[k] / nx)
	dg <- exp((log(g0) - log(g[1]))/(np - 1))
	df[k] <- 0
	ag <- g[k]
	new <- TRUE
	repeat{
		if(new){
			k <- k + 1
			if(k > np) break
			ag = ag * dg
			ba <- beta[, k - 1]
		} else ba <- beta[, k]
		ba[-A] <- 0
		Xsa <- Xs[, A, drop = FALSE]
		Xs2a <- Xs2[,A,drop = FALSE]
		Fim <- matrix(0, nav, nav)
		dnImm <- matrix(0, nav, nav)
		if(k > np) break
		for(i in 1:nstep){
			Za <- apply(Xsa, 2, function(z) cumsum(exp_eta * z))
			Za <- Za / wgh
			for(m in 1:nav){
				Cmm <- cumsum(exp_eta * Xs2a[, m])
				Cmm[D] <- Cmm[D] / wgh[D]
				for(n in 1:nav){
					Cmn <- cumsum(exp_eta * Xsa[, m] * Xsa[, n])
					Cmn[D] <- Cmn[D] / wgh[D]
					Cmmn <- cumsum(exp_eta * Xs2a[, m] * Xsa[, n])
					Cmmn[D] <- Cmmn[D] / wgh[D]
					Imn <- Cmn - Za[, m] * Za[, n]
					Fim[m, n] <- sum(Imn[D])
					dnImm[m, n] <- sum(Cmmn[D] - Cmm[D] * Za[D, n]) - 2 * sum(Imn[D] * Za[D,m])
				}
			}
			dl[A] <- Y[A] - apply(Za[D, , drop = FALSE], 2, sum)
			sqrtIm[A] <- sqrt(diag(Fim))
			rua <- dl[A] / sqrtIm[A]
			dba <- rua - ag * s[A]
			if(verbose) cat("||dba||^2 =", sum(dba^2), "g =", ag, "\n")
			if(is.nan(sum(dba^2))) break
			Jmn <- matrix(0.5 * rua / diag(Fim), nav, nav) * dnImm + Fim / matrix(sqrtIm[A], nav, nav)
			dba <- try(solve(Jmn,dba), silent = TRUE)
			if(class(dba) == "try-error")	break
			ba[A] <- ba[A] + dba
			eta <- drop(Xsa %*% ba[A])
			exp_eta <- exp(eta)
			wgh <- cumsum(exp_eta)
            if(sum(dba^2) <= eps) break
		}
		if(i == nstep | class(dba) == "try-error" | is.nan(sum(dba^2))){
			warning("NR does not converge at gamma = ", ag, "last values are reported")
			break
		}
		new <- TRUE
		Z <- apply(Xs, 2, function(z) cumsum(exp_eta * z))
		Z <- Z / wgh
		Cmm <- apply(Xs2, 2, function(z) cumsum(exp_eta * z))
		Cmm <- Cmm / wgh
		dl <- Y - apply(Z[D, ], 2, sum)
		s <- sign(dl)
		sqrtIm <- sqrt(apply(Cmm[D, ] - Z[D, ]^2, 2, sum))
		rum <- dl / sqrtIm
		R <- Fim + 0.5 * ag * dnImm * matrix(s[A] / sqrtIm[A], nav, nav)
		dli <- Xs[D, A] - Z[D, A]
		Q <- crossprod(dli) - ag * outer(s[A] * sqrtIm[A] / nd, dl[A], "*")
		Rinv_Q <- try(solve(R,Q), silent = TRUE)
		if(class(Rinv_Q) == "try-error"){
			warning("GIC complexity cannot be computed at gamma = ", ag, "last values are reported")
			break
		}
		id_new <- which(abs(rum[-A]) >= ag)
		A_new <- sort(c(A, id[-A][id_new]))
		if(length(A_new) != length(A)){
			A <- A_new
			nav <- length(A)
			new <- FALSE
		}
		steps[k] <- i
		g[k] <- ag
		df[k] <- nav
		gic_compl[k] <- sum(diag(Rinv_Q))
		beta[A, k] <- ba[A]
		ru[, k] <- rum
		logplik[k] <- sum(eta[D]) - sum(log(wgh[D]))
		dev[k] <- 2 * (logplik[k] - logplik[1])
		dev_ratio[k] <- 1 - logplik[k]/logplik[1]
		rsq[k] <- 1 - exp(-dev[k] / nx)
		if(verbose)
            cat("nav =", nav, "np =", k, "dev_ratio =", round(dev_ratio[k], 3),
            "R2 =", round(rsq[k], 3), "max R2 =", round(maxrsq, 3), "completed\n\n")
		if(rsq[k] / maxrsq > maxrsq_ratio) {
			k <- k + 1
			break
		}
	}
	np <- k - 1
	beta <- beta[, 1:np]
	ru <- ru[, 1:np]
	g <- g[1:np]
	logplik <- logplik[1:np]
	dev <- dev[1:np]
	dev_ratio <- dev_ratio[1:np]
	rsq <- rsq[1:np]
	df <- df[1:np]
	gic_compl <- gic_compl[1:np]
	steps <- steps[1:np]
	out <- list(beta = beta, ru = ru, g = g, logplik = logplik, dev = dev, dev_ratio = dev_ratio, rsq = rsq,
                maxrsq = maxrsq, df = df, gic_compl = gic_compl, steps = steps, time = time, status = status,
                g0 = g0, np = np, nstep = nstep, eps = eps, X = X)
	class(out) <- "dgcoxphR"
	out
}
