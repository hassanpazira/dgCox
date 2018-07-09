################################################################################################################
# Author: F. Abegaz
#
# Data: 2017-02-20
#
# Note: This simulation study is intended to analyse the behaviour of different measures of goodness-of-fit (gof)
# We considere the following gof:
#
# AIC = classical Akaike Information Criterion
# BIC = classical Bayesian Information Criterion
# Fan_13 = the BIC adjusted with the correction term proposed in JRSS-B (2013) Vol 75(3) pp. 531-552
# GIC_aic = the Generalized Information Criterion proposed in Biometrika (1996) Vol 83(4) pp. 875-890 with k = 2 and complexity = tr(R^{-1}Q)
# GIC_bic = the Generalized Information Criterion proposed in Biometrika (1996) Vol 83(4) pp. 875-890 with k = log(n) and complexity = tr(R^{-1}Q)
# GIC_bic_modified = the Generalized Information Criterion proposed in Biometrika (1996) Vol 83(4) pp. 875-890 with k = log(log(n)) * log(p) and complexity = tr(R^{-1}Q)
################################################################################################################

########################
## Libraries and R code
source("dgcoxphR.R")
source("gof.R")
library(MASS)

########################
## general setting

nsim <- 500                                 # number of simulations
n <- 50                                     # sample size
p <- c(50, 100, 1000)                       # number of variables
s <- 2^c(1, 3, 5)                           # number of non-zero regression coefficients
bs <- 0.5                                   # value of the non-zero regression coefficients
pi <- 0.2                                   # censoring probability
rho <- 0.9                                  # parameter of the autoregressive covariance matrix
save.point <- seq(100, nsim, by = 100)      # save point for tmp files
g0 <- c(0.5, 0.5, 1.5)

##########################
## Output
results <- array(0, dim = c(nsim, length(p), length(s), 5, 11),
	dimnames = list(nsim = 1:nsim, p = p, s = s,
	gof = c("AIC", "BIC", "Fan_13", "GIC_aic", "GIC_bic"),
	measures = c("value", "size", "TP", "FP", "TN", "FN", "TPR", "FPR", "FDR", "FNR", "F1")))

############################
## Starting simulation study

for(h in 1:length(p)){
    
    # parameters of the multivariate Gaussian distribution
    mu <- rep(0, p[h])
    Sigma <- outer(1:p[h], 1:p[h], function(i, j) rho^abs(i - j))
    # vector of regression coefficients
    b <- rep(0, p[h])
    
	for(k in 1:length(s)){
        
        A <- 1:s[k]                                 # set of true predictors
        X <- mvrnorm(n, mu = mu, Sigma = Sigma)     # design matrix
        b[A] <- bs                                  # setting the non-zero regression coefficients
        eta <- drop(X[, A] %*% b[A])                # linear predictor
        lambdas <- exp(eta)                         # parameter of the exponential distribution
        
		i <- 0
		repeat{
            
			tm <- rexp(n, lambdas)
			status <- rbinom(n, 1, 1 - pi)
			out.dgcoxph <- try(dgcoxphR(tm, status, X, g0 = g0[k], verbose = FALSE), TRUE)
			out.gof <- try(gof(out.dgcoxph, complexity = "gic"), TRUE)
            
			if(class(out.gof) != "try-error"){
                
				i <- i + 1
				if(i > nsim) break
				bh <- out.dgcoxph$beta
				logplik <- out.dgcoxph$logplik
				df <- out.dgcoxph$df
				gic_compl <- out.dgcoxph$gic_compl
				nd <- sum(out.dgcoxph$status == 1)
				
				######################################
				# AIC section
				aic <- gof(out.dgcoxph, k = "AIC", complexity = "df")
				id <- which.min(aic)
				results[i, h, k, "AIC", "value"] <- aic[id]
				results[i, h, k, "AIC", "size"] <- df[id]
				results[i, h, k, "AIC", "TP"] <- sum(bh[A, id] != 0)
				results[i, h, k, "AIC", "FP"] <- sum(bh[-A, id] != 0)
				results[i, h, k, "AIC", "TN"] <- (p[h] - s[k]) - sum(bh[-A, id] != 0)
				results[i, h, k, "AIC", "FN"] <- s[k] - sum(bh[A, id] != 0)
				results[i, h, k, "AIC", "TPR"] <- sum(bh[A, id] != 0) / s[k]
				results[i, h, k, "AIC", "FPR"] <- sum(bh[-A, id] != 0) / (p[h] - s[k])
				results[i, h, k, "AIC", "FDR"] <- sum(bh[-A, id] != 0) / (sum(bh[A, id] != 0)+sum(bh[-A, id] != 0))
				results[i, h, k, "AIC", "FNR"] <- 1 - (sum(bh[A, id] != 0) / s[k])
				results[i, h, k, "AIC", "F1"] <- 2 * sum(bh[A, id] != 0) / (s[k] + sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				
				######################################
				# BIC section
				bic <- gof(out.dgcoxph,k="BIC",complexity = "df")
				id <- which.min(bic)
				results[i, h, k, "BIC", "value"] <- bic[id]
				results[i, h, k, "BIC", "size"] <- df[id]
				results[i, h, k, "BIC", "TP"] <- sum(bh[A, id] != 0)
				results[i, h, k, "BIC", "FP"] <- sum(bh[-A, id] != 0)
				results[i, h, k, "BIC", "TN"] <- (p[h] - s[k]) - sum(bh[-A, id] != 0)
				results[i, h, k, "BIC", "FN"] <- s[k] - sum(bh[A,id] != 0)
				results[i, h, k, "BIC", "TPR"] <- sum(bh[A, id] != 0) / s[k]
				results[i, h, k, "BIC", "FPR"] <- sum(bh[-A, id] != 0) / (p[h] - s[k])
				results[i, h, k, "BIC", "FDR"] <- sum(bh[-A, id] != 0) / (sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				results[i, h, k, "BIC", "FNR"] <- 1 - (sum(bh[A, id] != 0) / s[k])
				results[i, h, k, "BIC", "F1"] <- 2 * sum(bh[A, id] != 0)/(s[k] + sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				
				######################################
				# Fan_13 section
                		k.gof <- log(log(nd)) * log(p[h])
				Fan_13 <- gof(out.dgcoxph, k = k.gof, complexity = "df")
				id <- which.min(Fan_13)
				results[i, h, k, "Fan_13", "value"] <- Fan_13[id]
				results[i, h, k, "Fan_13", "size"] <- df[id]
				results[i, h, k, "Fan_13", "TP"] <- sum(bh[A, id] != 0)
				results[i, h, k, "Fan_13", "FP"] <- sum(bh[-A, id] != 0)
				results[i, h, k, "Fan_13", "TN"] <- (p[h] - s[k]) - sum(bh[-A, id] != 0)
				results[i, h, k, "Fan_13", "FN"] <- s[k] - sum(bh[A, id] != 0)
				results[i, h, k, "Fan_13", "TPR"] <- sum(bh[A, id] != 0) / s[k]
				results[i, h, k, "Fan_13", "FPR"] <- sum(bh[-A, id] != 0) / (p[h] - s[k])
				results[i, h, k, "Fan_13", "FDR"] <- sum(bh[-A, id] != 0) / (sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				results[i, h, k, "Fan_13", "FNR"] <- 1 - (sum(bh[A,id] != 0) / s[k])
				results[i, h, k, "Fan_13", "F1"] <- 2 * sum(bh[A, id] != 0) / (s[k] + sum(bh[A, id] != 0)+sum(bh[-A, id] != 0))
				
				######################################
				# GIC_aic section
				GIC_aic <- gof(out.dgcoxph, k = "AIC", complexity = "gic")
				id <- which.min(GIC_aic)
				results[i, h, k, "GIC_aic", "value"] <- GIC_aic[id]
				results[i, h, k, "GIC_aic", "size"] <- df[id]
				results[i, h, k, "GIC_aic", "TP"] <- sum(bh[A, id] != 0)
				results[i, h, k, "GIC_aic", "FP"] <- sum(bh[-A, id] != 0)
				results[i, h, k, "GIC_aic", "TN"] <- (p[h] - s[k]) - sum(bh[-A, id] != 0)
				results[i, h, k, "GIC_aic", "FN"] <- s[k] - sum(bh[A, id] != 0)
				results[i, h, k, "GIC_aic", "TPR"] <- sum(bh[A, id] != 0) / s[k]
				results[i, h, k, "GIC_aic", "FPR"] <- sum(bh[-A, id] != 0) / (p[h] - s[k])
				results[i, h, k, "GIC_aic", "FDR"] <- sum(bh[-A, id] != 0)/(sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				results[i, h, k, "GIC_aic", "FNR"] <- 1 - (sum(bh[A, id] != 0) / s[k])
				results[i, h, k, "GIC_aic", "F1"] <- 2 * sum(bh[A, id] != 0) / (s[k] + sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				
				######################################
				# GIC_bic section
				GIC_bic <- gof(out.dgcoxph, k = "BIC", complexity = "gic")
				id <- which.min(GIC_bic)
				results[i, h, k, "GIC_bic", "value"] <- GIC_bic[id]
				results[i, h, k, "GIC_bic", "size"] <- df[id]
				results[i, h, k, "GIC_bic", "TP"] <- sum(bh[A, id] != 0)
				results[i, h, k, "GIC_bic", "FP"] <- sum(bh[-A, id] != 0)
				results[i, h, k, "GIC_bic", "TN"] <- (p[h] - s[k]) - sum(bh[-A, id] != 0)
				results[i, h, k, "GIC_bic", "FN"] <- s[k] - sum(bh[A, id] != 0)
				results[i, h, k, "GIC_bic", "TPR"] <- sum(bh[A, id] != 0) / s[k]
				results[i, h, k, "GIC_bic", "FPR"] <- sum(bh[-A, id] != 0) / (p[h] - s[k])
				results[i, h, k, "GIC_bic", "FDR"] <- sum(bh[-A, id] != 0) / (sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
				results[i, h, k, "GIC_bic", "FNR"] <- 1 - (sum(bh[A, id] != 0) / s[k])
				results[i, h, k, "GIC_bic", "F1"] <- 2 * sum(bh[A, id] != 0) / (s[k] + sum(bh[A, id] != 0) + sum(bh[-A, id] != 0))
								
				#*************************************
				cat("Simulation", i, "with p =", p[h], "and s =", s[k], "completed\n")
				if(i %in% save.point) save.image(file = "gof-n50_tmp.RData")
			} else
                warning("class(out.gof)=='try-error' in i=", i, ", p=" , p[h], ", s=", s[k])
		}
		file.name <- paste("n", n, "_p", p[h], "_s", s[k], ".RData", sep = "")
		save.image(file.name)
	}
}

save.image("gof-n50.RData")

