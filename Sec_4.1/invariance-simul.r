############################################################################################
#
# Author:   Luigi Augugliaro
# e-mail:   luigi.augugliaro@unipa.it
# home:     http://dssm.unipa.it/augugliaro/
#
# Description: R code used for the simulation study reported in "Sparse Relative Risk Regression Models"
# Section 4.1

# Libraris and R code
library(survival)
library(glmnet)
library(glmpath)
library(penalized)
source("dgcoxphR.R")
source("AUC.R")

# Setting
pi <- 0.5 # probability censorship
nsim <- 500 # number of simulatios
ntunpar	<- 100 # number of tuning parameters
n <- 100 # sample size
p <- 10 # number of predictor
s <- 3 # number of true predictors
A <- 1:s # true active set
sc <- seq(1, 4, length = 10) # scale factor
b <- rep(0.5, s)  # regression coefficients

# Output
data.glmnet <- array(NA, dim = c(length(sc), nsim, ntunpar, 2), 
					dimnames = list(scale = 1:length(sc), nsim = 1:nsim, ntunpar = 1:ntunpar, 
					rate = c("TPR", "FPR")))
auc.glmnet <- matrix(NA, length(sc), nsim)

data.coxpath <- array(NA, dim = c(length(sc), nsim, ntunpar, 2), 
					dimnames = list(scale = 1:length(sc), nsim = 1:nsim, ntunpar = 1:ntunpar, 
					rate = c("TPR", "FPR")))
auc.coxpath <- matrix(NA, length(sc), nsim)

data.penalized <- array(NA, dim = c(length(sc), nsim, ntunpar, 2), 
					dimnames = list(scale = 1:length(sc), nsim = 1:nsim, ntunpar = 1:ntunpar, 
					rate = c("TPR", "FPR")))
auc.penalized <- matrix(NA, length(sc), nsim)

data.dgcox <- array(NA, dim = c(length(sc), nsim, ntunpar, 2), 
					dimnames = list(scale = 1:length(sc), nsim = 1:nsim, ntunpar = 1:ntunpar, 
					rate = c("TPR", "FPR")))
auc.dgcox <- matrix(NA, length(sc), nsim)


# starting simulation study

for(h in 1:length(sc)){
	for(i in 1:nsim){
		
		X <- matrix(rnorm(n * p), n, p)
		X <- scale(X, center = TRUE, scale = FALSE)
		w <- sqrt(apply(X, 2, crossprod))
		X <- scale(X, FALSE, w)
		X[, -A] <- X[, -A] * sc[h]
		status  <- rbinom(n, 1, 1 - pi)
		lambdas <- exp(X[, A] %*% b)
		tm <- rexp(n, lambdas)
		y <- Surv(tm, event = status) 
		
        ##########################
		# glmnet section
        ##########################
		out.glmnet <- glmnet(X, y, family = "cox", nlambda = ntunpar, standardize = FALSE)
		b.est <- coef(out.glmnet)
		data.glmnet[h, i, 1:length(out.glmnet$lambda), "TPR"] <- apply(abs(b.est[A, ]) > 0, 2, mean)
		data.glmnet[h, i, 1:length(out.glmnet$lambda), "FPR"] <- apply(abs(b.est[-A, ]) > 0, 2, mean)
		auc.glmnet[h, i] <- AUC(data.glmnet[h, i, 1:length(out.glmnet$lambda), "FPR"], data.glmnet[h, i, 1:length(out.glmnet$lambda), "TPR"])
        
        ##########################
		# glmpath section
        ##########################
		out.coxpath <- coxpath(data = list(x = X, time = tm, status = status), standardize = FALSE)
		b.est <- t(out.coxpath$b.corr)
		data.coxpath[h, i, 1:length(out.coxpath$lambda), "TPR"] <- apply(abs(b.est[A, ]) > 0, 2, mean)
		data.coxpath[h, i, 1:length(out.coxpath$lambda), "FPR"] <- apply(abs(b.est[-A, ]) > 0, 2, mean)
		auc.coxpath[h, i] <-AUC(data.coxpath[h, i, 1:length(out.coxpath$lambda), "FPR"], data.coxpath[h, i, 1:length(out.coxpath$lambda), "TPR"])
		
        ##########################
		#penalized section
        ##########################
		out.penalized <- penalized(Surv(tm, status), penalized = X, lambda1 = 1.0e-4, steps = ntunpar, trace = FALSE, standardize = FALSE)
		b.est <- sapply(out.penalized, coefficients, "p")
		data.penalized[h, i, , "TPR"] <- apply(abs(b.est[A, ]) > 0, 2, mean)
		data.penalized[h, i, , "FPR"] <- apply(abs(b.est[-A, ]) > 0, 2, mean)
		auc.penalized[h, i] <- AUC(data.penalized[h, i, , "FPR"], data.penalized[h, i, , "TPR"])	
		
        ##########################
		# dgcox section
        ##########################
		out.dgcoxphR <- dgcoxphR(tm, status, X, g0 = 0.01, np = ntunpar, verbose = FALSE)
		b.est <- out.dgcoxphR$beta
		data.dgcox[h, i, , "TPR"] <- apply(abs(b.est[A, ]) > 0, 2, mean)
		data.dgcox[h, i, , "FPR"] <- apply(abs(b.est[-A, ]) > 0, 2, mean)
		auc.dgcox[h, i] <- AUC(data.dgcox[h, i, , "FPR"], data.dgcox[h, i, , "TPR"])
		
	}
	cat("simulation with schema", h, "completed!\n")
}

######################################################
# R code used  for the Figure 1 of the main document

pdf("simul_invariance.pdf", width = 5, height = 5)
plot(sc, apply(auc.dgcox, 1, mean), ylim = c(0.19, 0.65), ylab = "Average Area Under the Curve", xlab = "k", cex.lab = 1.3, col = "gray0", lwd = 2, lty = 1, type = "b", pch = 1)
points(sc, apply(auc.glmnet, 1, mean), lwd = 2, lty = 1, col = "gray0", type = "b", pch = 2)
points(sc, apply(auc.coxpath, 1, mean), lwd = 2, lty = 1, col = "gray0", type = "b", pch = 3)
points(sc, apply(auc.penalized, 1, mean), lwd = 2, lty = 1, col = "gray0", type = "b", pch = 4)
legend(x = 2.75, y = 0.5, c("dgCox", "CoxNet", "CoxPath", "CoxPen"), 
	lwd = 2, lty = 1, col = c("gray0", "gray0", "gray0", "gray0"),
	pch = 1:4, box.lwd = 0, box.col = "white", bg = "white", bty = "n", cex = 1.2)
dev.off()







































