#############################################################
## Description: R code of the simulation study reported in
## section 4.2 of the manuscript entitled
## 'Sparse Relative Risk Regression Models'
##
## Authors = H. Pazira, J. Gonzalez and F. Abegaz
##
## Setting
## s = 5, n = 50, rho = 0.6, pi = 0.2
#############################################################

# Libraries and R codes
library(MASS)
library(survival)
library(plyr)
library(glmnet)
library(glmpath)
library(penalized)
source("dgcoxphR.R")
source("ROC.R")

########################
# general setting
#
filename <- "rho06_pi02.RData"
nsim <- 100                             # number of simulations
n <- 50                                 # sanple size
p <- 100                                # number of variables
s <- 5                                  # number of non-zero coefficients
b <- c(rep(0.2, s), rep(0, p - s))      # subset of the non-zero regression coefficients
rho <- 0.6                              # parameter of the autoregressive convariance matrix
pi <- 0.2                               # probability of censorship
n.steps <- 100                          # number of points of the solution curve

###########################################
## matrices to save the results
FP.dgCox <- TP.dgCox <- matrix()
FP.penal <- TP.penal <- FP.coxpath <- TP.coxpath <- FP.coxnet <- TP.coxnet <- matrix(0, nsim, n.steps)

###########################################
## starting simulation study

mu <- rep(0, p) # expected value
R <- outer(1:p, 1:p, function(i,j) rho^abs(i-j)) # covariance matrix

for(k in 1:nsim){
    
    X <- mvrnorm(n, mu, R)
    ## generate the survival times using the exponential distribution
    lambdas <- exp(drop(X[, 1:s] %*% b[1:s]))
    tm <- rexp(n, lambdas)
    ## generate the censorship, uniformly across the data
    status <- rbinom(n, 1, 1-pi)

    ##################
    ## Methods
    ##

    ## glmnet
    out.glmnet <- glmnet(X, cbind(time = tm, status = status), family = "cox", nlambda = n.steps)
    b.coxnet <- t(out.glmnet$beta)

    # dgCox
    maxdev_ratio <- max(out.glmnet$dev.ratio)
    out.dgcoxphR <- dgcoxphR(tm, status, X, g0 = 0.2, maxrsq_ratio = maxdev_ratio, nstep = n.steps, verbose = FALSE)
    b.dgCox <- t(out.dgcoxphR$b)
        
    ##  Cox path
    out.coxpath <- coxpath(list(x = X, time = tm, status = status), max.steps = n.steps , method = "efron")
    b.coxpath <- out.coxpath$b.corrector
             
    # Cox penalized
    out.penal <- penalized(Surv(tm, status), penalized = X, lambda1 = 1, steps = n.steps)
    b.penal <- t(sapply(out.penal, coefficients, "p"))
     

    ####################
    ## ROC curves
    ##
    ROC.coxnet <- ROC(b, b.coxnet)
    ROC.dgCox <- ROC(b, b.dgCox)
    ROC.coxpath <- ROC(b, b.coxpath)
    ROC.penal <- ROC(b, b.penal)
    
    ## computing TRUE and FALSE POSITIVE RATE
    TP.dgCox <- rbind.fill.matrix(TP.dgCox, t(ROC.dgCox$TP))
    FP.dgCox <- rbind.fill.matrix(FP.dgCox, t(ROC.dgCox$FP))
    TP.coxnet[k, 1:length(ROC.coxnet$TP)] <- ROC.coxnet$TP
    FP.coxnet[k, 1:length(ROC.coxnet$FP)] <- ROC.coxnet$FP
    TP.penal[k, 1:length(ROC.penal$TP)] <- ROC.penal$TP
    FP.penal[k, 1:length(ROC.penal$FP)] <- ROC.penal$FP
    TP.coxpath[k, 1:length(ROC.coxpath$TP)] <- ROC.coxpath$TP
    FP.coxpath[k, 1:length(ROC.coxpath$FP)] <- ROC.coxpath$FP

    # print iteration
    print(k)
}

###################################
## average ROC curves + AUC

tpr <- apply(TP.dgCox, 2, mean, na.rm = TRUE)
fpr <- apply(FP.dgCox, 2, mean, na.rm = TRUE)
last_negative <- as.numeric(which(round(tpr[-1], 4) - round(tpr[-length(tpr)], 4) < 0)[1])
if(is.na(last_negative)) last_negative <- length(tpr)
new_fpr <- ifelse(is.na(as.numeric(fpr[last_negative])), max(fpr), as.numeric(fpr[last_negative]))

## AUC
tpr_coxnet <- apply(TP.coxnet, 2, mean, na.rm = TRUE)
fpr_coxnet <- apply(FP.coxnet, 2, mean, na.rm = TRUE)
tpr_coxpath <- apply(TP.coxpath, 2, mean, na.rm = TRUE)
fpr_coxpath <- apply(FP.coxpath, 2, mean, na.rm = TRUE)
tpr_penal <- apply(TP.penal, 2, mean, na.rm = TRUE)
fpr_penal <- apply(FP.penal, 2, mean, na.rm = TRUE)

new_max_mean_fpr_dgcox <- as.numeric(which(fpr == new_fpr))[1]
new_max_mean_fpr_net <- as.numeric(which(fpr_coxnet >= new_fpr))[1]
new_max_mean_fpr_path <- as.numeric(which(fpr_coxpath >= new_fpr))[1]
new_max_mean_fpr_penal <- as.numeric(which(fpr_penal >= new_fpr))[1]

ROC.dgCox_TP <- tpr[1:new_max_mean_fpr_dgcox]
ROC.dgCox_FP <- fpr[1:new_max_mean_fpr_dgcox]
ROC.coxnet_TP <- tpr_coxnet[1:new_max_mean_fpr_net]
ROC.coxnet_FP <- fpr_coxnet[1:new_max_mean_fpr_net]
ROC.coxpath_TP <- tpr_coxpath[1:new_max_mean_fpr_path]
ROC.coxpath_FP <- fpr_coxpath[1:new_max_mean_fpr_path]
ROC.penal_TP <- tpr_penal[1:new_max_mean_fpr_penal]
ROC.penal_FP <- fpr_penal[1:new_max_mean_fpr_penal]

AUC_dgCox <- sum((ROC.dgCox_FP[-1]-ROC.dgCox_FP[-length(ROC.dgCox_FP)])*(ROC.dgCox_TP[-1]+ROC.dgCox_TP[-length(ROC.dgCox_TP)])/2) + sum((1-ROC.dgCox_FP[length(ROC.dgCox_FP)])*(1+ROC.dgCox_TP[length(ROC.dgCox_TP)])/2)
AUC_coxnet <- sum((ROC.coxnet_FP[-1]-ROC.coxnet_FP[-length(ROC.coxnet_FP)])*(ROC.coxnet_TP[-1]+ROC.coxnet_TP[-length(ROC.coxnet_TP)])/2) + sum((1-ROC.coxnet_FP[length(ROC.coxnet_FP)])*(1+ROC.coxnet_TP[length(ROC.coxnet_TP)])/2)
AUC_coxpath <- sum((ROC.coxpath_FP[-1]-ROC.coxpath_FP[-length(ROC.coxpath_FP)])*(ROC.coxpath_TP[-1]+ROC.coxpath_TP[-length(ROC.coxpath_TP)])/2) + sum((1-ROC.coxpath_FP[length(ROC.coxpath_FP)])*(1+ROC.coxpath_TP[length(ROC.coxpath_TP)])/2)
AUC_penal <- sum((ROC.penal_FP[-1]-ROC.penal_FP[-length(ROC.penal_FP)])*(ROC.penal_TP[-1]+ROC.penal_TP[-length(ROC.penal_TP)])/2) + sum((1-ROC.penal_FP[length(ROC.penal_FP)])*(1+ROC.penal_TP[length(ROC.penal_TP)])/2)

round(c(AUC_dgCox, AUC_coxnet, AUC_coxpath, AUC_penal), 3)

###################################
## Figure + AUC

pdf(paste("fig_n50_rho06_pi02.pdf" ),width=5,height=5)

plot(fpr, tpr, type = "l", lwd=2, xlab = "False positive rate", ylim = c(0,1.2), xlim = c(0,new_fpr),
ylab = "True positive rate", main = expression(rho==0.6 ~~ and ~~ pi==0.2), cex.main = 1.75,
cex.lab = 1.3, col = "gray0", axes = FALSE)
axis(1)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), las = 1)


lines(apply(FP.coxnet,2,function(x)mean(na.omit(x))), apply(TP.coxnet,2,function(x)mean(na.omit(x))),
lwd = 2, lty = 2, col = "gray25")

lines(apply(FP.coxpath,2,function(x)mean(na.omit(x))), apply(TP.coxpath,2,function(x)mean(na.omit(x))),
lwd = 2, lty = 3, col = "gray50")

lines(apply(FP.penal,2,function(x)mean(na.omit(x))), apply(TP.penal,2,function(x)mean(na.omit(x))),
lwd = 2, lty = 4, col = "gray75")

legend(x = -0.005, y = 1.22, c("dgCox", "CoxNet", "CoxPath", "CoxPen"),
lwd = 2, lty = 1:4, col = c("gray0", "gray25", "gray50", "gray75"),
box.lwd = 0, box.col = "white", bg = "white", bty = "n", cex = 1.1)

legend(x = 0.10, y = 1.22, c("(AUC = 0.824)", "(AUC = 0.778)", "(AUC = 0.775)", "(AUC = 0.763)"),
box.lwd = 0, box.col = "white", bg = "white", bty = "n", cex = 1.1)

abline(a = 0, b = 1, lwd = 2, lty = 1, col = "gray")

dev.off()

###################################
## save results
save.image(file = filename)




