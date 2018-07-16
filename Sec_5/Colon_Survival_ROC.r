###############################################################
#
# Author:   Hassan Pazira
# e-mail:   h.pazira@rug.nl
# 
# Description: R code used for the application study reported in 
# "Sparse Relative Risk Regression Models"
# Section 5 and Supplementary Materials

# Libraris and R code
library(survival)
library(survivalROC)
library(glmnet)
library(glmpath)
library(penalized)
library(dglars)
source('dgcoxphR.R')


###############################################################
### Read the data: Colon
data = read.table("GSE28722_series_matrix.txt",row.names = 1, sep = "\t")
data = data.frame(t(data))
attach(data)
# select not NAs
sel = !is.na(time) & !is.na(status)
# define elemnets
X = as.matrix(data[sel,3:dim(data)[2]])
X = X[,colMeans(is.na(X)) == 0] 
time = as.numeric(data[sel,1])
censor = as.numeric(data[sel,2])
###############################################################


## Preselection of more expressed
sdevs = apply(X,2,sd)
gene.sel = sdevs > 0
p                  <- ncol(X[,gene.sel]);p
n                  <- nrow(X);n

# Setting
set.seed(11235)
id <- sample(2, n, replace=TRUE, prob=c(1.0, 0.0)) # c(1.0, 0.0) # c(0.6, 0.4)
train <- id == 1 ; sum(train)
n.train = sum(train)
test <- !train
n.test = sum(test)
#perc = 0.6
#train = sample(1:n,round(perc*n))
#test  = (1:n)[-train]
#n.train            <- length(train);n.train
#n.test             <- length(test);n.test
time_test = time[test]
censor_test = censor[test]
time_train = time[train]
censor_train = censor[train]


#############################
## dglars train data
##

## run procedure
out.dgcoxph.train  <- dgcoxphR(time_train, censor_train, X[train,],g0=1.7,verbose=T)
out.gof.train      <- gof(out.dgcoxph.train,complexity = "gic")
bh.train           <- out.dgcoxph.train$beta
logplik.train      <- out.dgcoxph.train$logplik
df.train           <- out.dgcoxph.train$df
gic_compl.train    <- out.dgcoxph.train$gic_compl

## Calculate Model Selection Criteria
aic.train                  <- -2 * logplik.train + 2 * df.train # Selected for the paper
bic.train                  <- -2 * logplik.train + log(n.train) * df.train # == gof(out.dgcoxph.train,complexity = "df")
GIC_bic.train              <- -2 * logplik.train + log(n.train) * gic_compl.train # == gof(out.dgcoxph.train,complexity = "gic")
GIC_aic.train              <- -2 * logplik.train + 2 * gic_compl.train # == gof(out.dgcoxph.train,complexity = "gic")


#### Separate by Weigth

## Train
W_train = X[train,gene.sel]%*%as.matrix(out.dgcoxph.train$beta[,which.min(GIC_aic.train)])

## Test
W_test = X[test,gene.sel]%*%as.matrix(out.dgcoxph.train$beta[,which.min(GIC_aic.train)])


####################
## Calculate inputs

tm<-12

# KM
ovarian.ROC<-survivalROC(time_train, censor_train, W_train, predict.time = tm, method="KM")
ovarian.ROC$AUC

# NNE & lambda=0.05
ovarian.ROC.NNE<-survivalROC(time_train, censor_train, W_train, predict.time = tm, lambda=0.05)
ovarian.ROC.NNE$AUC


##############################################
# OTHER METHODS

# FIT THE MODEL WITH ALL DATA
# SELECT OPTIMAL PENALTY
# CALCULATE THE "SCORE" (W_penalized) BASED ON SELECTED BETA
# APPLY survivalROC
# PLOT THE SURVIVAL ROC (together with the other methods)


# GLMNET
y<-cbind(time_train, censor_train)
colnames(y)<-c("time","status")
out.glmnet.train  <- glmnet(x= X[train,], y=y, family="cox")
out.cv.net<-cv.glmnet(x= X[train,], y=y, family="cox" )

## separate by weigth
W_glmnet = X[train,gene.sel]%*%as.matrix(out.glmnet.train$beta[,which(out.cv.net$lambda==out.cv.net$lambda.min)])

ovarian.ROC.glmnet<-survivalROC(time_train, censor_train, W_glmnet, predict.time = tm, lambda=0.05)
ovarian.ROC.glmnet$AUC

# GLMPATH
out.glmpath.train  <- coxpath(list(x=X[train,], time=time_train, status=(censor_train)), method="efron")
out.cv.path <- cv.coxpath(list(x=X[train,], time=time_train, status=(censor_train)), method="efron",plot.it = F,se=F)

## separate by weigth
W_glmpath = X[train,gene.sel]%*%as.matrix(t(out.glmpath.train$b.corrector)[,which.min(out.cv.path$cv.error[!is.infinite(out.cv.path$cv.error)])]) # out.cv.path$cv.se

ovarian.ROC.glmpath<-survivalROC(time_train, censor_train, W_glmpath, predict.time = tm, lambda=0.05)
ovarian.ROC.glmpath$AUC

# PENALIZED
out.pen.train  <- penalized(Surv(time_train, censor_train), penalized = X[train,], lambda1 = 1, model ="cox" )
#out.cv.pen<-cvl(Surv(time_train, censor_train), penalized = X[train,], lambda1 = 1, model ="cox" )
out.opt.pen.train<-optL1(Surv(time_train, censor_train), penalized = X[train,], model ="cox" )

## separate by weigth
#W_pen = X[train,gene.sel]%*%as.matrix(coefficients(out.pen.train, "penalized"))
W_pen = X[train,gene.sel]%*%as.matrix(coefficients(out.opt.pen.train$fullfit, "penalized"))

ovarian.ROC.pen<-survivalROC(time_train, censor_train, W_pen, predict.time = tm, lambda=0.05)
ovarian.ROC.pen$AUC


#### Final Plot
pdf("Colon-survival-roc.pdf",width=6,height=6)
plot(ovarian.ROC.NNE$FP,ovarian.ROC.NNE$TP,xlab="False positive", ylab="True positive",type="l")
ovarian.ROC.NNE$AUC
lines(ovarian.ROC.glmnet$FP,ovarian.ROC.glmnet$TP,lty=2)
ovarian.ROC.glmnet$AUC
lines(ovarian.ROC.glmpath$FP,ovarian.ROC.glmpath$TP,lty=3)
ovarian.ROC.glmpath$AUC
lines(ovarian.ROC.pen$FP,ovarian.ROC.pen$TP,lty=4)
ovarian.ROC.pen$AUC
title("Colon survival ROC")
legend("bottomright",c(paste( "dgCox ( AUC=",round(ovarian.ROC.NNE$AUC,3),")"), paste( "glmnet ( AUC=",round(ovarian.ROC.glmnet$AUC,3),")"), paste( "glmpath ( AUC=",round(ovarian.ROC.glmpath$AUC,3),")"), paste( "penalized ( AUC=",round(ovarian.ROC.glmpath$AUC,3),")")),lty=c(1,2,3,4))
dev.off()

save.image("Colon.RData")



#############################
## plots survival curves
##
#surv.pos = survfit(Surv(time_test[W_test>0], censor_test[W_test>0]) ~ 1)
#surv.neg = survfit(Surv(time_test[W_test<0], censor_test[W_test<0]) ~ 1)
#surv     = survfit(Surv(time_train, censor_train) ~ 1)

#pdf("ovarian.pdf",width=7,height=7)
#plot(surv.pos,main = "Ovarian Cancer",xlab= "Survival times (months)", ylab ="Survival prob.",col="red",lwd=3,cex.main=2,cex.lab=1.5,cex.axis=1.5)
#lines(surv.neg,col="blue",lwd=3,lty=2)
#lines(surv,col="grey")
#legend("topright", legend = c("Low surv., test", "High surv., test", "Surv. training"),lty = c(1,2,1),lwd=c(2,2,1),col=c("red","blue","grey"),cex=1.2)
#dev.off()

#survtest1 <- survdiff(Surv(time_test, censor_test) ~ sign(W_test),rho=0) 
#survtest2 <- survdiff(Surv(time_test, censor_test) ~ sign(W_test),rho=1) 

# number of active genes
#sum(sign(abs(as.matrix(out.dgcoxph.train$beta[,which.min(GIC_aic.train)]))))



