###############################################
# Function to compute the Area Under the Curve

AUC <- function(FPR, TPR){
	id <- order(FPR)
	FPR <- FPR[id]
	TPR <- TPR[id]
	TPR <- tapply(TPR, FPR, max)
	FPR <- unique(FPR)
	dFPR <- c(diff(FPR), 0)
	dTPR <- c(diff(TPR), 0)
	out <- sum(TPR * dFPR) + sum(dTPR * dFPR) / 2
	out
}
