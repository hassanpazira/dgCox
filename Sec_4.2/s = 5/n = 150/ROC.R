ROC <- function(B, B.est){
    notB <- abs(sign(abs(B)) - 1)
    B.true <- as.matrix(rep(1, nrow(B.est))) %*% B
    B.false <- as.matrix(rep(1, nrow(B.est))) %*% notB
    TP <- apply(abs(sign(B.est * B.true)), 1, sum) / sum(sign(abs(B)))
    FP <- apply(abs(sign(B.est * B.false)), 1, sum)/ sum(sign(abs(notB)))
    list(TP = TP, FP = FP)
}
