################################################################################################################
# Data:	8/10/2013
# Author: Luigi Augugliaro
#
# Description: print method for an object with class "dgcoxphR"
#
# Arguments
# x: the output of the dgcoxphR function

print.dgcoxphR <- function (x, digits = max(3, getOption("digits") - 3), ...){
	beta <- x$beta
	np <- ncol(beta)
	vn <- rownames(beta)
	id <- abs(beta) > 0
	mdl.list <- apply(id, 2, function(z) vn[z])
	action <- rep("", length = np)
	for(j in 2:np) action[j - 1] <- setDiff(mdl.list[[j]], mdl.list[[j - 1]])
	g <- x$g
	dev <- x$dev
	dev_ratio <- x$dev_ratio
	rsq <- x$rsq
	maxrsq <- x$maxrsq
	df <- x$df
	tbl <- data.frame(action, g, dev, dev_ratio, rsq, df)
	names(tbl) <- c("Sequence", "g", "Dev", "%Dev", "pseudo-R2", "df")
	id <- which(action != "")
	n.tbl <- dim(tbl)[1]
	n.space <- length(id)
	id.tbl <- vector(length = n.tbl + n.space, mode = "numeric")
	id.space <- id + seq(1, n.space)
	id.tbl[-id.space] <- seq(1:n.tbl)
	id.tbl[id.space] <- id
	tbl.format <- format(tbl[id.tbl,], digits = digits)
	tbl.format[id.space - 1, 1] <- ""
	tbl.format[id.space, -1] <- ""
	print.data.frame(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	invisible(tbl)
}
