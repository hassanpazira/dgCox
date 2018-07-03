################################################################################################################
# Data:	8/10/2013
# Author: Luigi Augugliaro
#
# Note: Internal function not intended for the user

setDiff <- function(new, old){
	if(setequal(new, old)){
		D <- ""
	} else {
		new_old <- union(new, old)
		plus <- setdiff(new_old, old)
		if(length(plus) ! =0) plus <- paste("  +", paste(plus, collapse = " +"), sep = "")
		minus <- setdiff(new_old,new)
		if(length(minus) != 0) minus <- paste("  -", paste(minus, collapse = " -"), sep = "")
		D <- paste(plus, minus, sep = "")
	}
	D
}
