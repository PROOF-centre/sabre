# Function to create a matrix of pairwise module similarities 
# between two lists of modules; the lists should be named
# in a meaningful way.
distm <- function(list1, list2, method="simpson") {
	dtable = matrix(NA,nrow=length(list1),ncol=length(list2))	
	for (i in 1:length(list1))
	for (j in 1:length(list2)) {
		dtable[i,j] = similarity(list1[[i]],list2[[j]],method)
	}
	rownames(dtable) = names(list1)
	colnames(dtable) = names(list2)
	
	return(dtable)
}
