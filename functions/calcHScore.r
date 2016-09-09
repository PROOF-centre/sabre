# Function to calculate h-score given a matrix of "best match"
# similarities for a reference module set across a number of
# bootstrap iterations (or other perturbations of the data)
calcHScore <- function(bestMatch, cutoffs) {
	hscore = sapply(colnames(bestMatch), function(m) {
		passThres = sapply(cutoffs, function(k) {
				return(mean(bestMatch[,m]>=k)>=k)
			})
		return(max(cutoffs[passThres]))
	})
	return(hscore)
}

