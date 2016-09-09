# Function to calculate similarity between two sets of modules
similarity <- function(m1, m2, method="simpson") {
	if (method=="jaccard") return(length(intersect(m1,m2))/length(union(m1,m2)))
	if (method=="simpson") return(length(intersect(m1,m2))/min(length(m1),length(m2)))
}
