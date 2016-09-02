
library(WGCNA)
enableWGCNAThreads()

source("calcHScore.r")
source("distm.r")
source("runWGCNA.r")
source("similarity.r")
data.all = read.table("sample gene expression.csv",sep=",",head=TRUE,row.names=1)

samples.all = colnames(data.all)

#################################################
############### Parameters to set ###############
#################################################

n = 10				# number of bootstrap iterations
seed = 273462987	# seed value, for reproducibility

#################################################
#################################################

set.seed(seed)

dir.create("modules",showWarnings=FALSE)

# Run WGCNA on the whole data set and bootstrap resamplings
for (i in 0:n) { 
	if (i!=0) samples = samples.all
	if (i==0) samples = sample(samples.all,replace=TRUE)

	data = t(data.all[,samples])
	rownames(data) = samples.all	# Rename samples so there are no duplicate rownames

	moduleList = runWGCNA(data,minModuleSize=50)

	if (i!=0) save(moduleList, file=paste0("modules/modules_run",i,".Rdata"))
	if (i==0) save(moduleList, file=paste0("modules/modules.Rdata"))
}

############################################

load("modules/modules.Rdata")
refmodules = moduleList

# This matrix will store the maximum (i.e. "best match") 
# Simpson index for each reference module, when compared
# to all the modules generated in each bootstrap iteration
bestMatch = matrix(NA,nrow=n,ncol=length(refmodules))
colnames(bestMatch) = names(refmodules)

for (i in 1:n) {
    load(paste("modules/modules_run",i,".Rdata",sep=""))
	cv = moduleList
	
	# Matrix of pairwise Simpson indices
	psim = data.matrix(distm(refmodules, cv, "simpson"))
	# Find the maximum value for each reference module
	bestMatch[i,rownames(psim)] = apply(psim,1,max)
}

############################################
# H-score calculation

# Granularity of the h-score depends on the
# number of bootstrap iterations run
cutoffs = seq(0,1,by=1/n)

hscore = calcHScore(bestMatch, cutoffs)
write.table(hscore,file="module h-scores.csv",sep=",",col.names=FALSE)


