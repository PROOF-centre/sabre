# Wrapper function to run WGCNA on a data matrix (rows correspond to samples, columns
# to genes) and process the output into a list containing vectors of genes belonging
# to the same module
runWGCNA <- function(datExpr, minConnectivity=0.85, TOMType="unsigned", blockSize=5000,
	powers=c(seq(4,10,by=1),seq(12,20,by=2)), minModuleSize=50, deepSplit=2, mergeCutHeight=0.15,
	minKMEtoStay=0.3) {
	
	# Call the network topology analysis function
	powerTables = pickSoftThreshold(datExpr, powerVector=powers, blockSize=blockSize, verbose = 2)[[2]]
	collectGarbage();

	# Pick the smallest power that gives us the required connectivity,
	# but if no such power exists, picked the power that gives us
	# the largest connectivity in the table
	sft = -sign(powerTables[,3])*powerTables[,2]
	pwr = powerTables[which(sft>=minConnectivity)[1],1]
	if (is.na(pwr)) pwr = powerTables[which.max(sft)[1],1]

	# Network construction
	net = blockwiseModules(datExpr,	TOMType=TOMType, power=pwr,
		maxblockSize=blockSize,	minModuleSize=minModuleSize, deepSplit=deepSplit,
		pamRespectsDendro=FALSE, mergeCutHeight=mergeCutHeight, numericLabels=TRUE,
		minKMEtoStay=minKMEtoStay, saveTOMs=FALSE, verbose=5)

	moduleLabels = net$colors;
	# Convert the numeric labels to color labels
	moduleColors = labels2colors(moduleLabels)
	consTree = net$dendrograms[[1]];
	index = order(consTree$height)
	consTree$height = consTree$height[index]

	# Calculate module eigenvalues with color labels
	MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)
	modNames = substring(names(MEs), 3)
	
	# Calculate gene module membership as correlation with module eigenvalues
	geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)));
	names(geneModuleMembership) = paste("MM", modNames, sep="");
	geneModuleMembership = data.frame(moduleColors=moduleColors,geneModuleMembership)
	
	# The grey module is removed as it is the trash module output by WGCNA
	geneModuleMembership = geneModuleMembership[geneModuleMembership$moduleColors!="grey",]
	
	mnames = levels(geneModuleMembership$moduleColors)
	moduleList = lapply(mnames,function(m) return(rownames(geneModuleMembership)[geneModuleMembership$moduleColors==m]))
	names(moduleList) = mnames
	
	return(moduleList)
}

