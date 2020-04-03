ROCcurve <- function (object, groups, probesetId = NULL, 
		geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...){
	
	if (!is.null(probesetId) & !is.null(geneSymbol))
		stop("Please provide either a 'probeset' or a 'gene'")
	
	### create gene expression vector
	if (is.null(geneSymbol)){ # probeset given
		probesetId <- as.character(probesetId)     # use names not position !!
		exprGene <- exprs(object)[probesetId, ]
	} else { # gene given
		probesetPos <- which(geneSymbol == featureData(object)$SYMBOL)
		if (!length(probesetPos))
			stop("gene 'gene' does not occur in ExpressionSet 'object'")
		
		probesetId <- featureNames(object)[probesetPos]
		if (length(probesetId) > 1)
			warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
							"probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
		exprGene <- exprs(object)[probesetId[1], ] # use names not position !!
	}
	
	labels <- factor(pData(object)[, groups])
	# sort levels of the labels so that the group with on average the lowest values comes first
	rankMeanLabels <- rank(by(exprGene, labels, mean))
	levels(labels) <- levels(labels)[rankMeanLabels]
	
	# prepare title
	if (probe2gene){
		gSymbol <- featureData(object)[probesetId[1],]$SYMBOL
	}
	
	mainTitle <- if (is.null(main)){
				if (probe2gene)
					paste(gSymbol, " (", probesetId[1], ")", sep = "")
				else probesetId[1]
			} else { 
				main 
			}
	
	pred <- prediction(predictions = exprGene, labels = labels)
	
	### plot ROC curve (x-axis: fpr, y-axis: tpr)
	perf <- performance(pred, "tpr", "fpr")
	ROCR::plot(perf, avg= "threshold", colorize = TRUE, main = mainTitle, lwd = 3)
	invisible(pred)
}
