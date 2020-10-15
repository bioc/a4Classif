#' Receiver operating curve
#' 
#' A ROC curve plots the fraction of true positives (TPR = true positive rate)
#'  versus the fraction of false positives (FPR = false positive rate) for a binary classifier
#'  when the discrimination threshold is varied. Equivalently, one can also plot
#'  sensitivity versus (1 - specificity).
#' @param object ExpressionSet object for the experiment
#' @param groups String containing the name of the grouping variable. This should be a 
#' the name of a column in the \code{pData} of the \code{expressionSet} object.
#' @param probesetId The probeset ID. These should be stored in the \code{featureNames}
#'  of the \code{expressionSet} object.
#' @param geneSymbol The gene symbol. These should be stored in the column \code{`Gene Symbol`}
#'  in the \code{featureData} of the \code{expressionSet} object.
#' @param main Main title on top of the graph
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#'  (used for the default title of the plot)
#' @param ... Possibility to add extra plot options. See \code{\link{par}}
#' @return a plot is drawn in the current device.
#' prediction object is returned invisibly.
#' @references Some explanation about ROC can be found on 
#' \url{http://en.wikipedia.org/wiki/ROC_curve}
#'  and \url{http://www.anaesthetist.com/mnm/stats/roc/Findex.htm}. 
#' The latter has at the bottom a nice interactive tool 
#' to scroll the cut-off and to see how it affects the FP/TP table and the ROC curve. 
#' @author Willem Talloen
#' @examples 
#' # simulated data set
#' esSim <- simulateData()
#' ROCcurve(probesetId = 'Gene.1', object = esSim, groups = 'type', addLegend = FALSE)
#' @importFrom Biobase exprs featureData featureNames
#' @importFrom ROCR prediction performance plot
#' @export
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
	
	pred <- ROCR::prediction(predictions = exprGene, labels = labels)
	
	### plot ROC curve (x-axis: fpr, y-axis: tpr)
	perf <- ROCR::performance(pred, "tpr", "fpr")
	ROCR::plot(perf, avg= "threshold", colorize = TRUE, main = mainTitle, lwd = 3)
	invisible(pred)
}
