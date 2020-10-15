##############################################
### Random Forest using package 'varSelRF' ###
##############################################

#' Classify using Random Forests
#' 
#' Classify using the Random Forest algorithm of Breiman (2001)
#' @param object object containing the expression measurements; currently the
#'  only method supported is one for ExpressionSet objects
#' @param groups character string indicating the column containing the class 
#'  membership
#' @param probe2gene logical; if \code{TRUE} Affymetrix probeset IDs are translated
#'  into gene symbols in the output object; if \code{FALSE} no such translation 
#'  is conducted
#' @return Object of class 'rfClass'
#' @references Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
#'  5-32.
#' @note topTable and plot methods are available for 'rfClass' objects.
#' @author Tobias Verbeke and Willem Talloen
#' @seealso \code{\link[randomForest]{randomForest}}
#' @keywords models
#' @importFrom Biobase pData exprs featureData featureNames
#' @importFrom varSelRF varSelRF
#' @examples 
#' if(require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'  # select only a subset of the data for computation time reason
#'  ALLSubset <- ALL[sample.int(n = nrow(ALL), size = 100, replace = TRUE), ]
#'  resultRf <- rfClass(object = ALLSubset, groups = "BTtype")
#'  plot(resultRf)
#'  topTable(resultRf, n = 15)
#' }
#' @export
rfClass <- function(object, groups, probe2gene = TRUE){
  labels <- factor(pData(object)[, groups])
  fit <- varSelRF(t(exprs(object)), labels, ntree = 200,
        ntreeIterat = 100,  vars.drop.frac = 0.2, verbose = FALSE)
  
  # transfer annotation
  if (probe2gene){
    fit$gSymbol <- featureData(object)$`SYMBOL`
    names(fit$gSymbol) <- featureNames(object)
  }

  class(fit) <- c("rfClass")# , class(fit))
  return(fit)
}

#' @importFrom methods setOldClass
setOldClass("rfClass")

#' Top table for \code{rfClass} object
#' @importFrom a4Core topTable
#' @inheritParams a4Core::topTable
#' @return \code{topTableRfClass} object
#' @importFrom methods setMethod
#' @export
setMethod("topTable", 
    "rfClass",
    function(fit, n){
      selGenes <- fit$selected.vars
      numberSelGenes <- length(selGenes)
      topProbes <- selGenes[seq_len(min(n, numberSelGenes))]
      topList <- data.frame(GeneSymbol=fit$gSymbol[topProbes])
      row.names(topList) <- topProbes 
      res <- list(topList = topList, numberSelGenes = numberSelGenes, n = n)
      class(res) <- "topTableRfClass"
      return(res)      
    }
)

#' @export
print.topTableRfClass <- function(x,  ...){
  cat("Random forest selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}

#' @export
#' @method plot rfClass
#' @importFrom graphics par lines segments axis points
plot.rfClass <- function(x, ...){
      size <- x$selec.history$Number.Variables
      MCR0 <- x$selec.history$OOB
      se0 <- x$selec.history$sd.OOB
      
      se <- rev(se0)
      MCR <- rev(MCR0)

      #create plot
      par(mar = c(5, 4, 4, 2) + 0.1)
	  graphics::plot(seq_along(size), MCR, xaxt = "n", col = "blue", cex = 1.2, ylim = c(0.05, (max(MCR)+max(se))),las=1,
          xlab = "number of genes", ylab = "Misclassification error", main = "")
      lines(seq_along(size), MCR, col = "grey")
      
      segments(seq_along(size), MCR - se, seq_along(size), MCR + se, col='blue')
      axis(1, seq_along(size), rev(size), las=2)
      # indicate chosen gene set size
      o <- rev(size) == length(x$selected.vars)
                    # alternatve:  # indicate minimum
                    #              o <- MCR == min(MCR)
      points(min(seq_along(size)[o]), MCR[o][1], pch = "x", cex = 2)
}
