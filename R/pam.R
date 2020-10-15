#' Classify using Prediction Analysis for MicroArrays
#' 
#' Classify using the Prediction Analysis for MicroArrays (PAM) algorithm as implemented 
#' in the pamr package
#' @param object object containing the expression measurements; currently the
#'  only method supported is one for ExpressionSet objects
#' @param groups character string indicating the column containing the class membership
#' @param probe2gene logical; if \code{TRUE} Affymetrix probeset IDs are translated
#' into gene symbols; if \code{FALSE} no such translation is conducted
#' @return object of class \code{pamClass}
#' @references
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and
#'  Gilbert Chu (1999). Diagnosis of multiple cancer types by shrunken
#'  centroids of gene expression.  PNAS 99: 6567-6572.    
#'  Available at \url{www.pnas.org}
#' @references Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'  Microarrays, Chapman \& Hall/CRC, p. 221.
#' @author Willem Talloen
#' @seealso \code{\link[pamr]{pamr.train}}
#' @keywords models
#' @importFrom Biobase featureNames pData featureData
#' @importFrom pamr pamr.train pamr.cv
#' @examples
#' if(require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'  resultPam <- pamClass(object = ALL, groups = "BTtype")
#'  plot(resultPam)
#'  topTable(resultPam, n = 5)
#'  confusionMatrix(resultPam)
#' }
#' @importFrom utils capture.output
#' @export
pamClass <- function(object, groups, probe2gene = TRUE){
  
  labels <- factor(pData(object)[, groups])
  
  gI <- featureNames(object)
  dat <- list(x = as.matrix(exprs(object)), y = labels, geneid = gI)
  
  co <- capture.output(model <- pamr.train(dat))
  co2 <- capture.output(modelCV <- pamr.cv(model, dat, nfold = 10))
  matModelCV <- data.frame(threshold = format(round(modelCV$threshold, 3)),
      nonzero = format(trunc(modelCV$size)), 
      errors = trunc(modelCV$error * nrow(modelCV$yhat)))
  selectSmallError <- which(matModelCV$errors == min(matModelCV$errors))
  selectSmallErrorFewGenes <- selectSmallError[length(selectSmallError)]
  Delta <- as.numeric(as.character(matModelCV[selectSmallErrorFewGenes, "threshold"]))
  res <- list(pamModel = model, pamCV = modelCV, delta = Delta, exprDat = dat)
  res$featureData <- pData(featureData(object))
  res$probe2gene <- probe2gene
  class(res) <- "pamClass"
  return(res)
}  

#' @importFrom pamr pamr.confusion
#' @importFrom a4Core confusionMatrix
#' @export
confusionMatrix.pamClass <- function(x, ...){
  res <- pamr.confusion(x$pamModel, x$delta, extra = FALSE)
  class(res) <- "pamClassConfusionTable"
  return(res)
}

#' @importFrom methods setOldClass
setOldClass("pamClass")

#' Top table for \code{pamClass} object
#' 
#' @importFrom a4Core topTable
#' @importFrom stats var
#' @importFrom pamr pamr.listgenes
#' @inheritParams a4Core::topTable
#' @return \code{topTablePam} object
#' @importFrom methods setMethod
#' @export
setMethod("topTable",
    "pamClass",
    function(fit, n){
      co <- capture.output(listGenes <- pamr.listgenes(fit$pamModel, fit$exprDat, fit$delta,
          fitcv = fit$pamCV, genenames = FALSE))
      if (fit$probe2gene){
        gSymbol <- fit$featureData[listGenes[,'id'], "SYMBOL"]
        listGenes <- data.frame(GeneSymbol = gSymbol, listGenes)
      } else {
        listGenes <- data.frame(listGenes)
      }
      
      rownames(listGenes) <- listGenes[,'id']
      listGenes <- listGenes[, !colnames(listGenes) == 'id']
      
      numberSelGenes <- nrow(listGenes)  
      topList <- listGenes[seq_len(min(n, numberSelGenes)),]
      
      res <- list(topList = topList, numberSelGenes = numberSelGenes, n = n, listGenes = listGenes)
      class(res) <- "topTablePam"
      return(res)
    }
)

#' @export
print.topTablePam <- function(x,  ...){
  cat("Pam selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n" )
  print(x$topList, ...)
}

#' @export
#' @method plot pamClass
#' @importFrom graphics par lines segments axis points
plot.pamClass <- function(x, ...){
      x <- x$pamCV
      n <- nrow(x$yhat)
      y <- rev(x$y)
      if (!is.null(x$newy)) {
        y <- x$newy[x$sample.subset]
      }
      nc <- length(table(y))
      nfolds <- length(x$folds)
      err <- matrix(NA, ncol = ncol(x$yhat), nrow = nfolds)
      temp <- matrix(y, ncol = ncol(x$yhat), nrow = n)
      ni <- rep(NA, nfolds)
      for (i in seq_len(nfolds)) {
        ii <- x$folds[[i]]
        ni[i] <- length(x$folds[[i]])
        err[i, ] <- apply(temp[ii, ] != x$yhat[ii, ], 2, sum)/ni[i]
      }
      se0 <- sqrt(apply(err, 2, var)/nfolds)
      se <- rev(se0)
      MCR <- rev(x$error)

      par(mar=c(5, 4, 4, 2) + 0.1)
      graphics::plot(seq_along(x$size),MCR,xaxt = "n",col="blue",cex=1.2, ylim = c(0.05, (max(MCR)+max(se))),las=1,
          xlab="number of genes",ylab="Misclassification error",main="")
      lines(seq_along(x$size),MCR,col="grey")
      
      segments(seq_along(x$size), MCR - se, seq_len(length(x$size)), MCR + se, col='blue')
      axis(1,seq_along(x$size), rev(x$size), las=2)
      # indicate minimum
      o <- MCR == min(MCR)
      points(min(seq_along(x$size)[o]), MCR[o][1], pch = "x", cex=2)
}

      

