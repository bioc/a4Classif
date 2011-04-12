# 
# PAM using package 'pamr'
#
###############################################################################

pamClass <- function(object, groups, probe2gene = TRUE){
  labels <- factor(pData(object)[, groups])
  
  gI <- featureNames(object)
  sI <- labels
  dat <- list(x = as.matrix(exprs(object)), y = labels, geneid = gI)
  
  co <- capture.output(model <- pamr.train(dat))
  co2 <- capture.output(modelCV <- pamr.cv(model, dat, nfold = 10))
  matModelCV <- data.frame(threshold = format(round(modelCV$threshold, 3)),
      nonzero = format(trunc(modelCV$size)), 
      errors = trunc(modelCV$error * nrow(modelCV$yhat)))
  selectSmallError <- which(matModelCV$errors == min(matModelCV$errors))
  selectSmallErrorFewGenes <- selectSmallError[length(selectSmallError)]
  Delta <- as.numeric(as.character(matModelCV[selectSmallErrorFewGenes,'threshold']))
  res <- list(pamModel = model, pamCV = modelCV, delta = Delta, exprDat = dat)
  res$featureData <- pData(featureData(object))
  res$probe2gene <- probe2gene
  class(res) <- "pamClass"
  return(res)
}  

confusionMatrix.pamClass <- function(x, ...){
  res <- pamr.confusion(x$pamModel, x$delta)
  class(res) <- "pamClassConfusionTable"
  return(res)
}

setOldClass("pamClass")

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
      topList <- listGenes[1:(min(n, numberSelGenes)),]
      
      res <- list(topList = topList, numberSelGenes = numberSelGenes, n = n, listGenes = listGenes)
      class(res) <- "topTablePam"
      return(res)
    }
)

print.topTablePam <- function(x,  ...){
  cat("Pam selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n" )
  print(x$topList, ...)
}


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
      for (i in 1:nfolds) {
        ii <- x$folds[[i]]
        ni[i] <- length(x$folds[[i]])
        err[i, ] <- apply(temp[ii, ] != x$yhat[ii, ], 2, sum)/ni[i]
      }
      se0 <- sqrt(apply(err, 2, var)/nfolds)
      se <- rev(se0)
      MCR <- rev(x$error)

      par(mar=c(5, 4, 4, 2) + 0.1)
      plot(seq(1:length(x$size)),MCR,xaxt = "n",col="blue",cex=1.2, ylim = c(0.05, (max(MCR)+max(se))),las=1,
          xlab="number of genes",ylab="Misclassification error",main="")
      lines(seq(1:length(x$size)),MCR,col="grey")
      
      segments(seq(1:length(x$size)), MCR - se, seq(1:length(x$size)), MCR + se, col='blue')
      axis(1,seq(1:length(x$size)), rev(x$size), las=2)
      # indicate minimum
      o <- MCR == min(MCR)
      points(min(seq(1:length(x$size))[o]), MCR[o][1], pch = "x", cex=2)
}

      

