##############################################
### Random Forest using package 'varSelRF' ###
##############################################

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

setOldClass("rfClass")

setMethod("topTable", 
    "rfClass",
    function(fit, n){
      selGenes <- fit$selected.vars
      numberSelGenes <- length(selGenes)
      topProbes <- selGenes[1:min(n, numberSelGenes)]
      topList <- data.frame(GeneSymbol=fit$gSymbol[topProbes])
      row.names(topList) <- topProbes 
      res <- list(topList = topList, numberSelGenes = numberSelGenes, n = n)
      class(res) <- "topTableRfClass"
      return(res)      
    }
)

print.topTableRfClass <- function(x,  ...){
  cat("Random forest selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}

plot.rfClass <- function(x, ...){
      size <- x$selec.history$Number.Variables
      MCR0 <- x$selec.history$OOB
      se0 <- x$selec.history$sd.OOB
      
      se <- rev(se0)
      MCR <- rev(MCR0)

      #create plot
      par(mar = c(5, 4, 4, 2) + 0.1)
      plot(seq(1:length(size)), MCR, xaxt = "n", col = "blue", cex = 1.2, ylim = c(0.05, (max(MCR)+max(se))),las=1,
          xlab = "number of genes", ylab = "Misclassification error", main = "")
      lines(seq(1:length(size)), MCR, col = "grey")
      
      segments(seq(1:length(size)), MCR - se, seq(1:length(size)), MCR + se, col='blue')
      axis(1, seq(1:length(size)), rev(size), las=2)
      # indicate chosen gene set size
      o <- rev(size) == length(x$selected.vars)
                    # alternatve:  # indicate minimum
                    #              o <- MCR == min(MCR)
      points(min(seq(1:length(size))[o]), MCR[o][1], pch = "x", cex = 2)
}
