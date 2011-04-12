\name{ROCcurve}
\alias{ROCcurve}
\title{Receiver operating curve}
\description{
A ROC curve plots the fraction of true positives (TPR = true positive rate)
 versus the fraction of false positives (FPR = false positive rate) for a binary classifier
 when the discrimination threshold is varied. Equivalently, one can also plot
 sensitivity versus (1 - specificity).
}
\usage{
ROCcurve(object, groups, probesetId = NULL, geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...)
}
\arguments{
  \item{object}{ExpressionSet object for the experiment}
  \item{groups}{String containing the name of the grouping variable. This should be a 
  the name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{probesetId}{The probeset ID. These should be stored in the \code{featureNames}
   of the \code{expressionSet} object.}
  \item{geneSymbol}{The gene symbol. These should be stored in the column \code{`Gene Symbol`}
   in the \code{featureData} of the \code{expressionSet} object.}
  \item{main}{Main title on top of the graph}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}
\references{Some explanation about ROC can be found on \url{http://en.wikipedia.org/wiki/ROC_curve}
 and \url{http://www.anaesthetist.com/mnm/stats/roc/Findex.htm}. The latter has
 at the bottom a nice interactive tool to scroll the cut-off and to see how it affects the FP/TP table and the ROC curve. 
}
\author{Willem Talloen}
\examples{
# simulated data set
esSim <- simulateData()
ROCcurve(probesetId = 'Gene.1', object = esSim, groups = 'type', addLegend = FALSE)

# ALL data set
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
  ROCres <- ROCcurve(gene = "ABL1", object = ALL, groups = "BTtype")
}
}
