\name{pamClass}
\alias{pamClass}
\title{Classify using Prediction Analysis for MicroArrays}
\description{
  Classify using the Prediction Analysis for MicroArrays (PAM) algorithm as implemented 
  in the pamr package
}
\usage{
pamClass(object, groups, probe2gene = TRUE)
}
\arguments{
  \item{object}{object containing the expression measurements; currently the
    only method supported is one for ExpressionSet objects}
  \item{groups}{character string indicating the column containing the class membership}
  \item{probe2gene}{logical; if \code{TRUE} Affymetrix probeset IDs are translated
    into gene symbols; if \code{FALSE} no such translation is conducted}
}
\value{
  object of class \code{pamClass}
}
\references{
  Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and
     Gilbert Chu (1999). Diagnosis of multiple cancer types by shrunken
     centroids of gene expression.  PNAS 99: 6567-6572.  
     
     Available at \url{www.pnas.org}
  
  Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
     Microarrays, Chapman \& Hall/CRC, p. 221.
}
\author{Willem Talloen}
\seealso{\code{\link[pamr]{pamr.train}}}
% TODO add examples section
\keyword{models}
