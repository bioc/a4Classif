\name{rfClass}
\alias{rfClass}
\alias{plot.rfClass}
\title{Classify using Random Forests}
\description{
  Classify using the Random Forest algorithm of Breiman (2001)
}
\usage{
rfClass(object, groups, probe2gene = TRUE)
}
\arguments{
  \item{object}{object containing the expression measurements; currently the
    only method supported is one for ExpressionSet objects}
  \item{groups}{character string indicating the column containing the class 
    membership}
  \item{probe2gene}{logical; if \code{TRUE} Affymetrix probeset IDs are translated
    into gene symbols in the output object; if \code{FALSE} no such translation 
    is conducted}
}
\value{
  Object of class 'rfClass'
}
% TODO add reference book Goehlmann and Talloen
\references{
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
     5-32.
}
\author{Tobias Verbeke and Willem Talloen}
\note{
  topTable and plot methods are available for 'rfClass' objects.
}
\seealso{\code{\link[randomForest]{randomForest}}}
% TODO add examples section
\keyword{models}
