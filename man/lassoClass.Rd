\name{lassoClass}
\alias{lassoClass}
\title{Classify using the Lasso}
\description{
  Classify using the Lasso algorithm as implemented in the glmnet package
}
\usage{
lassoClass(object, groups)
}
\arguments{
  \item{object}{object containing the expression measurements; currently the
    only method supported is one for ExpressionSet objects}
  \item{groups}{character string indicating the column containing the class membership}
}
\value{
  object of class \code{glmnet}
}
\references{Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
    Microarrays, Chapman \& Hall/CRC, pp. 183, 205 and 212.}
\author{Willem Talloen}
\seealso{\code{\link[glmnet]{glmnet}}}
\examples{

if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

  resultLasso <- lassoClass(object = ALL, groups = "BTtype")
  plot(resultLasso, label = TRUE,
    main = "Lasso coefficients in relation to degree of
  penalization.")
\dontrun{
  featResultLasso <- topTable(resultLasso, n = 15)
}
}
}
\keyword{models}
