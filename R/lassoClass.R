#' Classify using the Lasso
#' @param object object containing the expression measurements; currently the
#' only method supported is one for ExpressionSet objects
#' @param groups character string indicating the column containing the class membership
#' @return object of class \code{glmnet}
#' @references Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'  Microarrays, Chapman \& Hall/CRC, pp. 183, 205 and 212.
#' @author Willem Talloen
#' @seealso \code{\link[glmnet]{glmnet}}
#' @examples 
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#' 
#'   resultLasso <- lassoClass(object = ALL, groups = "BTtype")
#'   plot(resultLasso, label = TRUE,
#'     main = "Lasso coefficients in relation to degree of
#'   penalization.")
#'   topTable(resultLasso, n = 15)
#' }
#' @importFrom Biobase exprs pData featureData
#' @importFrom glmnet glmnet
#' @keywords models
#' @export
lassoClass <- function(object, groups){
  labels <- factor(pData(object)[,groups])
  object <- object[,!is.na(labels)]
  labels <- factor(pData(object)[,groups])
  
  fit <- glmnet(t(exprs(object)), labels, family="binomial", alpha = 1)
  fit$featureData <- pData(featureData(object))
  return(fit) 
}
