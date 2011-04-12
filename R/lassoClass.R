lassoClass <- function(object, groups){
  labels <- factor(pData(object)[,groups])
  object <- object[,!is.na(labels)]
  labels <- factor(pData(object)[,groups])
  
  fit <- glmnet(t(exprs(object)), labels, family="binomial", alpha = 1)
  fit$featureData <- pData(featureData(object))
  return(fit) 
}
