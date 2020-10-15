#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname){
  options(error = NULL)
  packageStartupMessage(paste("\na4Classif version ", packageDescription("a4Classif")$Version, 
          "\n", sep = ""))
}

