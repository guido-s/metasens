summary.limitmeta <- function(object, ...){
  
  if (!inherits(object, "limitmeta"))
    stop("Argument 'object' must be an object of class \"limitmeta\"")
  
  res <- object
  class(res) <- c("summary.limitmeta", "limitmeta")
  
  res
}
