updateversion <- function(x) {
  
  if (inherits(x, "meta")) {
    ##
    ## Update older meta objects
    ##
    if (update_needed(x$version, pkg = "meta"))
      x <- update(x, warn = FALSE, warn.deprecated = FALSE)
    ##
    return(x)
  }
  ##
  if (inherits(x, "orbbound")) {
    if (update_needed(x$version, 1, 5, pkg = "metasens")) {
      x$common <- x$fixed
      x$x <- updateversion(x$x)
    }
  }
  ##
  x
}


update_needed <- function(version,  major = 1, minor = 5,
                          verbose = FALSE,
                          pkg = NULL) {
  if (is.null(version)) {
    version <- 0.1
    major.cur <- 0
    minor.cur <- 1
  }
  else {
    version <- unlist(strsplit(version, "-")[1])
    major.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[1])
    minor.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[2])
  }
  ##
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < minor))
  if (res & verbose) {
    if (!is.null(pkg))
      addtext <- paste0(pkg, " ")
    else
      addtext <- ""
    ##
    message(paste0("Update to ", addtext, "version ", major, ".", minor))
  }
  ##
  res
}
