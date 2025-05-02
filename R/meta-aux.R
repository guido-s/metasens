## Auxiliary functions
##
## Package: metasens
## Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
## License: GPL (>= 2)
##

allNA <- function(x)
  all(is.na(x))

catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)

replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    return(replace)
  x
}

replaceNA <- function(x, replace = NA) {
  if (is.null(x))
    return(x)
  else
    x[is.na(x)] <- replace
  x
}

warnarg <- function(x, y, fun, cl, otherarg) {
  if (x %in% y)
    if (!missing(cl))
      warning("Argument '", x, "' has been removed from R function ", fun,
              ".\nThis argument can be used in R function ", cl, ".",
              call. = FALSE)
    else if (!missing(otherarg))
      warning("Argument '", x, "' has been replaced by argument '", otherarg,
              "' in R function ", fun, ".\nSee help page of R function ",
              fun, " for information on the use of the new argument.",
              call. = FALSE)
  ##
  invisible(NULL)
}

deprecated <- function(newvar, newmiss, args, old, warn = TRUE) {
  ##
  new <- deparse(substitute(newvar))
  ##
  if (length(args) == 0)
    return(newvar)
  ##
  if (is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (!is.na(charmatch(old, additional.arguments)))
    if (!newmiss) {
      if (warn)
        warning("Deprecated argument '", old, "' ignored as ",
                "'", new, "' is also provided.",
                call. = FALSE)
      return(newvar)
    }
    else {
      if (warn)
        warning("Use argument '", new, "' instead of '",
                old, "' (deprecated).",
                call. = FALSE)
      return(args[[charmatch(old, additional.arguments)]])
    }
  else
    return(newvar)
}

deprecated2 <- function(newvar, newmiss, oldvar, oldmiss, warn = TRUE,
                        oldtxt = NULL) {
  ##
  new <- deparse(substitute(newvar))
  if (is.null(oldtxt))
    oldtxt <- deparse(substitute(oldvar))
  ##
  if (newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & !oldmiss) {
    if (warn)
      warning("Deprecated argument '", oldtxt, "' ignored as ",
              "'", new, "' is also provided.",
              call. = FALSE)
    return(newvar)
  }
  else if (newmiss & !oldmiss) {
    if (warn)
      warning("Use argument '", new, "' instead of '",
              oldtxt, "' (deprecated).",
              call. = FALSE)
    return(oldvar)
  }
}

cond <- function(x, only.finite = TRUE, digits = 2, big.mark = "") {
  if (is.null(x))
    return(x)
  ##
  if (only.finite)
    x <- x[is.finite(x)]
  ##
  paste(formatN(unique(round(x, digits = digits)), digits = digits,
                big.mark = big.mark), collapse = ", ")
}

expandvar <- function(x, n, length = NULL) {
  res <- x
  if (!is.null(length))
    lenOK <- length(x) == length
  else
    lenOK <- TRUE
  ##
  if (lenOK & length(x) != n)
    res <- rep(x, rep_len(n, length(x)))
  ##
  res
}
