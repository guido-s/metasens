.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'metasens' package (version ",
           utils::packageDescription("metasens")$Version,
           ").",
           "\nReaders of 'Meta-Analysis with R (Use R!)' should install",
           "\nolder version of 'metasens' package: ",
           "https://tinyurl.com/jhv33bfm")           
  packageStartupMessage(msg)
}
