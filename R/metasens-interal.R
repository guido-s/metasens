.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'metasens' package (version ",
           utils::packageDescription("metasens")$Version,
           ").",
           "\nType 'help(metasens)' for a brief overview.")
  packageStartupMessage(msg)
}
