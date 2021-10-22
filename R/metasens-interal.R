.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'metasens' package (version ",
           utils::packageDescription("metasens")$Version,
           ").",
           "\nSupporting 'Meta-Analysis with R (Use R!)', ",
           "first edition.")
  packageStartupMessage(msg)
}
