.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'copas' package (version ",
               utils::packageDescription("copas")$Version,
               ").", sep="")
  packageStartupMessage(msg)
}
