##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../metasens")


##
## (3) Build R package and PDF file with help pages
##
build("../metasens")
build_manual("../metasens")


##
## (4) Install R package
##
install("../metasens")


##
## (5) Check R package
##
check("../metasens")


##
## (6) Check examples
##
setwd("..")
library(numDeriv)
run_examples("metasens", run_dontrun = TRUE, run_donttest = TRUE)
warnings()
