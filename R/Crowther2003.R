#' Aspirin after Myocardial Infarction
#'
#' Meta-analysis on phenobarbital prior to preterm birth for
#' preventing neonatal periventricular haemorrhage
#'
#' @docType data
#' 
#' @usage Crowther2003
#'
#' @format A data frame with the following columns:
#' 
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{pvh.e}}\tab number of periventricular haemorrhages in experimental group \cr
#' \bold{\emph{n.e}}\tab number of observations in experimental group \cr
#' \bold{\emph{pvh.c}}\tab number of periventricular haemorrhages in control group \cr
#' \bold{\emph{n.c}}\tab number of observations in control group \cr
#' }
#'
#' @source
#'
#' Crowther CA, Henderson-Smart DJ (2003):
#' Phenobarbital prior to preterm birth for preventing neonatal
#' periventricular haemorrhage.
#' \emph{Cochrane Database of Systematic Reviews}, CD000164
#'
#' @examples
#' data(Crowther2003)
#' metabin(pvh.e, n.e, pvh.c, n.c, data = Crowther2003, studlab = study)


"Crowther2003"
