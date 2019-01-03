#' NSAIDS in acute pain
#'
#' Meta-analysis on the effectiveness of topical non-steroidal
#' anti-inflammatory drugs (NSAIDS) in acute pain.
#'
#' Treatment success is defined as a reduction in pain of at least
#' 50\%.
#'
#' @docType data
#' 
#' @usage Moore1998
#'
#' @format A data frame with the following columns:
#' 
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study number \cr
#' \bold{\emph{succ.e}}\tab number of treatment successes in NSAIDS group \cr
#' \bold{\emph{nobs.e}}\tab number of patients in NSAIDS group \cr
#' \bold{\emph{succ.c}}\tab number of treatment successes in control group \cr
#' \bold{\emph{nobs.c}}\tab number of patients in control group \cr
#' }
#' 
#' @source
#'
#' Moore RA, Tramer MR, Carroll D, Wiffen PJ, McQuay HJ (1998):
#' Quantitive systematic review of topically applied non-steroidal
#' anti-inflammatory drugs.
#' \emph{British Medical Journal},
#' \bold{316}, 333--8
#'
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'               data = Moore1998, sm = "OR", method = "Inverse")
#'
#' print(limitmeta(m1), digits = 2)


"Moore1998"
