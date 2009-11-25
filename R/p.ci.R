p.ci <- function(lower, upper){
  lower <- rmSpace(lower)
  upper <- rmSpace(upper)
  ##
  ifelse (lower=="NA" & upper=="NA",
          "",
          paste(" [", format(lower, justify="right"),
                "; ", format(upper, justify="right"), "]", sep=""))
}
