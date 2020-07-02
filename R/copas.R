#' Copas selection model analysis
#' 
#' Perform a Copas selection model analysis for selection bias in
#' meta-analysis.
#' 
#' The program takes an object of class \code{meta}, which is most
#' easily created by an analysis using one of the functions
#' \code{metabin}, \code{metacont} and \code{metagen} in the package
#' meta, performs a 'Copas selection model analysis' and presents a
#' graphical and tabular summary of the results. An object of class
#' \code{copas} is created and this can be used to recreate the
#' results table and graphs subsequently, without re-running the
#' analysis, using the \code{print}, \code{summary} and \code{plot}
#' function.
#' 
#' Conduct a Copas selection model analysis to investigate, and
#' attempt to correct for, selection / publication bias in a
#' meta-analysis.
#' 
#' The Copas selection model consists of two models, which are fitted
#' jointly.  The first is the usual random effects meta-analysis
#' model, and the second is a selection model, where study i is
#' selected for publication if Z>0, where
#' 
#' Z = gamma0 + gamma1 / (SE(i)) + delta(i)
#' 
#' The error delta(i) is correlated with the error in the random
#' effects meta-analysis, with correlation rho. If rho=0, the model
#' corresponds to the usual random effects meta-analysis. As rho moves
#' from 0 to 1, studies with larger treatment estimates are more
#' likely to be selected/published.
#' 
#' The software chooses a grid of gamma0 and gamma1 values,
#' corresponding to a range of selection / publication probabilities
#' for the study with the largest treatment effect standard error
#' (often the smallest study). For each value in this grid, the
#' treatment effect is estimated using the function \code{optim}. This
#' information is used to produce the contour plot (top right panel of
#' output from \code{plot.copas}).
#' 
#' Contours of constant treatment effect are usually locally
#' parallel. The software estimates the slope of these contours, and
#' combines this information with other parameter estimates from the
#' model to explore (i) how the treatment estimate, and its standard
#' error, change with increasing selection (bottom left panel,
#' \code{plot.copas}) and (ii) how much selection needs to be
#' accounted for before any remaining asymmetry in the funnel plot is
#' likely to have occurred by chance (bottom right panel,
#' \code{plot.copas}).
#' 
#' A table of results can be produced by the function
#' \code{summary.copas}. A more detail output is provided by the
#' function \code{print.copas}.
#' 
#' For a fuller description of the model, our implementation and
#' specifically our approach to estimating the locally parallel
#' contours, see Carpenter et al. (2009) and Schwarzer et al. (2010).
#' 
#' @param x An object of class \code{meta}, obtained from one of the
#'   functions \code{metabin}, \code{metacont} and \code{metagen} in
#'   the package meta.
#' @param gamma0.range (Advanced users only) A numerical vector of
#'   length two specifying the range of gamma0 values the program will
#'   explore.
#' 
#' The parameter gamma0 is the constant in the probit selection model
#' for study publication. Thus, the cumulative normal of gamma0 is
#' approximately the probability that a small study is published (in
#' non-technical terms gamma0 relates to the probability of publishing
#' a small study, although its values are not restricted to the range
#' [0,1]; larger values correspond to higher probabilities of
#' publishing a small study). Most users will not need to specify a
#' range for this parameter. When no argument is specified, the
#' program uses an algorithm to determine a suitable range. This is
#' based on the range of treatment effect standard errors in the
#' meta-analysis, and is described in more detail below.
#' @param gamma1.range (Advanced users only) A numerical vector of
#'   length two specifying the range of gamma1 values the program will
#'   explore.
#' 
#' The parameter gamma1 is the coefficient of study precision
#' (1/standard error) in the probit selection model for study
#' publication (in non-technical terms gamma1 relates to the rate at
#' which the probability of publishing a study increases as the
#' standard error of the treatment effect it reports decreases; larger
#' values correspond to higher probabilities of publishing a small
#' study). Most users will not need to specify a range for this
#' parameter. When no argument is specified, the program uses an
#' algorithm to determine a suitable range.  This is based on the
#' range of treatment effect standard errors in the meta-analysis, and
#' is described in more detail below.
#' @param ngrid The program fits the Copas selection model over a grid
#'   defined by the range of values of gamma0 and gamma1 specified in
#'   the previous two arguments. This parameter fixes the square-root
#'   of the number of points in the grid.
#' @param nlevels (Advanced users only). Fitting the Copas model over
#'   the grid specified by the previous three arguments results in a
#'   treatment estimate at every point in the grid. These can then be
#'   displayed on a contour plot where contours of treatment effect
#'   (z-axis) are shown by gamma0 (x-axis) and gamma1 (y-axis). This
#'   argument specifies the number of contour lines that will be
#'   drawn.
#' 
#' \bold{Note}
#' 
#' (i) Calculations for the contour plot are performed by the function
#' \code{copas}, so this argument has no effect in the \code{plot}
#' function.
#' 
#' (ii) If a large number of contour lines are desired, then you may
#' wish to consider increasing the grid size (argument \code{ngrid}
#' above).
#' 
#' Leave this option unspecified if you are using the option
#' \code{levels} below.
#' @param levels A numerical vector of treatment values for which
#'   contour lines will be drawn. In more detail, fitting the Copas
#'   model over the grid specified by the arguments
#'   \code{gamma0.range}, \code{gamma1.range} and \code{ngrid} results
#'   in a treatment estimate at every point in the grid.  These are
#'   then displayed on a contour plot where contours of treatment
#'   effect (z-axis) are shown by gamma0 (x-axis) and gamma1
#'   (y-axis). This argument is a numerical vector which specifies the
#'   treatment effects for which contour lines will be drawn.
#' 
#' It is usually not a good idea to set this argument for initial
#' runs, as one does not know the range of treatment values that the
#' contour plot will cover, and treatment values which do not
#' correspond to values in the contour plot (defined by the range of
#' gamma0 and gamma1) will not be plotted.
#' 
#' \bold{Note}
#' 
#' (i) Calculations for the contour plot are performed by the function
#' \code{copas}, so this argument has no effect in the \code{plot}
#' function.
#' 
#' (ii) Contours will not be drawn if a large number of contour lines
#' are desired, then you may wish to consider increasing the grid size
#' (argument \code{ngrid} above).
#' 
#' Leave this option unspecified if you are using the option
#' \code{nlevels} above.
#' @param slope A numeric providing the slope of the line
#'   approximately orthogonal to contours in the contour plot. If the
#'   argument \code{slope} is \code{NULL} (default) the program seeks
#'   to estimate the slope of the contours in the region of the
#'   maximum, which are usually approximately parallel. Most users
#'   will leave the argument \code{slope} unspecified, at least for
#'   the first analysis of a data set, but in certain cases setting it
#'   manually can improve the results.
#' @param left A logical indicating whether the cause of any selection
#'   bias is due to missing studies on the left or right of the funnel
#'   plot: left hand side if \code{left=TRUE}, right hand side if
#'   \code{left=FALSE}. This information is needed in order to be sure
#'   the test for presence of residual selection bias is calculated
#'   correctly. If not set, the linear regression test for funnel plot
#'   asymmetry (i.e., function \code{metabias(..., meth="linreg")}) is
#'   used to determine whether studies are missing on the left or
#'   right hand side. In the majority of cases this will work
#'   correctly.
#' @param rho.bound (Advanced users only) A number giving the upper
#'   bound for the correlation parameter \code{rho} (see details
#'   below). This must be < 1, and usually > 0.95. The lower bound is
#'   calculated as -(the upper bound).
#' @param sign.rsb The significance level for the test of residual
#'   selection bias (between 0 and 1).
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param silent A logical indicating whether information on progress
#'   in fitting the Copas selection model should be printed:
#'   \code{silent=TRUE}, do not print information (the default);
#'   \code{silent=FALSE}, print information.
#' @param warn A number setting the handling of warning messages. It
#'   is not uncommon for numerical problems to be encountered during
#'   estimation over the grid of (gamma0, gamma1) values. Usually this
#'   does not indicate a serious problem. This option specifies what
#'   to do with warning messages.  \code{warn=-1}: ignore all
#'   warnings; \code{warn=0} (the default): store warnings till
#'   function finishes; if there are less than 10, print them,
#'   otherwise print a message saying warning messages were generated;
#'   \code{warn=1}: print warnings as they occur; \code{warn=2}: stop
#'   the function when the first warning is generated. For further
#'   details see \code{help(options)}.
#' 
#' @return An object of class \code{copas} with corresponding
#'   \code{print}, \code{summary}, and \code{plot} function. The
#'   object is a list containing the following components:
#' \item{TE}{Vector of treatment effects plotted in treatment effect
#'   plot}
#' \item{seTE}{Vector of standard error of \code{TE}}
#' \item{TE.random}{Usual random effects estimate of treatment effect}
#' \item{seTE.random}{Usual standard error of \code{TE.random}}
#' \item{left}{Whether selection bias expected on left or right}
#' \item{rho.bound}{Bound on \code{rho}}
#' \item{gamma0.range}{Range of gamma0 (see help on \code{copas}
#'   arguments above)}
#' \item{gamma1.range}{Range of gamma1 (see help on \code{copas}
#'   arguments above)}
#' \item{slope}{Slope of line approximately orthogonal to contours in
#'   contour plot}
#' \item{regr}{A list containing information on regression lines
#'   fitted to contours in contour plot}
#' \item{ngrid}{Square root of grid size}
#' \item{nlevels}{Number of contour lines}
#' \item{gamma0}{Vector of gamma0 values at which model fitted
#'   (determined by gamma0.range and grid). x-axis values for contour
#'   plot}
#' \item{gamma1}{vector of gamma1 values at which model fitted
#'   (determined by gamma1.range and grid). y-axis values for contour
#'   plot}
#' \item{TE.contour}{Treatment values (ie z-axis values) used to draw
#'   contour plot.}
#' \item{x.slope}{x coordinates for 'orthogonal line' in contour plot}
#' \item{y.slope}{y coordinates for 'orthogonal line' in contour plot}
#' \item{TE.slope}{Vector of treatment values plotted in treatment
#'   effect plot}
#' \item{seTE.slope}{Standard error of \code{TE.slope}} 
#' \item{rho.slope}{Vector of estimated rho values corresponding to
#'   treatment estimates in \code{TE.slope}}
#' \item{tau.slope}{Vector of estimated heterogeneity values
#'   corresponding to treatment estimates in \code{TE.slope}}
#' \item{loglik1}{Vector of log-likelihood values corresponding to
#'   treatment estimates in \code{TE.slope}}
#' \item{conv1}{Numerical vector indicating convergence status for
#'   each treatment estimate in \code{TE.slope} - see parameter
#'   \code{convergence} in function \code{optim}}
#' \item{message1}{Character vector - translation of \code{conv1}}
#' \item{loglik2}{Vector of log-likelihoods from fitting model to
#'   evaluate presence of residual selection bias}
#' \item{conv2}{Numerical vector indicating convergence status for
#'   models to evaluate presence of residual selection bias - see
#'   parameter \code{convergence} in function \code{optim}}
#' \item{message2}{Character vector - translation of \code{conv2}}
#' \item{publprob}{Vector of probabilities of publishing the smallest
#'   study, used in x-axis of bottom two panels in function
#'   \code{plot.copas}}
#' \item{pval.rsb}{P-values for tests on presence of residual
#'   selection bias, plotted in bottom right panel in
#'   \code{plot.copas}}
#' \item{sign.rsb}{The significance level for the test of residual
#'   selection bias}
#' \item{N.unpubl}{Approximate number of studies the model suggests
#'   remain unpublished}
#' \item{sm}{Effect measure (e.g., for binary data, OR - odds ratio,
#'   RR - risk ratio, RD - risk difference, AS - arcsin difference)}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.} \item{call}{Call to \code{copas}
#'   function}
#' \item{version}{Version of R package metasens used to create
#'   object.}
#' \item{x}{Details of meta-analysis object used as input into
#'   \code{copas} function}
#' 
#' @author James Carpenter \email{James.Carpenter@@lshtm.ac.uk}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{plot.copas}}, \code{\link{summary.copas}},
#' \code{\link[meta]{metabias}}, \code{\link[meta]{metagen}},
#' \code{\link[meta]{funnel}}
#' 
#' @references
#'
#' Carpenter JR, Schwarzer G, Rücker G, Künstler R (2009):
#' Empirical evaluation showed that the Copas selection model provided
#' a useful summary in 80\% of meta-analyses.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{62}, 624--31
#' 
#' Copas J (1999):
#' What works?: Selectivity models and meta-analysis.
#' \emph{Journal of the Royal Statistical Society, Series A},
#' \bold{162}, 95--109
#' 
#' Copas J, Shi JQ (2000):
#' Meta-analysis, funnel plots and sensitivity analysis.
#' \emph{Biostatistics}, \bold{1}, 247--62
#' 
#' Copas JB, Shi JQ (2001):
#' A sensitivity analysis for publication bias in systematic reviews.
#' \emph{Statistical Methods in Medical Research},
#' \bold{10}, 251--65
#' 
#' Schwarzer G, Carpenter J, Rücker G (2010):
#' Empirical evaluation suggests Copas selection model preferable to
#' trim-and-fill method for selection bias in meta-analysis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{63}, 282--8
#' 
#' @examples
#' data(Fleiss93)
#' 
#' # Perform meta-analysis
#' #  (Note event.e indicates events, n.e total in exposed arm;
#' #        event.c indicates events, n.c total in control arm)
#' #
#' m1 <- metabin(event.e, n.e, event.c, n.c, data = Fleiss93, sm = "OR")
#' summary(m1)
#' 
#' # Perform a basic Copas selection model analysis
#' #
#' cop1 <- copas(m1)
#' plot(cop1)
#' summary(cop1)
#' #
#' # Interpretation: 
#' #
#' # a. The initial meta-analysis shows the fixed and random effects
#' #    pooled ORs differ; consistent with asymmetry in the funnel
#' #    plot and possible selection bias. Both fixed effect and random
#' #    effects model show a significant treatment effect in this
#' #    dataset.
#' #
#' # b. Plotting the copas analysis shows
#' #
#' # (i) funnel plot: asymmetry indicates possible selection bias.
#' #
#' # (ii) contour plot treatment effect declines steadily as selection
#' #      increases (no selection, top right, log OR < -0.12;
#' #      increasing selection as move to left of plot, log OR rises
#' #      to -0.03.
#' #
#' 
#' # (iii) Treatment effect plot suggests that even with no selection,
#' #       p-value for treatment effect is larger than 0.05 which is
#' #       different from the result of the usual random effects model
#' #       (see output of summary(cop1). This difference is due to the
#' #       use of different methods to estimate the between-study
#' #       variance: maximum-likelihood in Copas analysis compared to
#' #       method-of-moments in usual random effects model.  The
#' #       p-value for treatment effect is increasing with increasing
#' #       selection.
#' #
#' # (iv) P-value for residual selection bias plot: this shows that
#' #      even with no selection bias, the p-value for residual
#' #      selection bias is non-significant at the 10% level. As
#' #      expected, as selection increases the p-value for residual
#' #      selection bias increases too.
#' 
#' 
#' # Repeat the same example, setting several arguments of the copas
#' # function:
#' #
#' cop2 <- copas(m1,
#'               gamma0.range = c(-0.5, 2.1), # range of gamma0 parameter
#'               gamma1.range = c(0, 0.08),   # range of gamma1 parameter
#'               ngrid = 20,                  # specify a 20x20 grid (finer than default)
#'               levels = c(-0.13, -0.12, -0.1, -0.09,
#'                          -0.07, -0.05, -0.03), # specify contour lines
#'               slope = 0.2,       # specify slope of 'orthogonal' line in contour plot
#'               left = FALSE,      # as any selection bias due to missing studies on right
#'               rho.bound = 0.998, # constrain rho between [-0.998, 0.998]
#'               silent = FALSE,    # update user on progress
#'               warn = -1          # suppress warning messages
#'              )
#' plot(cop2)
#' #
#' # Print table of results used to draw treatment effect plot:
#' #
#' summary(cop2)
#' @export copas
#'
#' @importFrom meta metabias metagen
#' @importFrom grDevices contourLines
#' @importFrom stats lm optim pchisq pnorm vcov


copas <- function(x,
                  gamma0.range = NULL,
                  gamma1.range = NULL,
                  ngrid = 20,
                  nlevels = 10,
                  levels = NULL,
                  slope = NULL,
                  left = NULL,
                  rho.bound = 0.9999,
                  sign.rsb = 0.1,
                  backtransf = x$backtransf,
                  silent = TRUE,
                  warn = options()$warn) {
  
  meta:::chkclass(x, "meta")
  
  if (!is.numeric(rho.bound) && (rho.bound <=0 | rho.bound > 1))
    stop("no valid value for 'rho.bound'")
  
  if (!is.null(slope)) {
    if (!is.numeric(slope))
      stop("Argument 'slope' must be numeric")
    if (length(slope) > 1) {
      warning(paste("Argument 'slope' must be of length 1;",
                    "first element of vector is used"))
      slope <- slope[1]
    }
  }
  
  
  ## Check significance level for test of residual selection bias
  ##
  meta:::chklevel(sign.rsb)
  
  
  oldopt <- options(warn = warn)
  on.exit(options(oldopt))
  
  
  TE <- x$TE
  seTE <- x$seTE
  sel <- !is.na(TE) & !is.na(seTE)
  if (length(TE) != sum(sel))
    warning(paste(length(TE) - sum(sel),
                  "observation(s) dropped due to missing values"))
  TE <- TE[sel]
  seTE <- seTE[sel]
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  ##
  tau <- x$tau
  ##
  seTE.min <- min(seTE)
  seTE.max <- max(seTE)
  
  
  if (!silent) {
    cat("\n\n")
    cat("====================================\n")
    cat("========== COPAS ANALYSIS ==========\n")
    cat("====================================\n\n")
    ##
    ## print:
    ## (a) fixed effect analysis
    ## (b) random effect analysis
    ## (c) test for heterogeneity, using appropriate function
    ##  
    cat("\n")
    cat("1) Summary statistics and test for heterogeneity\n")
    cat("================================================\n")
    print(summary(x))
  }
  
  
  ## calculate and display the range of gamma0,
  ## gamma1 used subsequently,
  ## unless options given by user override these
  ##
  if (is.null(gamma1.range)) {
    gamma1.range <- c(0, 1.29 / (1 / seTE.min - 1 / seTE.max))
  }
  if (is.null(gamma0.range))
    gamma0.range <- c(-0.25 - gamma1.range[2] / seTE.max, 2)
  ##
  ##
  gamma0 <- seq(gamma0.range[1], gamma0.range[2], length = ngrid)
  gamma1 <- seq(gamma1.range[1], gamma1.range[2], length = ngrid)
  ##
  gamma0.min <- min(gamma0)
  gamma0.max <- max(gamma0)
  gamma1.min <- min(gamma1)
  gamma1.max <- max(gamma1)
  ##
  ##
  if (!silent) {
    ##
    cat("\n")
    cat("2) Ranges of gamma0, gamma1 used for calculating contour plot\n")
    cat("=============================================================\n")
    cat("gamma0 ranges from ",
        round(gamma0.range[1], 2),
        " to ",
        round(gamma0.range[2],2),
        "\n")
    cat("gamma1 ranges from ",
        round(gamma1.range[1], 2),
        " to ",
        round(gamma1.range[2],2),
        "\n")
    ##
    ##
    publprob.seTE.max <- range(pnorm(gamma0 + gamma1 / seTE.max))
    publprob.seTE.min <- range(pnorm(gamma0 + gamma1 / seTE.min))
    ##
    cat(paste("Range of probability publishing trial with largest SE:  (",
              round(publprob.seTE.max[1], 3),", ",
              round(publprob.seTE.max[2], 3),")",sep = "") ,"\n")
    cat(paste("Range of probability publishing trial with smallest SE: (",
              round(publprob.seTE.min[1], 3),", ",
              round(publprob.seTE.min[2], 3),")",sep = "") ,"\n\n")
  }
  
  
  ## calculate the contour plot
  ##
  if (!silent) {
    cat("\n")
    cat("3) Starting calculations for contour plot\n")
    cat("=========================================\n")
  }
  ##
  if (is.null(left))
    left <- as.logical(sign(metabias(x, meth = "linreg", k.min = 3)$estimate[1]) == 1)
  ##
  if (left)
    rho0 <-  rho.bound / 2
  else
    rho0 <- -rho.bound / 2
  
  
  ##
  TE.contour <- matrix(NA, nrow = ngrid, ncol = ngrid)
  ##
  for (i in seq(along = gamma0)) {
    for (j in seq(along = gamma1)) {
      ##
      try(junk0 <- optim(c(TE.random, rho0, tau),
                         fn = copas.loglik.without.beta,
                         gr = copas.gradient.without.beta,
                         lower = c(-Inf, -rho.bound,   0),
                         upper = c( Inf,  rho.bound, Inf),
                         gamma = c(gamma0[i], gamma1[j]),
                         TE = TE, seTE = seTE,
                         method = "L-BFGS-B"),
          silent = TRUE)
      ##
      TE.contour[i, j] <- junk0$par[1]
    }
    if (!silent) {
      cat(paste(round(100 * i * j / (ngrid * ngrid), 0), "%, ",sep = "" ))
    }
  }
  ##
  if (!silent) {
    cat("Done\n\n")
  }
  
  
  ## calculations to approximate route of orthogonal line
  ##
  if (!silent) {
    cat("\n")
    cat("4) Calculating approximate orthogonal line\n")
    cat("==========================================\n")
  }
  ##
  ## the approach of lowess smoothing the smallest distances
  ## gives very bendy curves: try instead to calculate gradient of
  ## each contour, then average (-1 / .) and then draw straight line through
  ## top right value.
  ##
  gamma0.rescale <- (gamma0 - gamma0.min) / (gamma0.max - gamma0.min)
  gamma1.rescale <- (gamma1 - gamma1.min) / (gamma1.max - gamma1.min)
  ##
  if (is.null(levels))
    junk <- contourLines(x = gamma0.rescale,
                         y = gamma1.rescale,
                         z = TE.contour,
                         nlevels = nlevels)
  else
    junk <- contourLines(x = gamma0.rescale,
                         y = gamma1.rescale,
                         z = TE.contour,
                         levels = levels)
  ##
  levels <- as.numeric(unlist(junk)[names(unlist(junk)) == "level"])
  nlevels <- length(levels)
  ##
  nobs <- rep(NA, nlevels)
  adj.r.squareds <- rep(NA, nlevels)
  slopes <- rep(NA, nlevels)
  se.slopes <- rep(NA, nlevels)
  intercepts <- rep(NA, nlevels)
  ##
  for(l in 1:(nlevels)) {
    lm.op <- lm(junk[[l]]$y ~ junk[[l]]$x)
    nobs[l] <- length(junk[[l]]$x)
    adj.r.squareds[l] <- summary(lm.op)$adj.r.squared
    slopes[l] <- lm.op$coef[2]
    se.slopes[l] <- sqrt(diag(vcov(lm.op))[2])
    intercepts[l] <- lm.op$coef[1]
  }

  ## calculate crossings of contours with approximate orthogonal line
  ## this provides the points to estimate the effect
  ##
  adj.r.squareds[is.nan(adj.r.squareds)] <- -100
  ##
  sel <- adj.r.squareds > 0
  ##
  if (is.null(slope)) {
    ##
    if (all(!sel)) {
      warning("No contour line with corresponding adjusted r.square ",
              "larger than zero")
      slope <- NA
    }
    else
      slope <- -1 / metagen(slopes[sel], 1 / sqrt(nobs)[sel],
                            method.tau.ci = "")$TE.fixed
  }
  ##
  ##x.slope <- ((1 - slope - intercepts ) / (slopes - slope))[sel]
  ##
  x.slope <- ((1 - slope - intercepts ) / (slopes - slope))
  y.slope <- 1 + slope * (x.slope - 1)
  
  
  if (!silent) {
    cat("Done\n\n")
  }
  
  
  ## calculations for plot how mean / se changes as
  ## come down from the maximum
  ##
  if (!silent) {
    cat("\n")
    cat("5) Calculating TE and seTE, as come down slope\n")
    cat("==============================================\n")
  }
  ##
  gamma0.slope <- x.slope * (gamma0.max - gamma0.min) + gamma0.min
  gamma1.slope <- y.slope * (gamma1.max - gamma1.min) + gamma1.min
  ##
  ## Select only those values within the range of gamma0 and gamma1
  ##
  sel0 <- (gamma0.slope >= min(gamma0.range) &
           gamma0.slope <= max(gamma0.range))
  sel1 <- (gamma1.slope >= min(gamma1.range) &
           gamma1.slope <= max(gamma1.range))
  ##
  x.slope <- x.slope[sel0&sel1]
  y.slope <- y.slope[sel0&sel1]
  ##
  gamma0.slope <- gamma0.slope[sel0&sel1]
  gamma1.slope <- gamma1.slope[sel0&sel1]
  ##
  ## Reorder x.slope, y.slope, gamma0.slope, gamma1.slope
  ## according to publprob.seTE.max
  ##
  ord <- rev(order(pnorm(gamma0.slope + gamma1.slope / seTE.max)))
  ##
  x.slope      <- x.slope[ord]
  y.slope      <- y.slope[ord]
  gamma0.slope <- gamma0.slope[ord]
  gamma1.slope <- gamma1.slope[ord]
  ##
  ## Add the point (gamma0 = 10, gamma1 = 0) representing the usual random
  ## effects model with no allowance for selection
  ## (Copas, Shi, 2001, p.256)
  ##
  gamma0.slope <- c(10, gamma0.slope)
  gamma1.slope <- c( 0, gamma1.slope)
  ##
  ##
  sel2 <- !is.na(x.slope) & !is.na(y.slope)
  sel3 <- !is.na(gamma0.slope) & !is.na(gamma1.slope)
  ##
  x.slope <- x.slope[sel2]
  y.slope <- y.slope[sel2]
  ##
  gamma0.slope <- gamma0.slope[sel3]
  gamma1.slope <- gamma1.slope[sel3]
  ##
  n.gamma0.slope <- length(gamma0.slope)
  ##
  publprob <- pnorm(gamma0.slope + gamma1.slope / seTE.max)
  ##
  ##
  ##
  TE.slope   <- rep(NA, n.gamma0.slope)
  seTE.slope <- rep(NA, n.gamma0.slope)
  rho.slope  <- rep(NA, n.gamma0.slope)
  tau.slope  <- rep(NA, n.gamma0.slope)
  beta.slope <- rep(NA, n.gamma0.slope)
  ##
  loglik1  <- rep(NA, n.gamma0.slope)
  conv1    <- rep(NA, n.gamma0.slope)
  message1 <- rep("", n.gamma0.slope)
  ##
  for (i in seq(along = gamma1.slope)) {
    try(junk1 <- optim(c(TE.random, rho0, tau),
                       fn = copas.loglik.without.beta,
                       gr = copas.gradient.without.beta,
                       lower = c(-Inf, -rho.bound,   0),
                       upper = c( Inf,  rho.bound, Inf),
                       gamma = c(gamma0.slope[i], gamma1.slope[i]),
                       TE = TE, seTE = seTE,
                       method = "L-BFGS-B"),
        silent = TRUE)
    ##
    TE.slope[i]  <- junk1$par[1]
    rho.slope[i] <- junk1$par[2]
    tau.slope[i] <- junk1$par[3]
    ##
    loglik1[i]  <- junk1$value
    conv1[i]    <- junk1$convergence
    message1[i] <- junk1$message
    ##
    ## in case of singular hessian, do the best we can:
    ##
    try(junk2 <- optim(c(TE.random, rho0, tau),
                       fn = copas.loglik.without.beta,
                       gr = copas.gradient.without.beta,
                       lower = c(-Inf, -rho.bound,   0),
                       upper = c( Inf,  rho.bound, Inf),
                       gamma = c(gamma0.slope[i], gamma1.slope[i]),
                       TE = TE, seTE = seTE,
                       method = "L-BFGS-B", hessian = TRUE),
        silent = TRUE)
    ##
    ##print(junk2$hessian)
    ##
    try(seTE.slope[i] <-
        sqrt(solve(junk2$hessian + 0.00000001)[1, 1]),
        silent = TRUE)
    ##
    ## if this fails, take previous sd
    ##
    if ((i > 1 & is.na(seTE.slope[i])) ||
        (i > 1 & seTE.slope[i] == 0))
      seTE.slope[i] <- seTE.slope[i - 1]
    ##
    ## if that fails, try this!
    ##
    if (is.na(seTE.slope[i]) || seTE.slope[i] == 0)
      try(seTE.slope[i] <- sqrt(1 / junk2$hessian[1, 1]),
          silent = TRUE)
  }
  
  
  if (!silent) {
    cat("Done\n\n")
  }
  
  
  ## calculations for goodness of fit plot along orthogonal line
  ## (above), calculate log-likelihood in model containing sd
  ##
  if (!silent) {  
    cat("\n")
    cat("6) Calculating goodness of fit, as come down orthogonal line\n")
    cat("============================================================\n")
  }
  ##
  if (left)
    rho.lim <- c(0, rho.bound)
  else
    rho.lim <- c(-rho.bound, 0)
  ##
  ## beta constrained:
  ##
  TE.slope.bc   <- rep(NA, n.gamma0.slope)
  rho.slope.bc  <- rep(NA, n.gamma0.slope)
  tau.slope.bc  <- rep(NA, n.gamma0.slope)
  beta.slope.bc <- rep(NA, n.gamma0.slope)
  ##
  loglik2  <- rep(NA, n.gamma0.slope)
  conv2    <- rep(NA, n.gamma0.slope)
  message2 <- rep("", n.gamma0.slope)
  ##
  for (i in seq(along = gamma1.slope)) {
    try(junk3 <- optim(c(TE.random, rho0, tau, 0),
                       fn = copas.loglik.with.beta,
                       lower = c(-Inf, rho.lim[1], 0  , -Inf),
                       upper = c( Inf, rho.lim[2], Inf,  Inf),
                       gamma = c(gamma0.slope[i], gamma1.slope[i]),
                       TE = TE, seTE = seTE,
                       method = "L-BFGS-B"),
        silent = TRUE)
    ##
    TE.slope.bc[i]   <- try(junk3$par[1])
    rho.slope.bc[i]  <- try(junk3$par[2])
    tau.slope.bc[i]  <- try(junk3$par[3])
    beta.slope.bc[i] <- try(junk3$par[4])
    ##
    loglik2[i]  <- try(junk3$value)
    conv2[i]    <- try(junk3$convergence)
    message2[i] <- try(junk3$message)
  }
  ##
  ## P-value for test of residual selection bias:
  ##
  pval.rsb <- 1 - pchisq(2 * (loglik1 - loglik2), df = 1)

  
  if (!silent) {
    cat("Done\n\n")
  }
  
  
  ##
  ## compute no of unpublished studies
  ##
  N.unpubl <- rep(NA, length(publprob))
  ##
  for (i in seq(along = N.unpubl)) {
    p.si <- pnorm(gamma0.slope[i] + gamma1.slope[i] / seTE)
    N.unpubl[i] <- sum((1 - p.si) / p.si)
  }
  
  
  res <- list(TE = TE,
              seTE = seTE,
              TE.random = TE.random,
              seTE.random = seTE.random,
              left = left,
              rho.bound = rho.bound,
              gamma0.range = gamma0.range,
              gamma1.range = gamma1.range,
              slope = slope,
              regr = list(
                levels = levels,
                nobs = nobs,
                adj.r.squareds = adj.r.squareds,
                slopes = slopes,
                se.slopes = se.slopes,
                intercepts = intercepts),
              ngrid = ngrid,
              nlevels = nlevels,
              gamma0 = gamma0,
              gamma1 = gamma1,
              TE.contour = TE.contour,
              ##
              x.slope = x.slope,
              y.slope = y.slope,
              ##
              TE.slope = TE.slope,
              seTE.slope = seTE.slope,
              rho.slope = rho.slope,
              tau.slope = tau.slope,
              ##
              loglik1 = loglik1,
              conv1 = conv1,
              message1 = message1,
              ##
              ##TE.slope.bc = TE.slope.bc,
              ##rho.slope.bc = rho.slope.bc,
              ##tau.slope.bc = tau.slope.bc,
              ##beta.slope.bc = beta.slope.bc,
              ##
              loglik2 = loglik2,
              conv2 = conv2,
              message2 = message2,
              ##
              publprob = publprob,
              pval.rsb = pval.rsb,
              sign.rsb = sign.rsb,
              N.unpubl = N.unpubl,
              sm = x$sm, call = match.call(), x = x)
  
  res$version <- utils::packageDescription("metasens")$Version
  
  if (!is.null(x$title))
    res$title <- x$title
  if (!is.null(x$complab))
    res$complab <- x$complab
  if (!is.null(x$outclab))
    res$outclab <- x$outclab
  
  res$backtransf <- backtransf
  
  class(res) <- c("copas")
  
  res
}
