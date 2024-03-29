## metasens, version 1.5-2 (2023-02-28)

### Major changes

* For Copas selection model, show value of between-study variance and
  standard deviation in printouts

### Bug fixes

* Arguments 'lfkindex' and 'xlim' were ignored in doiplot()

### User-visible changes

* doiplot():
  - print more informative label on horizontal axis for meta-analysis
    objects

### Internal changes

* lfkindex():
  - new list element 'x' with meta-analysis object used as input

* New branch 'release' on GitHub starting with **metasens**, version
  1.5-2


## metasens, version 1.5-1 (2022-12-21)

### User-visible changes

* Change maintainer's email address

### Internal changes

* Remove 'meta:::' call


## metasens, version 1.5-0 (2022-07-11)

### Major changes

* Use term 'common effect model' instead of 'fixed effect model' in
  the documentation and argument 'common' instead of 'fixed' to (not)
  show results for common effect model


## metasens, version 1.0-1 (2022-01-28)

### Bug fixes

* Calculation of adjusted standard error in limitmeta() for adjustment
  method beta0 (default):
  - use random effects instead of fixed effect weights
  - consider covariance between alpha and beta

### User-visible changes

* plot.copas():
  - new argument 'main' replacing deprecated argument 'caption'


## metasens, version 1.0-0 (2021-10-22)

### Major changes

* Behaviour of print and print.summary functions switched (to be in
  line with other print and print.summary functions in R)

* Renamed arguments:
  - 'fixed' (instead of 'comb.fixed')
  - 'random' (instead of 'comb.random')
  - 'level.ma' (instead of 'level.comb')

* In metamiss(), an infinite IMOR value is set to 1e10 instead of 9999

### User-visible changes

* metamiss():
  - arguments 'IMOR.e' and 'IMOR.c' can be specified for each
    individual study


## metasens, version 0.6-0 (2021-01-15)

### Major changes

* New functions doiplot() and lfkindex() implementing the Doi plot for
  asymmetry and the LFK index to test for asymmetry ([Furuya-Kanamori
  et al., 2018](https://doi.org/10.1097/XEB.0000000000000141))

* Confidence interval of limit meta-analysis estimate can be shown in
  funnel plot

### User-visible changes

* funnel.limitmeta():
  - new argument 'show.ci.adjust' to show confidence interval of
    adjusted estimate as a diamond


## metasens, version 0.5-0 (2020-09-29)

### Major changes

* Calculation of adjusted treatment estimate for Copas model in R
  function copas() instead of summary.copas()

### User-visible changes

* copas():
  - new argument 'level.comb' to calculate confidence interval for
    pooled estimates
  - new arguments 'title', 'complab' and 'outclab' to print
    information on systematic review / meta-analysis

* summary.copas() and plot.copas():
  - arguments 'sign.rsb' and 'level' removed

* print.copas():
  - argument 'sign.rsb' removed

* print.limitmeta(), print.orbbound(), print.summary.limitmeta():
  - argument 'digits.zval' renamed to 'digits.stat'

### Internal changes

* copas():
  - new list elements (with information for adjusted overall effect):
    'TE.adjust', 'seTE.adjust', 'lower.adjust', 'upper.adjust',
    'statistic.adjust' and 'pval.adjust'


## metasens, version 0.4-1 (2020-07-02)

### User-visible changes

* Use Markdown for NEWS
    
### Internal changes

* Call funnel.meta() instead of funnel.default() which will be removed
  from R package **meta**, version 4.13-0


## metasens, version 0.4-0 (2019-08-06)

### Major changes

* New function metamiss() implementing imputation methods for missing
  binary data ([Gamble & Hollis,
  2005](https://doi.org/10.1016/j.jclinepi.2004.09.013); [Higgins et
  al., 2008](https://doi.org/10.1097/XEB.0000000000000141))
  

* Use **roxygen2** for development of R package **metasens**

* Dataset 'nsaids' renamed to 'Moore1998'


## metasens, version 0.3-2 (2017-12-06)

### Major changes

* Version of R package **meta** must be larger or equal 4.9-0

* For limit meta-analysis, R_b measure of between-study heterogeneity
  added ([Crippa et al., 2016](https://doi.org/10.1002/sim.6980))

* P-values can be printed in scientific notation

* P-values equal to 0 are actually printed as "0" instead of
  "< 0.0001"

* Thousands separator can be used in printouts and forest plots for
  large numbers

### User-visible changes

* print.summary.limitmeta():
  - new argument print.Rb to specify if heterogeneity measure should
    be shown in output

### Internal changes

* Code customisations due to changes in R package **meta**, version
  4.9-0


## metasens, version 0.3-1 (2016-10-15)

### User-visible changes

* plot.copas():
  - default label on x-axis in two bottom plots should read
    '... largest se' instead of '... largest sd'


## metasens, version 0.3-0 (2016-02-16)

### Major changes

* Version of R package **meta** must be larger or equal 4.0-0

* Checks implemented which are available in R package **meta**

### User-visible changes

* copas():
  - new argument 'sign.rsb' (which has been available in function
    summary.copas since version 0.6-3 of R package **copas**)

* plot.copas(), print.copas(), and summary.copas():
  - consider default value for argument 'sign.rsb' from copas object

* limitmeta():
  - argument 'sm' removed (not necessary)
  - set z- and p-value of adjusted effect equal to NA for metaprop
    objects

* Help pages updated

### Bug fixes

* print.copas():
  - other values than 0.1 considered for argument 'sign.rsb'
  
* print.summary.limitmeta():
  - no error for metaprop objects

### Internal changes

* copas() and print.copas():
  - use internal function chklevel() from R package **meta** to check
    significance level


## metasens, version 0.2-0 (2014-12-06)

### Major changes
 
* Argument 'backtransf' added to R package **meta**, version 3.8-0,
  considered in R package **metasens**.

### User-visible changes

* copas(), limitmeta(), and funnel.limitmeta():
  - new argument 'backtransf'

* orbbound():
  - new argument 'backtransf'
  - new argument 'left' to chose whether selection bias is expected
    on the left or right side of the funnel plot

* forest.orbbound(), print.copas(), print.limitmeta(),
  print.orbbound(), print.summary.copas(),
  print.summary.limitmeta():
  - new argument 'backtransf' which replaces argument 'logscale'

* Several help pages updated


## metasens, version 0.1-0 (2014-06-24)

### First version released on CRAN - replacement of R package **copas**

### Major changes

* Limit meta-analysis implemented ([Rücker et al.,
  2011](https://doi.org/10.1093/biostatistics/kxq046))

### User-visible changes

* New functions:
  - limitmeta(), print.limitmeta()
  - summary.limitmeta(), print.summary.limitmeta()
  - funnel.limitmeta()

* print.copas() and print.summary.copas():
  - new argument 'logscale' to print results for relative effect
    measures on log scale

* New dataset 'nsaids'

* New help pages added and update of existing help pages

### Internal changes
 
* New internal function radialregression()
