# metasens: Advanced Statistical Methods to Model and Adjust for Bias in Meta-Analysis
Official GitHub repository of R package **metasens**

[![Build Status](https://travis-ci.org/guido-s/metasens.svg?branch=master)](https://travis-ci.org/guido-s/metasens)
[![CRAN Version](http://www.r-pkg.org/badges/version/metasens)](https://cran.r-project.org/package=metasens)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/metasens)](http://cranlogs.r-pkg.org/badges/metasens)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/metasens)](http://cranlogs.r-pkg.org/badges/grand-total/metasens)
[![Research software impact](http://depsy.org/api/package/cran/metasens/badge.svg)](http://depsy.org/package/r/metasens)


## Description

The following methods are implemented to evaluate how sensitive the results of a meta-analysis are to potential bias in meta-analysis and to support Schwarzer et al. (2015), Chapter 5 "Small-Study Effects in Meta-Analysis":
 - Copas selection model (Copas and Shi, 2001; Schwarzer et al., 2010);
 - limit meta-analysis (Rücker et al., 2011);
 - upper bound for outcome reporting bias (Copas and Jackson, 2004).

Furthermore, R package **metasens** provides functions and datasets to
support Schwarzer et al. (2015), Chapter 5 "Small-Study Effects in
Meta-Analysis", http://meta-analysis-with-r.org/ .

### References

[Copas J, Jackson D (2004), A bound for publication bias based on the fraction of unpublished studies.  *Biometrics*, **60**, 146-153.](https://scholar.google.de/scholar?q=Copas+Jackson+2004+A+bound+for+publication+bias+based+on+the+fraction+of+unpublished+studies)

[Copas JB, Shi JQ (2001), A sensitivity analysis for publication bias in systematic reviews.  *Statistical Methods in Medical Research*, **10**, 251-265.](https://scholar.google.de/scholar?q=Copas+Shi+2001+A+sensitivity+analysis+for+publication+bias+in+systematic+reviews)

[Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011), Treatment-effect estimates adjusted for small-study effects via a limit meta-analysis. *Biostatistics*, **12**, 122-142.](https://scholar.google.de/scholar?q=Rücker+Schwarzer+Carpenter+Binder+Schumacher+2011+Treatment-effect+estimates+adjusted+for+small-study+effects+via+a+limit+meta-analysis)

[Schwarzer G, Carpenter J, Rücker G (2010), Empirical evaluation suggests Copas selection model preferable to trim-and-fill method for selection bias in meta-analysis.  *Journal of Clinical Epidemiology*, **63**, 282-288.](https://scholar.google.de/scholar?q=Schwarzer+Carpenter+Rücker+2010+Empirical+evaluation+suggests+Copas+selection+model+preferable+to+trim-and-fill+method+for+selection+bias+in+meta-analysis)

[Schwarzer G, Carpenter JR and Rücker G (2015), *Meta-Analysis with R (Use-R!)*. Springer International Publishing, Switzerland](http://www.springer.com/gp/book/9783319214153)


## Installation

### Current official [![CRAN Version](http://www.r-pkg.org/badges/version/metasens)](https://cran.r-project.org/package=metasens) release:
```r
install.packages("metasens")
```

### Current beta / GitHub release:

Installation using R package
[**ghit**](https://cran.r-project.org/package=ghit) (without
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows):
```r
install.packages("ghit")
ghit::install_github("guido-s/metasens")
```

Installation using R package
[**devtools**](https://cran.r-project.org/package=devtools) (with
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows):
```r
install.packages("devtools")
devtools::install_github("guido-s/metasens")
```


### Bug Reports:

```r
bug.report(package = "metasens")
```
