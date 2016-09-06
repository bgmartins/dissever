[![Travis-CI Build Status](https://travis-ci.org/pierreroudier/dissever.svg?branch=master)](https://travis-ci.org/pierreroudier/dissever)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/dissever)](http://cran.r-project.org/web/packages/dissever)

# dissever
Dissever: an R package implementing methods for spatial downscaling/dissagregation

`dissever` is a general method for spatial downscaling introduced by Malone *et al.* in their [2012 *Computers & Geosciences* paper](http://www.sciencedirect.com/science/article/pii/S0098300411002895).

An R package from Pierre Roudier implemented dissever in R, modifying the general procedure so that, through the  [the `caret` package](https://topepo.github.io/caret), numerous regression methods could be tested. The code in this repository further extends the R package from Pierre Roudier, in order to support spatial dissagregation and geographically weighted regression.

The package is not on CRAN yet, but can be installed easily using the `devtools` package. If you don't have `devtools` installed:

```
install.packages('devtools')
```

Then to install `dissever`:

```
devtools::install_github('bgmartins/dissever')
```


