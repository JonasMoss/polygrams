polygrams
================

The polygrams package includes function for handling of *Bernstein polygrams*, a class of non-parametric density estimators.

Installation
------------

From inside R, the package can be installed from GitHub by using devtools.

``` r
library("devtools")
install_github("JonasMoss/polygrams")
```

Example usage
-------------

Load the library as is usually done.

``` r
library("polygrams")
```

As an application, we study the monthly sunspot data.

``` r
dat = datasets::sunspot.month
summary(dat)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.00   15.70   42.00   51.96   76.40  253.80

Based on the summary we choose the support \[0,270\]. Let's try |s| = 10:

``` r
plot(polygram(dat, support = c(0, 270), s = 10))
```

<img src="README_files/figure-markdown_github/plot-1.png" width="750px" />

This density could be convex and monotone:

``` r
plot(polygram(dat, support = c(0, 270), s = 10, shape = "convex", monotone = "decreasing"))
```

<img src="README_files/figure-markdown_github/plot_monotone-1.png" width="750px" />
