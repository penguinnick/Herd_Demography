To Cull or Not to Cull? Using Population Projection Modeling to Reveal
Risk-Sensitive Herd Culling Strategies
================
Nicholas Triozzi
2025-02-25

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Procedures for Herd Demography JCA Article

<!-- badges: start -->
<!-- badges: end -->

# Livestock Population Dynamics and Optimziation of Culling Rates

## Introduction

This document outlines the procedures used to simulate goat and sheep
herd dynamics under the constraints of various culling strategies. The
content of this document and associated scripts are supporting materials
for \[TITLE OF ARTICLE\].

## Install HerdDynamics Package

The first step is to install the `HerdDynamics` package. This package
contains all the functions and data used in this analysis.

``` r
devtools::install_github("penguinnick/HerdDynamics")
library(HerdDynamics)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub.
