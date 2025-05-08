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

We will also be using several other packages, including `ggplot2`,
`dplyr`, and `purrr`. You can install these packages using the following
command:

``` r
library(ggplot2)
library(dplyr)
```

## Age Classes

The next step is to establish age classes for the simulation. Since the
goal is to evaluate whether theoretical culling profiles allow herd
sizes to be be maintained given intrinsic mortality and the reproductive
biology of sheep and goats, the age classes and lengths (or widths) of
each age class are informed by studies focusing on culling profiles. For
example, Payne’s (1973) work separates ages into classes A through I and
Grant’s (1982) work distinguishes age by month. Synthesized mandible
ontogenic stages provided by Vigne and Helmer (2007) demonstrate theses
groupings and their corresponding ages in years. Although Payne’s (1973)
Age class grouping is the most common implementation Marom and Bar-Oz
(2009) compare 10 different the culling profiles by restructuring age
groups to nine age classes ranging from 0.17 to 8 years.

| Age.Class | months | years |
|:----------|-------:|------:|
| A (0-2m)  |      2 |  0.17 |
| B (2-6m)  |      4 |  0.33 |
| C (6-12m) |      6 |  0.50 |
| D (1-2y)  |     12 |  1.00 |
| E (2-3y)  |     12 |  1.00 |
| F (3-4y)  |     12 |  1.00 |
| G (4-6y)  |     24 |  2.00 |
| H (6-8y)  |     24 |  2.00 |
| I (\>8y)  |    Inf |   Inf |

Table S1.1. Age classes and widths

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub.

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Grant1982" class="csl-entry">

Grant, Annie. 1982. “<span class="nocase">The Use of Tooth Wear as a
Guide to the Age of Domestic Ungulates</span>.” In *Ageing and Sexing
Animal Bones from Archaeological Sites - British Archaeological Records
International Series 109.*, edited by Bob Wilson, Caroline Grigson, and
Sebastian Payne, 91–108. Oxford: Archaeopress.

</div>

<div id="ref-Marom2009" class="csl-entry">

Marom, Nimrod, and Guy Bar-Oz. 2009. “<span class="nocase">Culling
profiles: the indeterminacy of archaeozoological data to survivorship
curve modelling of sheep and goat herd maintenance strategies</span>.”
*Journal of Archaeological Science* 36 (5): 1184–87.
<https://doi.org/10.1016/j.jas.2009.01.007>.

</div>

<div id="ref-Payne1973" class="csl-entry">

Payne, Sebastian. 1973. “<span class="nocase">Kill-off Patterns in Sheep
and Goats: The Mandibles from Aşvan Kale</span>.” *Anatolian Studies* 23
(September): 281–303. <https://doi.org/10.2307/3642547>.

</div>

<div id="ref-Vigne2007" class="csl-entry">

Vigne, Jean-Denis, and Daniel Helmer. 2007. “<span class="nocase">Was
milk a ‘secondary product’ in the Old World Neolithisation process? Its
role in the domestication of cattle, sheep and goats</span>.”
*Anthropozoologica* 42 (2): 9–40.

</div>

</div>
