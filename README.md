To Cull or Not to Cull? Using Population Projection Modeling to Reveal
Risk-Sensitive Herd Culling Strategies in Neolithic Dalmatia
================
Nicholas Triozzi
2025-02-25

# Livestock Population Dynamics and Optimziation of Culling Rates

## Introduction

This document outlines the procedures used to simulate goat and sheep
herd dynamics under the constraints of various culling strategies. The
content of this document and associated scripts are supporting materials
for “To Cull or Not to Cull? Using Population Projection Modeling to
Reveal Risk-Sensitive Herd Culling Strategies in Neolithic Dalmatia”.

## Overview

The population projection model presented here uses the
[HerdDynamics](https://github.com/penguinnick/HerdDynamics) package for
R, which features wrapper functions for the [mmage
package](https://gitlab.cirad.fr/selmet/livtools/mmage) ([Lesnoff
2024a](#ref-Lesnoff2024)) and customized scripts for comparing different
herd culling strategies ([Triozzi 2024](#ref-Triozzi2024)). In the
`HerdDynamics` package, some procedures described in the [mmage
reference
manual](https://gitlab.cirad.fr/selmet/livtools/mmage/-/blob/master/doc/mmage_v1.7.pdf?ref_type=heads)
have been consolidated into wrapper functions such as “build_tcla”, and
“build.param”.

The `HerdDynamics` package is designed to be used with the `mmage`
package, which is a powerful tool for analyzing age-structured
populations. The `HerdDynamics` package provides additional
functionality for simulating herd dynamics and optimizing culling
strategies.

## Setup - Population Projection Matrix

Here, a Lefkovitch population projection model ([Lefkovitch
1965](#ref-Lefkovitch1965)) is used to project sheep and goat herd
population dynamics. The probability of survival and fecundity are
influenced by the competing risks of intrinsic mortality and offtake
(i.e., slaugther). The transition matrix is created using published data
on fecundity, intrinsic mortality, and offtake rates which are derived
from theoretical culling profiles in archaeological literature (See
Chapter Four).

# Age Classes

The first step is to establish age classes for the simulation. Since the
goal is to evaluate whether theoretical culling profiles allow herd
sizes to be be maintained given intrinsic mortality and the reproductive
biology of sheep and goats, the age classes and lengths (or widths) of
each age class are informed by studies focusing on culling profiles. For
example, Payne’s ([1973](#ref-Payne1973)) work separates ages into
classes A through I and Grant’s ([1982](#ref-Grant1982)) work
distinguishes age by month. Synthesized mandible ontogenic stages
provided by Vigne and Helmer ([2007](#ref-Vigne2007)) demonstrate theses
groupings and their corresponding ages in years. Although Payne’s
([1973](#ref-Payne1973)) Age class grouping is the most common
implementation Marom and Bar-Oz ([2009](#ref-Marom2009)) compare 10
different the culling profiles by restructuring age groups to nine age
classes ranging from 0.17 to 8 years.

    ##   ageClasses  ages lclass
    ## 1   A (0-2m)  0.17      2
    ## 2   B (2-6m)  0.50      4
    ## 3  C (6-12m)  1.00      6
    ## 4   D (1-2y)  2.00     12
    ## 5   E (2-3y)  3.00     12
    ## 6   F (3-4y)  4.00     12
    ## 7   G (4-6y)  6.00     24
    ## 8   H (6-8y)  8.00     24
    ## 9    I (>8y) 10.00    Inf

    ##   Age.Class months years
    ## 1  A (0-2m)      2  0.17
    ## 2  B (2-6m)      4  0.33
    ## 3 C (6-12m)      6  0.50
    ## 4  D (1-2y)     12  1.00
    ## 5  E (2-3y)     12  1.00
    ## 6  F (3-4y)     12  1.00
    ## 7  G (4-6y)     24  2.00
    ## 8  H (6-8y)     24  2.00
    ## 9   I (>8y)    Inf   Inf

| ageClasses |  ages | lclass |
|:-----------|------:|-------:|
| A (0-2m)   |  0.17 |      2 |
| B (2-6m)   |  0.50 |      4 |
| C (6-12m)  |  1.00 |      6 |
| D (1-2y)   |  2.00 |     12 |
| E (2-3y)   |  3.00 |     12 |
| F (3-4y)   |  4.00 |     12 |
| G (4-6y)   |  6.00 |     24 |
| H (6-8y)   |  8.00 |     24 |
| I (\>8y)   | 10.00 |    Inf |

Table S1.1. Age classes and widths

The model presented here is an “untruncated” model since the final age
class includes animals 8 years old and up, with no terminal age
specified for old animals. The age classes according to Payne
([1973](#ref-Payne1973)) are included in `HerdDynamics` as `Payne_ages`.
Note the last value in `lclass` is set to `Inf`.

``` r
HerdDynamics::Payne_ages
```

    ##   ageClasses  ages lclass
    ## 1   A (0-2m)  0.17      2
    ## 2   B (2-6m)  0.50      4
    ## 3  C (6-12m)  1.00      6
    ## 4   D (1-2y)  2.00     12
    ## 5   E (2-3y)  3.00     12
    ## 6   F (3-4y)  4.00     12
    ## 7   G (4-6y)  6.00     24
    ## 8   H (6-8y)  8.00     24
    ## 9    I (>8y) 10.00    Inf

``` r
# # Payne's Class
# ageClasses = c("", "A (0-2m)",  "B (2-6m)", "C (6-12m)", "D (1-2y)", "E (2-3y)", "F (3-4y)", "G (4-6y)", "H (6-8y)", "I (>8y)")
# 
# #-- Age classes in years
# # ages = c(  0.17,  0.5,   1, 2, 3, 4, 5, 6, 7)
# ages = c(  0.17,  0.5,   1, 2, 3, 4, 6, 8, 10)
# 
# #-- length of each age class
# # lclass = c(0.17, 0.33, 0.5, 1, 1,   1, 1, 1, Inf ) # in years
# # lclass = c( 2,   4,    6, 12, 12, 12, 12, 12, Inf) # in months
# lclass =  c( 2,      4,    6,  12, 12, 12, 24, 24, Inf)
```

Once the age classes are set, the wrapper function `build_tcla()` from
`HerdDynamics` uses the `mmage::fclass` function to create the
data.frame `tcla`. Parameter `nbphase` sets the timestep interval, which
should be either annually (12) or monthly (1). `tcla` table includes for
each age group, the length of each group (i.e., 2 months, 4 months, 6
months, 12 months, etc.), and the minimum and maximum ages of the group.
This procedure is streamlined by the customized `build_tcla` function
which is designed to adjust the values provided as input (lclass)
according to `nbphase`, and recognize whether `lclass` is formatted for
a “truncated” or “untruncated” model. We use the `lclass` variable in
the `Payne_ages` table provided by `HerdDynamics`, which has already
been set up to track age in months. So we can set `nbphase` to 1.

``` r
#-- Build tcla table using HerdDynamics::build_tcla function which uses mmage package. Here we use the same male and female age classes. The parameter nbphase is set to 1 since lclass is already formatted as months.
nbphase = 1
tcla = build_tcla(female.ages = Payne_ages$lclass, male.ages = Payne_ages$lclass, nbphase = nbphase)
flextable(tcla) %>% 
  theme_vanilla() %>% 
  autofit() %>%
  set_caption( "Table S1.2. Initiated age and sex table (tcla), containing age classes, lengths, and minimum and maximum ages in months for each age class." )
```

<img src="README_files/figure-gfm/tcla-1.png" width="564" />

## Herd Culling Strategies

`HerdDynamics` includes survivorship probabilities associated with
Redding’s ([1981](#ref-Redding1981)) **Energy** and **Security**,
Payne’s ([1973](#ref-Payne1973)) **Meat**, **Milk**, and **Wool**
strategies, and Vigne and Helmer’s ([2007](#ref-Vigne2007)) Meat and
Milk types A and B (**MeatA**, **MeatB**, **MilkA**, **MilkB**) and
**Fleece** strategies. Marom and Bar-Oz ([2009](#ref-Marom2009)) provide
survivorship probabilities that were standardized across nine age
classes. These data are stored as a list in the HerdDynamics package as
`offtake_models`.

The offtake rates given here are generated from zooarchaeological
research which does not account for sex-based differences in culling
practices. But these offtake rates are often interpreted on the level of
the herd regardless of sex. Here they are used to model the offtake
rates of males. These rates are reduced for females in a later step.

<img src="README_files/figure-gfm/theoretical-survivorship-table-1.png" width="1661" />

    ## [1] "C:/Users/nt25j496/OneDrive - Universitaet Bern/Writing_Projects/Herd_Demography"

## Age-at-death data from four Neolithic sites

Age-at-death data for sheep and goat mandibles collected from
**Benkovac-Barice**, **Islam Grčki**, **Smilčić**, and **Zemunik Donji**
are used to create survival probabilities for each age class ([Triozzi
2024](#ref-Triozzi2024)). First, the data is called from
`age_at_death_data.csv`. The data has already been formatted and
filtered, combining all goats and sheep as ovicaprids, and dropping Late
Neolithic Islam Grcki, for which there are only 10 mandibles. The
absolute frequencies are calculated using the `HerdDynamics` function,
`correct_counts()` (see Appendix B) that accounts for mandibles that
were assigned multiple age classes. Then the survivorship probabilities
are calculated following Price et al. ([2016](#ref-Price2016)), using
the function `survivorship()` from `HerdDynamics` package. Lastly, the
data is formatted as a list, similar to the `offtake.models` list
created previously.

<img src="README_files/figure-gfm/Table-3-site-survivorship-1.png" width="1287" />

### Plot survivorship curves

The data in Tables C.2 and C.3 can now be visualized as survivorship
curves.

![Figure S1.1. Survivorship curves for threoretical harvest
strategies.](README_files/figure-gfm/Figure-S1-1-suvivorship-curves-1.png)![Figure
S1.1. Survivorship curves for threoretical harvest
strategies.](README_files/figure-gfm/Figure-S1-1-suvivorship-curves-2.png)![Figure
S1.1. Survivorship curves for threoretical harvest
strategies.](README_files/figure-gfm/Figure-S1-1-suvivorship-curves-3.png)

## Herd Population Growth Parameters

Offtake represents the number of animals slaughtered in each age group.
These numbers will affect herd size changes by restricting the number of
animals available to reproduce. Two other important parameters are
*intrinsic mortality* and *fertility*.

### Mortality

The probability that an animal will survive from one time step to the
next is affected by the competing risks of being slaughtered (i.e.,
offtake) and intrinsic mortality. We set intrinsic mortality for males
and females separately. According to Redding ([1981](#ref-Redding1981)),
there is no conclusive evidence that there is a difference in mortality
between male and female lambs and kids. Redding applies a mortality rate
of 32% for lambs and 45% for kids in the first year of life. These rates
can be proportionally allocated to the first three age classes, which
include animals under 1 year. Redding ([1981](#ref-Redding1981)) assigns
ewe and doe mortality from ages 1-9 to 18% per year; 10% for rams and
15% for bucks aged 1 to 2 years, 5% for both species ages 2-6 years and
100% for ages six and up. Male mortality is modeled slightly lower than
female mortality since rams and bucks selected for breeding should be
more robust and have better resistance to disease, whereas this
artificial selection should be weaker for females. We simulate this via
higher mortality for older females.

In this program age classes for males and females may be specified
separately. However, culling profiles constructed from archaeological
remains cannot distinguish male from female mandibles. Therefore,
offtake rates modeled here will be applied to the entire herd while
intrinsic mortality rates will be defined separately for males and
females.

### Fertility - Reproduction Parameters

Several parameters are important regarding the reproductive biology of
goats and sheep. `pfbirth` specifies probability of giving birth to a
female. This value is set as 0.5 by default in the code generating the
parameter table. `part.age` specifies the age of first parturition. , as
an index corresponding to the age classes defined by `lclassf` and
`lclassm`. By setting this variable to 4, we set the age of first
parturition to 1 year. `parturition` specifies the number of
parturitions per female per year. Redding ([1981](#ref-Redding1981))
models a breeding rate of once per year but notes that under conditions
of very good pasturage goats may breed twice per year (1981:59).
Parturition rate is also referred to as the annual reproductive rate
(ARR) and is calculated as: $$ARR=\frac{l*365}{i}$$ where `l` is the
litter size (i.e., prolificacy rate) multiplied by the number of days in
the year and `i` is the parturition interval ([Upton
1984](#ref-Upton1984ModelsOI); [R. T. Wilson, Peacock, and Sayers
1984](#ref-Wilson1984)). To keep variation in fertility consistent,
parturition varies from year to year, depending on the prolificacy
rates, which are generated as described below. `parturitionInterval` is
set to 300 days, reflecting a mean of 10 months for bot goats and sheep
([R. T. Wilson 1989](#ref-Wilson1989)).

`prolificacy` specifies the prolificacy rate, defined as the number of
live offspring per parturition per year. Goats are capable of giving
birth to twins, triplets, and quadruplets, while kidding rates increase
with age up to five or six years (Redding, 1981:70). To simulate
stochastic inter-annual variation in fertility, a function is used to
sample prolificacy rates from a normal distribution given mean kidding
and lambing rates obtained from the literature. The `prolificacy.csv`
file contains this data.

The net fecundity rate is the product of the parturition and net
prolificacy rates. The parturition rate is calculated differently for
birth-pulse (all births occur at the same time) and birth-flow models
(births occur continuously over the time interval (t, t+1)). The
Parturition rate `vpar` is dependent on the number of individuals
removed from the herd, or `ptot` which is the sum of natural death rates
and offtake rates: $$vpar=\frac{1-ptot}{2}*hpar$$ where `hpar` is the
number of parturitions expected per female-time unit.

Using the mortality and prolificacy data called in the previous chunk
two lists are created (one for goats, one for sheep) containing the
fertility and mortality parameters used to compute **lambda**
($\lambda$), reproductive values, and project herd growth. The function
`generate_prolificacy_rates()` is used to vary inter-annual fertility
values. This function and `generate_infant_mortality_rates` are included
in HerdDynamics package.

The above function is used in this chunk to generate a series of
prolificacy rates based on fertility data. Also defined is the function
`Part.rate` which calculates $ARR$.

``` r
# param.dat$Age = ages

#-- function for auto-generating prolificacy rates, used in build param function when prolificacyRate="auto"
NetPro = pro.dat %>% group_by(Taxon) %>% reframe(mn.pro = mean(LitterSize), sd.pro = sd(LitterSize))
NetPro
```

    ## # A tibble: 2 × 3
    ##   Taxon mn.pro sd.pro
    ##   <chr>  <dbl>  <dbl>
    ## 1 Goats   1.49  0.275
    ## 2 Sheep   1.15  0.122

``` r
#-- Ages, parturition rate, and age of first parturition set here
ages = unique(param.dat$Age)

#-- function to calculate annual reproduction rate (ARR), included in HerdDynamics package 
ARR = function(mean_litter_size, parturition_interval){
  (mean_litter_size * 365) / parturition_interval
}

ARR(mean_litter_size = 1.3, parturition_interval = 300)
```

    ## [1] 1.581667

``` r
#-- age of first parturition set to animals in the 1-2 year age class
parturition.Interval = 300 
part.age = 2

#-- get vector of age structured mortality rates for male and female goats
goat.mort.f = param.dat %>% filter(Taxon == "Goat" & Sex == "Female") %>% select(Mortality)
goat.mort.m = param.dat %>% filter(Taxon == "Goat" & Sex == "Male") %>% select(Mortality)
#-- get vector of age structured mortality rates for male and female sheep
sheep.mort.f = param.dat %>% filter(Taxon == "Sheep" & Sex == "Female") %>% select(Mortality)
sheep.mort.m = param.dat %>% filter(Taxon == "Sheep" & Sex == "Male") %>% select(Mortality)

#-- create parms lists
goat.parms = list(
  ages = ages,
  parturition = ARR(NetPro$mn.pro[1], parturition.Interval),  # calculates annual reproduction rate
  parturition.Interval = 300, 
  part.age = part.age,
  MeanProlificacy = NetPro$mn.pro[1], 
  sdProlificacy = NetPro$sd.pro[1], 
  f.mortality = as.vector(goat.mort.f$Mortality),
  m.mortality = as.vector(goat.mort.m$Mortality)
)

sheep.parms = list(
  ages = ages,
  parturition = ARR(NetPro$mn.pro[2], parturition.Interval), # calculates annual reproduction rate
  parturition.Interval = 300, 
  part.age = part.age,
  MeanProlificacy = NetPro$mn.pro[2], 
  sdProlificacy = NetPro$sd.pro[2], 
  f.mortality = as.vector(sheep.mort.f$Mortality),
  m.mortality = as.vector(sheep.mort.m$Mortality)
)

#-- clean up
rm( goat.mort.m, goat.mort.f, sheep.mort.m, sheep.mort.f, pro.dat)
```

## Build Lefkovitch matrix

Now that the mortality, offtake, and fertility parameters have been
defined, it’s time to create the Lefkovitch matrix. The first step is to
convert the survivorship percentages to mortality probabilities. This is
done by dividing each value by 100 and subtracting from 1. A
**Baseline** strategy is also created to track herd population dynamics
free from the offtake constraint by setting survivorship probability to
100% which translates to an offtake rate of 0 for each age group.

In the chunk below, `Baseline.offtake` is created as a list of 100%
survivorship for all age classes to model herd growth with no offtake.
The list is appended to the list of all offtake models to be run in the
subsequent projections.

``` r
#-- put all harvest strategies into a single list
Baseline.offtake = list(rep(100, length(offtake_models$Energy)))
names(Baseline.offtake) = "Baseline"
all.offtake = append(offtake_models, culling.profiles)
all.offtake = append(offtake_models, culling.profiles2)
all.offtake = append(all.offtake, Baseline.offtake)

#-- convert survivorship to mortality
offtake.mortality = lapply(all.offtake, function(x){1-(x/100)})
```

Next, a list variable containing all the parameters defined above, the
original `tcla` table, `nbphase`, and a female offtake modifier,
`female.offtake` is created and named `param.props`. This list is
supplied to the `build.param` function which streamlines the creation of
the hazards table and Lefkovitch matrix used by mmage to project
population changes. The variable `female.offtake` is supplied (or not,
if set to `NULL`). This value adjusts the female offtake rate, setting
the female offtake probabilities to a product of offtake, female
mortality, and `female.offtake`. Setting this value to 0 implies that
females exit the herd acording to female mortality probabilities. Since
those values are global, setting `female.offtake` to 0 eliminates all
variation between different offtake models with respect to herd growth
and defeats the purpose of this experiment. `female.offtake` is
therefore set to `15` reflecting a culling rate of females due to
infertility reported by Malher et al. ([2001](#ref-Malher2001)).

``` r
#-- This number is used to account for culling of females due to infertility. If Null, no female offtake assumed.
female.offtake = 15

sheep.param.props = list(
  tcla = tcla, parms = sheep.parms, nbphase = nbphase, female.offtake = female.offtake)

goat.param.props = list(
  tcla = tcla, parms = goat.parms, nbphase = nbphase, female.offtake = female.offtake)

param.props = list(goat = goat.param.props, sheep = sheep.param.props)
```

This chunk shows the definition of the `vary.fert.mort()` function,
which is used to generate a new set of prolificacy and fertility rates
based on mean and standard deviations provided in `parms`.

``` r
set.seed(123)
#-- function to vary prolificacy and mortality
# source("../R/varyFertMort.R")
vary.fert.mort = function( parms, n = 6 ){
  parms$prolificacy = generate_prolificacy_rates(
    meanPro = parms$MeanProlificacy,
    sdPro = parms$sdProlificacy,
    n = n
    )
  #-- vary parturition rate
  parms$parturition = ARR(mean(parms$prolificacy), parms$parturition.Interval)
  
  
  # generate_InfantMortalityRates <- function( Mort = f.mortality, n = 3) {
  #   inf.mort = Mort[1:3]
  #   sort( generate_ProlificacyRates( meanPro = mean( inf.mort ), sdPro = sd( inf.mort ), n = 3), decreasing = TRUE )
    # sort( generate_ProlificacyRates( meanPro = sum( inf.mort ), sdPro = sd( inf.mort ), n = 3), decreasing = TRUE )
  # }
  
  generate_AdultMortalityRates <- function( Mort = f.mortality, n = 6) {
    adult.mort = Mort[4:length(Mort)]
    sort( generate_prolificacy_rates( meanPro = mean( adult.mort ), sdPro = sd( adult.mort ), n = 6), decreasing = FALSE )
    # sort( generate_ProlificacyRates( meanPro = sum( inf.mort ), sdPro = sd( inf.mort ), n = 3), decreasing = TRUE )
  }
  
  #-- update male and female mortality rates
  parms$f.mortality = c(
    generate_infant_mortality_rates(parms$f.mortality),
    generate_AdultMortalityRates(parms$f.mortality)
    )
  
  parms$m.mortality = c(
    generate_infant_mortality_rates(parms$m.mortality),
    generate_AdultMortalityRates(parms$m.mortality)
    )
  return(parms)
}

#-- example of parameters created for a single year.
vary.fert.mort(parms = goat.parms)
```

    ## $ages
    ## [1]  0.17  0.50  1.00  2.00  3.00  4.00  6.00  8.00 10.00
    ## 
    ## $parturition
    ## [1] 2.033403
    ## 
    ## $parturition.Interval
    ## [1] 300
    ## 
    ## $part.age
    ## [1] 2
    ## 
    ## $MeanProlificacy
    ## [1] 1.492941
    ## 
    ## $sdProlificacy
    ## [1] 0.2750401
    ## 
    ## $f.mortality
    ## [1] 0.28836290 0.20936609 0.12055848 0.06515195 0.11205936 0.11550642 0.20617116
    ## [8] 0.40912531 0.41804890
    ## 
    ## $m.mortality
    ## [1] 0.2413435 0.1466634 0.1264610 0.0398888 0.1226288 0.1461968 0.1596290
    ## [8] 0.2064455 0.2366741
    ## 
    ## $prolificacy
    ## [1] 1.320242 1.612634 1.869362 1.984414 1.837791 1.403297

Here a single param table is created for each survivorship profile in
**two steps**.

**First**, the parms object is updated with prolificacy and infant
mortality rates generated by sampling from a normal distribution using
the `vary.fert.mort()` function defined above.

**Second**, build_param function is run which creates a table that
contains offtake, mortality, fecundity, and survival probabilities for
each age class

``` r
#-- create param table for each offtake strategy with varying fertility and mortality
set.seed(542)
param = lapply(param.props, function( p ) {
  p$parms = vary.fert.mort( p$parms )
  lapply( offtake.mortality, function( o ) {
    with( p , {
      build_param( tcla, parms, nbphase, female.offtake,  correctionfec = TRUE, offtake = o)
      })
    })
  })

str(param$goat$Energy$param)
```

    ## 'data.frame':    20 obs. of  11 variables:
    ##  $ sex    : chr  "F" "F" "F" "F" ...
    ##  $ class  : int  0 1 2 3 4 5 6 7 8 9 ...
    ##  $ lclass : num  1 2 4 6 12 ...
    ##  $ cellmin: num  0 1 3 7 13 25 37 49 73 97 ...
    ##  $ cellmax: num  0 2 6 12 24 ...
    ##  $ nupar  : num  0 0 0 0.147 1.807 ...
    ##  $ ff     : num  0 0 0 0.101 1.241 ...
    ##  $ fm     : num  0 0 0 0.101 1.241 ...
    ##  $ pdea   : num  0.0836 0.1494 0.1163 0.097 0.0146 ...
    ##  $ poff   : num  0.00687 0.01319 0.01567 0.03483 0.07632 ...
    ##  $ g      : num  1 0.4558 0.1997 0.1137 0.0468 ...

``` r
rm(sheep.param.props, goat.param.props, Baseline.offtake)
```

# Part II - Simulating Environmental Variation

In the previous chunk, the function `build.param()` was run to produce a
parameter set using auto-generated prolificacy rates and infant
mortality rates using mean and standard deviation for each taxon.

In this chunk, environmental variation is simulated by creating a new
`param` table for each time-step containing unique sets of fertility and
mortality parameters using the `vary.fert.mort()` function described
above. This is done using the `make.listpar()` function which has the
`vary.fert.mort()` function built into it. The function simulates a
single environment with unique fertility and mortality rates while
keeping the offtake rates defined earlier constant. The environment is
called `listpar`. The number of years over which the simulation will run
is specified by `nbcycle`.

``` r
source("./R/makeListpar.R")
# 
#-- makeListpar definition
make.listpar = function( param.props, nbcycle, offtake.mortality  ){
  # source("../R/buildparam.R")
  # source("../R/varyFertMort.R")
  with( param.props, {

    #-- get number of timesteps
    nbstep = nbcycle * nbphase

    #-- create list of parms
    parms.rep = replicate(n = nbstep, vary.fert.mort( parms = parms ), simplify = FALSE)

    #-- create listpar
    lapply(offtake.mortality, function(o){
      lapply(parms.rep, function(p){
        build_param(
          tcla = tcla,
          parms = p, offtake = o,
          female.offtake = female.offtake,
          correctionfec = TRUE,
          nbphase = nbphase
        )$param
      })
    })
  })
}
```

# Create Stochastic Environment

Using the functions defined above (‘vary.fert.mort()’ and
‘make.listpar()’), a stochastic environment is created for each offtake
strategy. The function `make.listpar()` is run with the parameters
defined above. The number of cycles is set to `200` and the initial
population size is set to `150`.

``` r
set.seed(600)
#-- set number of cycles and initial population size
nbcycle = 200
p0 = 150

#-- makeListpar definition
make.listpar = function( param.props, nbcycle, offtake.mortality  ){
  # source("../R/buildparam.R")
  # source("../R/varyFertMort.R")
  with( param.props, {

    #-- get number of timesteps
    nbstep = nbcycle * nbphase

    #-- create list of parms
    parms.rep = replicate(n = nbstep, vary.fert.mort( parms = parms ), simplify = FALSE)

    #-- create listpar
    lapply(offtake.mortality, function(o){
      lapply(parms.rep, function(p){
        build_param(
          tcla = tcla,
          parms = p, offtake = o,
          female.offtake = female.offtake,
          correctionfec = TRUE,
          nbphase = nbphase
        )$param
      })
    })
  })
}
#-- uncomment lines below to re-create simulation environment or use existing by loading from data folder
listpar = lapply(param.props, function(p){
  make.listpar(
    param.props = p, 
    nbcycle = nbcycle, 
    offtake.mortality = offtake.mortality
    )
  })

#-- save stochastic environment to replicate results
save(listpar, file = "data/listpar.RData")
load(file = "data/listpar.RData")
```

## Basic Herd structure and growth

With the fertility and mortality rates have been established the next
step is to extract the structure of the herd, lambda, and the proportion
of male to females for each model. To do this, an initial herd size,
`p0` is set to `150`. In the following chunk, these traits are
calculated for comparison of the impacts of culling strategies on herd
demography.

While lambda can be computed from a single `param` table, the herd
multiplication rate will change from one time-step to the next as
fertility and mortality rates change. To estimate lambda given this
variation a bootstrapping procedure is used.

Using the list of mortality and fertility parameters created in the
previous chunk, the custom function `getLambda()` is run using `lapply`
to obtain lambda, sex proportion, and number of individuals in each age
class for the initial population. This function is a wrapper function
that uses functions in the `mmage` package to obtain the dominant
eigenvalue (lambda) to get the multiplication rate for the population.

With the environment simulated as `listpar`, the chunk below calculates
lambda and proportion of females for each year under each offtake
strategy.

``` r
#-- a function that gets lambda, sex proportion, and herd structure.
# source("../R/getLambda.R")

#-- wrapper function to get demography information from every listpar table
wrapper.repro = function( listpar, out = c("lambda", "sex") ){
  out <- match.arg(out)
  lapply(seq_along(listpar), function(l) {
  p = listpar[[l]]
  lapply(p, function(s){
    sapply(s, function(x){
      r = get_lambda(x, tcla=tcla, p0=p0) #$lambda
      if(out=="lambda"){
        return(r$lambda)
      } else {
        if(out=="sex"){
          return(r$sex.proportion[1,2])
        }
      }
      })
    })
  })
}

#-- put all lambda and female proportions in lists
lambda.list = wrapper.repro(listpar, out = "lambda")
sex.prop.list = wrapper.repro(listpar, out = "sex")

names(lambda.list) = names(listpar)
names(sex.prop.list) = names(listpar)

lambda.list = unlist(lambda.list, recursive = F)
sex.prop.list = unlist(sex.prop.list, recursive = F)
```

## Bootstrap lambda and sex proportions

This chunk defines the function used in the bootstrap estimation of mean
$\lambda$ ($\lambda_{boot}$) and proportion of females. The bootstrap is
replicated 1000 times. The output is `repro.boot.df`.

``` r
set.seed(665)
#-- Define function to calculate the mean
mean_function <- function(data, indices) {
  # This function will be applied to resampled data
  return(mean(data[indices]))
}

#-- function to compute bootstrapped mean and 95% confidence intervals 
boot.fun = function( s, n ){
  # Perform bootstrapping
  bootstrap_results = boot(data = s, statistic = mean_function, R = n )
  # Obtain bootstrapped confidence interval
  boot_conf_interval = boot.ci(bootstrap_results, type = "perc")
  #-- returns mean, lower and upper ci
  return(
    data.frame(
      t0 = boot_conf_interval$t0, 
      low = boot_conf_interval$percent[4], 
      up = boot_conf_interval$percent[5]
      )
    )
}

#-- set number of bootstrap replicates
n.boot = 1000

lambda.boot.df = do.call(rbind.data.frame, lapply(lambda.list, FUN = boot.fun, n = n.boot))
sex.prop.boot.df = do.call(rbind.data.frame, 
                           lapply(sex.prop.list, function(l){ 
                             s=unlist(l); boot.fun(s, n.boot) 
                             }))

#-- create taxon and strategy columns
taxon.strat = str_split(rownames(lambda.boot.df), "\\.", simplify = T)

#-- bind columns from each bootstrapped results table
repro.boot.df = cbind.data.frame(
  lambda.boot.df %>% 
    reframe(
      Taxon = taxon.strat[,1], 
      strategy = taxon.strat[,2], 
      Lambda = t0, 
      low.lambda = low, 
      up.lambda = up
      ),
  sex.prop.boot.df %>% 
    reframe(
      "Proportion.Female" = t0, 
      "Proportion.Male" = 1 - t0, 
      low.sex.F = low, 
      up.sex.F = up
      )
  )

head(repro.boot.df)
```

    ##   Taxon strategy   Lambda low.lambda up.lambda Proportion.Female
    ## 1  goat   Energy 1.006795  0.9999808 1.0136126         0.6119602
    ## 2  goat Security 0.999645  0.9930901 1.0062981         0.6256652
    ## 3  goat     Meat 1.002505  0.9958053 1.0089477         0.6667756
    ## 4  goat     Milk 0.962240  0.9565801 0.9681252         0.7859124
    ## 5  goat     Wool 1.009340  1.0030086 1.0159178         0.6657699
    ## 6  goat    MeatA 0.933090  0.9273255 0.9381032         0.7397474
    ##   Proportion.Male low.sex.F  up.sex.F
    ## 1       0.3880398 0.6031908 0.6200410
    ## 2       0.3743348 0.6182360 0.6322519
    ## 3       0.3332244 0.6592271 0.6741369
    ## 4       0.2140876 0.7821001 0.7896977
    ## 5       0.3342301 0.6592174 0.6718739
    ## 6       0.2602526 0.7346593 0.7447408

``` r
write.table(repro.boot.df, file = "output/reproboot.csv", sep = ",")

# rm(lambda.list, sex.prop.list, sex.prop.boot.df)
```

``` r
#-- Plot bootstrapped lambda results
Fig3 = repro.boot.df %>%
  mutate(strategy = factor(strategy, levels = repro.boot.df$strategy[1:16])) %>%
  ggplot(aes(x = strategy, y = Lambda, color = Taxon)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = cbPalette[10] ) +
  geom_point(size = 2, position = position_dodge(width = 0.5), alpha = 0.8 ) +
  # geom_point(aes(size = Proportion.Female), position = position_dodge(width = 0.5), alpha = 0.8 ) +
  # geom_point(aes(size = `Proportion Female`), position = position_dodge(width = 0.5), alpha = 0.8 ) +
  geom_errorbar(aes(ymin = low.lambda, ymax = up.lambda), width = 0.5, position = "dodge") + 
  scale_color_manual(values = cbPalette[c(1,2)]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Fig3
```

![](README_files/figure-gfm/lambda-plot-1.png)<!-- -->

``` r
#-- save figure
ggsave(Fig3, filename = "./Figures/Fig3_lambda-boot.jpg", dpi = 300, width = 6, height = 4)
```

In this chunk, $\lambda_{boot}$ and proportion of females in the herd
are plotted to show the negative correlation between these two values.

``` r
#-- check correlation between lambda and sex proportions
repro.boot.df %>% filter(strategy!="Baseline") %>% with(cor.test( Lambda, `Proportion Female` ))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  Lambda and Proportion Female
    ## t = -8.9905, df = 28, p-value = 9.55e-10
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.9325402 -0.7273791
    ## sample estimates:
    ##        cor 
    ## -0.8618092

## Herd Structure

This chunk calculates herd structure which is used to produce Figure 4.

``` r
#-- set initial population size
p0 = 150

#-- call function getLambda to extract reproduction traits
repro = lapply(listpar, function(l){
  lapply(l, function(x){
    get_lambda(x[[1]], tcla=tcla, p0=p0) 
  })
})

#-- wrapper function to get a dataframe with herd structure from repro object
ageClasses = c("", Payne_ages$ageClasses)
xini.to.data.frame = function( repro, ageClasses, ages ) {
  xini.list = lapply(repro, function( r ) {  # given repro, a list of reproductive params for n species,
    lapply( r, function( x ) {
      with(x, {
        # summarize initial herds (total males and females for each age class)
        xini = initial.herd %>%    
          group_by( class ) %>% 
          reframe( n = as.integer( sum( xini ) ) )
        }) 
      })
    })
  xini.list = unlist( xini.list, recursive = F ) # unlist species groupings
  taxon.strat = str_split( names( xini.list ), "\\." ) # create list of species and strategies
  # create Taxon, strategy, and AgeClass fields for each df in xini.list
  for ( i in 1: length( xini.list ) ) {
    ts = taxon.strat[[ i ]]
    taxon = rep( unlist( ts )[ 1 ], 9 )
    strat = rep( unlist( ts )[ 2 ], 9 )
    xini.list[[ i ]]$Taxon = taxon
    xini.list[[ i ]]$strategy = strat
    xini.list[[ i ]]$AgeClass = ageClasses[ -1 ]
    xini.list[[ i ]]$Age = as.character(ages)
  }
  # call all dataframes into a single df
  xini.df = do.call( rbind.data.frame, xini.list )
  # make strategy field a factor
  xini.df$strategy = factor( xini.df$strategy, levels = c(  names(offtake.mortality) ))
  xini.df$AgeClass = factor( xini.df$AgeClass, levels = ageClasses[-1] )
  return( xini.df )
}

xini.df = xini.to.data.frame(repro, ageClasses = ageClasses, ages = ages)
head(xini.df)
```

    ## # A tibble: 6 × 6
    ##   class     n Taxon strategy AgeClass  Age  
    ##   <int> <int> <chr> <fct>    <fct>     <chr>
    ## 1     1    49 goat  Energy   A (0-2m)  0.17 
    ## 2     2    56 goat  Energy   B (2-6m)  0.5  
    ## 3     3    31 goat  Energy   C (6-12m) 1    
    ## 4     4    10 goat  Energy   D (1-2y)  2    
    ## 5     5     1 goat  Energy   E (2-3y)  3    
    ## 6     6     0 goat  Energy   F (3-4y)  4

# Part III - Population Projection

The next step is to simulate herd demographic change through time given
the parameters defined above and the reproductive traits calculated for
each model and culling profile. The initial herd size was established
already as `p0` and the simulated environment as `listpar`. The number
of times-teps is equal to the length of `listpar`.

The projection is run using the script “projectHerd2.R”. This function
takes fewer inputs than the `projectHerd()` function. First, `listpar`
is created using the makeListpar.R function. This function simulates
environmental variation of one herd under inter-annual variation in
mortality and fertility. One herd, multiple offtake strategies can be
assessed this way. Results are in the same format as projectHerd.R. This
function is a wrapper function that streamlines the herd projection
procedure described in the mmage documentation ([Lesnoff
2024b](#ref-mmage)) but is used when the environment is predefined. The
alternative function `projectHerd` has additional functionality such as
varying environments from one step to the next. This is useful when
there is only one `param` table but environmental variation is generated
for each time-step. This process is explained in another document.

The output of the projection is saves as `results`. When the projection
is run for every strategy using the mortality and fertility rates for
each taxon using `lapply`, the output is a list of length n.strategies
\* n.taxa.

``` r
#-- load projectHerd function
source("./R/projectHerd2.R")

listpar = unlist(listpar, recursive = F)
results = lapply(listpar, function(l) { projectHerd2( listpar = l, p0 = p0) } )
# results = lapply(listpar, function(l) { project_herd( listpar = l, p0 = p0) } )
```

## Output Metrics

Here summary statistics of interest are calculated to compare how the
different models perform (i.e., values that will be used to compare
differences in herd size, etc.).

First, a wrapper function is defined as `pop.summary`, which takes the
results list (`results`) and creates a list of data.frames each
containing two columns: the time step (`time`) and the number of animals
in the herd (`pop`). Then another wrapper function, `res.to.df` is used
to create a single data.frame of all results which is used for plotting.

``` r
#-- function to summarize dynamics of the total population through time
pop_summary2 = function(x, sex = FALSE, interval = "year"){
  vecx=x$vecx
  if(interval=="month"){
    if( sex ) {
      out = aggregate( x ~ tim + sex , data=vecx, FUN=sum)
      colnames(out) = c("time", "sex", "pop")
    } else {
      out = aggregate(x ~ tim, data=vecx, FUN=sum)
      colnames(out) = c("time", "pop")
    }
  }
  if(interval=="year"){
    if( sex ) {
      out = aggregate( x ~ cycle + sex , data=vecx, FUN=sum)
      colnames(out) = c("time", "sex", "pop")
    } else {
      out = aggregate(x ~ cycle , data=vecx, FUN=sum)
      colnames(out) = c("time", "pop")
    }
  }
  return(out)
}

tot.pop.res = lapply(results, function(r){pop_summary2(r, sex = FALSE, interval = "year")})
# tot.pop.res.sex = lapply(results, function(r){pop.summary(r, sex = TRUE, interval = "year")})

# function to gather results into a single dataframe
res.to.df = function(tot.pop.res){
  df = do.call(rbind.data.frame, args=c(tot.pop.res, make.row.names=FALSE))
  strats = names(tot.pop.res)
  df$strategy = unlist(lapply(strats, FUN=rep, (nrow(df)/length(strats))))
  ts = str_split(df$strategy, "\\.", simplify = T)
  df$Taxon = ts[, 1]
  df$strategy = ts[, 2]
  return(df)
}
tot.pop.df = res.to.df(tot.pop.res)
# tot.pop.df.sex = res.to.df(tot.pop.res.sex)
head(tot.pop.df)
```

    ##   time       pop strategy Taxon
    ## 1    1 150.00000   Energy  goat
    ## 2    2 147.29744   Energy  goat
    ## 3    3 122.53690   Energy  goat
    ## 4    4  96.93422   Energy  goat
    ## 5    5  93.12667   Energy  goat
    ## 6    6  80.79782   Energy  goat

``` r
# head(tot.pop.df.sex)
```

<figure>
<img src="README_files/figure-gfm/Fig5-1.png"
alt="Figure. 5. Plots showing the change in sheep and goat herd sizes through time under the various survivorship profiles used." />
<figcaption aria-hidden="true">Figure. 5. Plots showing the change in
sheep and goat herd sizes through time under the various survivorship
profiles used.</figcaption>
</figure>

## Stochastic variation

Since there is an element of stochastic variation, it is useful to
examine multiple iterations of the simulation. Here `nbrep` is set as
the number of simulations to replicate. A wrapper function is created,
`stochastic.rep` which combines `make.listpar`, `projectHerd2`, and
`pop.summary` functions. Then, the simulation is replicated `nbrep`
times. Results are combined into a data.frame for plotting.

``` r
set.seed(1056)
#-- function to replicate
stochastic.rep = function( p0, param.props, offtake.mortality, nbcycle ){
  #-- create new set of environmental parameters
  lp = make.listpar(
    param.props = param.props, 
    nbcycle = nbcycle, 
    offtake.mortality = list(offtake.mortality)
    )
  #-- reproject herd growth
  result = projectHerd2( listpar = unlist(lp, recursive = F) , p0 = p0 )
  #-- summarize results
  pop_summary2(result, sex = FALSE, interval = "year")$pop
}

#-- set number
nbrep = 10
nbcycle = nbcycle

#-- run simulation for all strategies
stochastic.sim.res = replicate(
  n = nbrep, 
  lapply(param.props, function( p ){
    lapply(offtake.mortality, function( o ){
      r = stochastic.rep(p0 = p0, param.props = p, offtake.mortality = o, nbcycle = nbcycle)
      })
    }), 
  simplify = "array"
  )

#-- create df from matrix output of replicated stochastic.rep func
matrix.to.df = function( mat ){
  df = as.data.frame(mat)
  df$time = as.numeric(str_replace_all(rownames(df), "[:alpha:]|[:punct:]", ""))
  df$strategy = str_replace_all(rownames(df), "[:digit:]", "")
  df %>% pivot_longer( cols = c(-time, -strategy) )
}

#-- function to transform lists into matrix
list.to.mat = function( l ){ 
  n.col = length( l )
  l.mat = matrix( unlist( l ), ncol = n.col)
  rownames(l.mat) = unique(names(unlist(l)))
  l.mat
}

#-- create matrix, then df from results
sto.res.mat = apply(stochastic.sim.res, 1, FUN = list.to.mat, simplify = F)
sto.res.df = lapply(sto.res.mat, FUN = matrix.to.df)

#-- set taxon fields
sto.res.df$goat$taxon = "goat"
sto.res.df$sheep$taxon = "sheep"

#-- combine into one df
sto.res.df = do.call(rbind.data.frame, sto.res.df)
sto.res.df$strategy = factor(sto.res.df$strategy, levels = c(names(offtake.mortality)))
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

<figure>
<img src="README_files/figure-gfm/Figure7-1.png"
alt="Figure. 7. Plots showing results of replicated projections of change in sheep and goat herd sizes through time under the various survivorship profiles used." />
<figcaption aria-hidden="true">Figure. 7. Plots showing results of
replicated projections of change in sheep and goat herd sizes through
time under the various survivorship profiles used.</figcaption>
</figure>

# Part IV - Optimizing offtake

Modification of offtake rates may allow for herd size to remain stable.
An optimization routine is needed to determine the factor, $\varphi$, by
which female offtake needs to be adjusted to maintain herd size. The
goal is to achieve a steady state multiplication rate under stochastic
changes in fertility and mortality. The first step is to determine the
multiplication rate of the herd from one timestep to the next. This is
done using the modified `mmage::fm` function, `mmagefm()`. In this
chunk, a dataframe, `lambda.df` is created by obtaining the
multiplication rate from the results of the herd projection which used
the `listpar` object that simulates inter-annual variation in fertility
and mortality.

After creating `lambda.df`, an index of the years with lowest lambda is
created. This index will be used to determine `phi` ($\varphi$) in a
later step.

## Adjust offtake function

In this chunk, a function is defined that will take as input a `param`
table from the list of params associated with the worst years for herd
growth and return a value of $(\lambda - m)^2$, which should be close to
0. This function will be used in the optimization routine in the
subsequent chunk.

``` r
#-- function to adjust offtake based on poff in param table
get.off.adjust.poff = function( phi, param.ref, m = 1 ){
  u = param.ref
  zf = u$poff[ u$sex == "F" ]
  zm = u$poff[ u$sex == "M" ]
  u$poff = c( phi * zf, zm )
  Lf = length( u$sex[ u$sex == "F" ]) - 1
  Lm = length( u$sex[ u$sex == "M" ]) - 1
  param = u
  A = mmage::fmat( param, Lf, Lm )$A
  ( mmage::feig( A )$lambda - m ) ^ 2
}
```

## Optimize

In this chunk, the function defined above is optimized to find
$\varphi_{opt}$, which is the value that will be used to adjust female
offtake in a new projection.

``` r
#-- optimize to find the value that minimizes (lambda-m)^2 (i.e., difference between projected herd growth rate and objective growth)
#-- an example of optimized phi
optimize(f=get.off.adjust.poff, param.ref = listpar$goat.Energy[[1]], interval = c(0,5))$minimum
```

    ## [1] 0.6648578

``` r
#-- optimize for all
optimize.res = lapply( listpar$goat.Energy, function( b ){
  optimize( f = get.off.adjust.poff, param.ref = b, interval = c(0,5))$minimum
  })

#-- phi
optimize.res[1:5]
```

    ## [[1]]
    ## [1] 0.6648578
    ## 
    ## [[2]]
    ## [1] 0.4286676
    ## 
    ## [[3]]
    ## [1] 7.802868e-05
    ## 
    ## [[4]]
    ## [1] 0.332951
    ## 
    ## [[5]]
    ## [1] 7.802868e-05

Now that the optimization procedure has been demonstrated for one set of
parameters, a function is defined that will adjust offtake rates by
$\varphi_{opt}$ when \$\< \$ `low.threshold` or \$\$ `high.threshold`.

### Dynamic optimization function

Here a function is defined that will multiply female offtake
probabilities by $\varphi_{opt}$ if $\lambda$ falls outside of the
defined thresholds. This will run inside of the project herd wrapper
function and return a new `param` table for every cycle

``` r
#-- function for adjusting p. female offtake
adjust.offtake = function( in.param, low.threshold, high.threshold, p0 = 150 ){
  #-- calculate phi.opt
  phi.opt = optimize(f = get.off.adjust.poff, param.ref = in.param, interval = c(0,5))$minimum
  
  #-- get unadjusted female offtake
  f.off = in.param$poff[ in.param$sex == "F" & in.param$class > 0 ]
  
  #-- extract tcla
  tcla1 = in.param[ , c( 1:5 ) ]
  
  #-- extract lclass
  # lclass1 = in.param$lclass[ in.param$sex=="F" & in.param$class > 0 ]
  
  #-- get lambda
  lambda = HerdDynamics::get_lambda( param = in.param,  tcla = tcla1, p0 = p0)$lambda
  
  #-- compare lambda to threshold value
  if( lambda < low.threshold | lambda >= high.threshold ){
    #-- if lambda below threshold, adjust female offtake
    in.param$poff[in.param$sex=="F" & in.param$class > 0] = f.off * phi.opt
  } else {
    in.param$poff[in.param$sex=="F" & in.param$class > 0] = f.off
  }
  return(in.param)
}
```

### Reproject herd with optimal offtake adjustments

Now the simulation is run after setting the upper and lower thresholds
as the 25% and 75% quantiles of $\lambda_{boot}$ stored in in
`repro.boot.df`.

Based on $\lambda_{boot}$ estimates, the optimization function will
adjust offtake if $\lambda$ is below 0.948 or above 1.016. After
optimizing, the original projection results are combined with the new
results and can be plotted.

``` r
#-- reproject with optimization
new.results = lapply( listpar, function( l ) {
  l = lapply( l, 
             FUN = adjust.offtake, 
             low.threshold = Lambda.threshold.low, 
             high.threshold = Lambda.threshold.high )
  projectHerd2( listpar = l, p0 = p0 ) 
  })

#-- send results to df
new.tot.pop.res = lapply(new.results, function(r){pop_summary2(r, sex = FALSE, interval = "year")})
# new.tot.pop.res = lapply(new.results, function(r){HerdDynamics::pop_summary(r, interval = "year")})
# res2df(new.tot.pop.res) #-- this function is in the HerdDynamics package
new.pop.res.df = res.to.df(new.tot.pop.res)

#-- create column for offtake adjustment 
tot.pop.df$offtake = "unadjusted"
new.pop.res.df$offtake = "adjusted"

#-- merge into a single data.frame
all.res.df = rbind.data.frame(tot.pop.df, new.pop.res.df)
```

![Figure C.5. Projected sheep and goat populations after applying
optimized offtake procedure for the Energy and Benkovac-Barice
strategies.](README_files/figure-gfm/plot-optimized-results-1.png)![Figure
C.5. Projected sheep and goat populations after applying optimized
offtake procedure for the Energy and Benkovac-Barice
strategies.](README_files/figure-gfm/plot-optimized-results-2.png)

# Part V - Evaluate optimization results

This chunk uses a modified version of `mmage::fprod()` function from the
`mmage` package ([Lesnoff 2024b](#ref-mmage)) to calculate mean herd
size as shown in Figure 9.

``` r
source("./R/mmagefprod.R")
new.vprod = do.call(rbind.data.frame, lapply(new.results, function(l){ fprod(formula = ~ 1, l$vecprod)  }))
new.vprod$Taxon = str_split_i(rownames(new.vprod), "\\.", i = 1)
new.vprod$strategy = str_split_i(rownames(new.vprod), "\\.", i = 2)
new.vprod = new.vprod %>%
  mutate(
    strategy = factor( strategy, 
      levels = rev(names(offtake.mortality)),
      labels = rev(c(
        names(offtake.mortality)[1:5], 
        "Meat A", "Meat B", "Milk A", "Milk B", 
        names(offtake.mortality)[10:16])
        )),
    Taxon = factor(Taxon, levels = c( "goat", "sheep")))
```

<figure>
<img src="README_files/figure-gfm/xmean-1.png"
alt="Figure 9. Mean herd size after optimized culling rates applied" />
<figcaption aria-hidden="true">Figure 9. Mean herd size after optimized
culling rates applied</figcaption>
</figure>

This chunk uses a modified version of the `fm()` function in the `mmage`
package to obtain the annual multiplication rates of the optimized
offtake (adjusted) and unadjusted simulations.

``` r
#-- modified mmage function to get multiplication rate of the herd
source("R/mmagefm.R")

#-- summarize multiplication rate for old and new results
m.df = lapply( results, function( r ){
  o = mmage.fm( formula = ~ cycle , vecprod = r$vecprod ); return( data.frame( cycle = o$cycle, m = o$m ))
  })
m.df = res.to.df(m.df)

new.m.df = lapply(new.results, function(r){
  o = mmage.fm( formula = ~ cycle , vecprod = r$vecprod ); return( data.frame( cycle = o$cycle, m = o$m ))
  })
new.m.df = res.to.df(new.m.df)

#-- create column for offtake adjustment 
m.df$offtake = "unadjusted"
new.m.df$offtake = "adjusted"

#-- merge into a single data.frame
all.m.df = rbind.data.frame(m.df, new.m.df)
```

10-year moving averages of the multiplication rates were calculated to
simplify illustration. This chunk defines functions to calculate the
moving averages and detect detect maximum and minimum multiplication
rates. These functions are used to produce Figures 5.7 and 5.8 in
Chapter 5.

``` r
# Function to calculate moving average
calculate_moving_average <- function(time_series, window_size, type = c("simple", "exponential"), alpha = 0.1) {
  type <- match.arg(type)
  
  if (type == "simple") {
    # Calculate simple moving average
    moving_average <- zoo::rollmean(time_series, k = window_size, fill = NA, align = "right")
  } else if (type == "exponential") {
    # Calculate exponential moving average
    moving_average <- zoo::rollapply(time_series, width = window_size, FUN = function(x) {
      n <- length(x)
      weights <- (1 - alpha)^(n:1)
      sum(weights * x) / sum(weights)
    }, fill = NA, align = "right")
  }
  
  return(moving_average)
}


# Function to detect peaks
detect_peaks <- function(x, span = 3) {
  z <- embed(x, span)
  s <- span %/% 2
  v <- max.col(z, ties.method = "first") == 1 + s
  c(rep(FALSE, s), v, rep(FALSE, s))
}

# Function to detect valleys
detect_valleys <- function(x, span = 3) {
  z <- embed(x, span)
  s <- span %/% 2
  v <- apply(z, 1, which.min) == 1 + s
  c(rep(FALSE, s), v, rep(FALSE, s))
}
```

### Levene’s Test for Equality of Variance

This section performs a Levene’s Test for Equality of variances between
the herd multiplication rate after offtake was optimized and the
predicted $\lambda$ value (i.e., unadjusted rates).

``` r
all.m.list = all.m.df %>% 
  mutate(strategy = factor(strategy, levels = names(offtake.mortality),
                           labels = c(names(offtake.mortality)[1:5], "Meat A", "Meat B", "Milk A", "Milk B",
                                      names(offtake.mortality)[10:16]))) %>% 
  group_by(strategy, Taxon) %>% 
  group_split()

m.lev.res = lapply(all.m.list, function(l){
  lev.res = l %>% with(car::leveneTest(m, offtake))
  data.frame(strategy = first(l$strategy), Taxon = first(l$Taxon), F.value = lev.res$`F value`[1], p = lev.res$`Pr(>F)`[1])
})

lev.res.df = do.call(rbind.data.frame, m.lev.res)
head(lev.res.df)
```

    ##   strategy Taxon   F.value           p
    ## 1   Energy  goat  5.127998 0.024080117
    ## 2   Energy sheep  6.287879 0.012553969
    ## 3 Security  goat  5.501712 0.019488071
    ## 4 Security sheep  9.022214 0.002835679
    ## 5     Meat  goat  4.640999 0.031816422
    ## 6     Meat sheep 10.979832 0.001005430

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Grant1982" class="csl-entry">

Grant, Annie. 1982. “<span class="nocase">The Use of Tooth Wear as a
Guide to the Age of Domestic Ungulates</span>.” In *Ageing and Sexing
Animal Bones from Archaeological Sites - British Archaeological Records
International Series 109.*, edited by Bob Wilson, Caroline Grigson, and
Sebastian Payne, 91–108. Oxford: Archaeopress.

</div>

<div id="ref-Lefkovitch1965" class="csl-entry">

Lefkovitch, L. P. 1965. “<span class="nocase">The study of population
growth in organisms grouped by stages</span>.” *Biometrics* 21 (1):
1–18.

</div>

<div id="ref-mmage" class="csl-entry">

Lesnoff, Matthieu. 2024b. *Mmage: A r Package for Sex-and-Age Population
Matrix Models*.

</div>

<div id="ref-Lesnoff2024" class="csl-entry">

———. 2024a. *<span class="nocase">mmage: A R package for sex-and-age
population matrix models</span>*.

</div>

<div id="ref-Malher2001" class="csl-entry">

Malher, X., H. Seegers, and F. Beaudeau. 2001.
“<span class="nocase">Culling and mortality in large dairy goat herds
managed under intensive conditions in western France</span>.” *Livestock
Production Science* 71 (1): 75–86.
<https://doi.org/10.1016/S0301-6226(01)00242-1>.

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

<div id="ref-Price2016" class="csl-entry">

Price, Max, Jesse Wolfhagen, and Erik Otárola-Castillo. 2016.
“<span class="nocase">Confidence Intervals in the Analysis of Mortality
and Survivorship Curves in Zooarchaeology</span>.” *American Antiquity*
81 (1): 157–73. <https://doi.org/10.7183/0002-7316.81.1.157>.

</div>

<div id="ref-Redding1981" class="csl-entry">

Redding, Richard William. 1981. “<span class="nocase">Decision making in
subsistence herding of sheep and goats in the Middle East</span>.”
Doctoral Dissertation, University of Michigan.

</div>

<div id="ref-Triozzi2024" class="csl-entry">

Triozzi, Nicholas Peter. 2024. “<span class="nocase">Herding is Risky
Business: Agropastoral Strategies in Neolithic Northern
Dalmatia</span>.” Ph.D. Dissertation, University of California, Santa
Barbara.

</div>

<div id="ref-Upton1984ModelsOI" class="csl-entry">

Upton, Martin. 1984. “<span class="nocase">Models of improved production
systems for small ruminants</span>.” In *Proceedings of the Workshop on
Small Ruminant Production Systems in the Humid Zone of West Africa*,
edited by J. E. Sumberg and K. Cassaday, 55–67. Ibadan, Nigeria:
International Livestock Centre for Africa.
<https://api.semanticscholar.org/CorpusID:131461953>.

</div>

<div id="ref-Vigne2007" class="csl-entry">

Vigne, Jean-Denis, and Daniel Helmer. 2007. “<span class="nocase">Was
milk a ‘secondary product’ in the Old World Neolithisation process? Its
role in the domestication of cattle, sheep and goats</span>.”
*Anthropozoologica* 42 (2): 9–40.

</div>

<div id="ref-Wilson1989" class="csl-entry">

Wilson, R T. 1989. “<span class="nocase">Reproductive performance of
African indigenous small ruminants under various management systems: a
review</span>.” *Animal Reproduction Science* 20 (4): 265–86.

</div>

<div id="ref-Wilson1984" class="csl-entry">

Wilson, R. T., Christie Peacock, and A. R. Sayers. 1984.
“<span class="nocase">Aspects of reproduction in goats and sheep in
south-central Kenya</span>.” *Animal Science* 38 (3): 463–67.
<https://doi.org/10.1017/S0003356100041660>.

</div>

</div>
