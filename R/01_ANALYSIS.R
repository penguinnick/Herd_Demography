###########################################################################
#                                                                         #
#         Full Analysis Script for JCAA Article                           #
#                                                                         #   
# How to use:                                                             #                       
#         Execute code from parent folder                                 #
#             to reproduce results and figures in the article.            #  
#                                                                         #
###########################################################################
#-- change this to your local path to the project folder
# r_path = "~/Herd_Demography_UBELIX/" # HPC cluster
r_path = "./" # local

#-- load functions
source(paste0(r_path,"R/00_FUNS.R"))
# Load packages
using( "tidyverse", "purrr", "boot",  "stringr", "rcartocolor", "officer", "ggplot2", "ggh4x", "car", "parallel", "devtools", "sensitivity",
      "flextable", "parallelly", "pbapply", "mmage"
      )
# install.packages("devtools")
# devtools::install_github("penguinnick/HerdDynamics")
library(HerdDynamics)

library(sensitivity)
#-- load parameter values
# source("./00_PARAMETERS.R")
source(paste0(r_path,"R/00_PARAMETERS.R"))


#####################--SETUP--######################
cat("Building tcla table...\n")
#-- Build tcla table using HerdDynamics::build_tcla() 
nbphase = 1
tcla = build_tcla(
  female.ages = Payne_ages$lclass,
  male.ages = Payne_ages$lclass,
  nbphase = nbphase
)

#####################--MORTALITY DATA--######################
cat("Processing mortality data and building archaeological culling profiles...\n")
#-- read in age-at-death data
wear.df = read.table(file = paste0(r_path,"data/age_at_death_data.csv"), sep = ",")

#-- create list for each site/period
wear.df.list = wear.df %>%
  group_by(Site, Period) %>%
  group_split(.keep = T)

#-- get names
l.names = wear.df %>%
  group_by(Site, Period) %>%
  group_keys() %>%
  reframe(name = paste(Site, Period))

#-- assign names to each list
names(wear.df.list) = l.names$name

cprof <- lapply(wear.df.list, function(w){
  t <- correct.counts(w$Payne.Group, probability.correction = FALSE)
  t <- cbind.data.frame(t, s = survivorship(t$n))
  HI = 1 - (sqrt (1 -  t$s[8]))
  G = 1 - (sqrt(1 - t$s[7]))
  t$s[7] = t$s[8] = G
  t$s[9] = HI
  list(t = t, culling.profile =  c(100, t$s * 100))
})

cprof <- lapply(cprof, function(x){
  x$culling.profile})

culling.profiles <- cprof

#-- create Baseline offtake (no offtake, 100% survivorship)
Baseline.offtake = list(Baseline = rep(100, length(offtake_models$Energy)))

#-- put all culling strategies into a single list
all.offtake <-  c(offtake_models, culling.profiles, Baseline.offtake)

#-- convert survivorship to mortality probabilities
offtake.mortality = lapply(all.offtake, function(x) {
  1 - (x / 100)
})

#####################--PARAMETERS--######################
cat("Processing parameters for mortality, prolificacy, and parturition rates...\n")
#-- read in lambing and kidding rates
pro.dat = pro.dat
#-- summarize prolificacy rates from the literature
NetPro = pro.dat %>% 
  group_by(Taxon) %>% 
  reframe(mn.pro = mean(LitterSize),
          sd.pro = sd(LitterSize))

#-- read in parameters file for mortality, prolificacy, and parturition rates
param.dat = param.dat


#-- Ages, parturition interval, and age of first parturition set here
ages = HerdDynamics::Payne_ages$ages  # unique(param.dat$Age)
parturition.Interval = parturition.Interval # 300 days between births
part.age = part.age


set.seed(388)
#-- generate Prolificacy rates for Table 1 using function above
set_prolificacy("Goat",  1)
set_prolificacy("Sheep", 2)

#-- calculate Parturition Rates based on prolificacy for Table 1
param.dat$ParturitionRate <-  ARR(param.dat$Prolificacy, parturition.Interval)

#-- create parms lists
goat.parms = list(
  ages = ages,
  parturition = ARR(NetPro$mn.pro[1], parturition.Interval),
  # calculates annual reproduction rate
  parturition.Interval = parturition.Interval,
  part.age = part.age,
  MeanProlificacy = NetPro$mn.pro[1],
  sdProlificacy = NetPro$sd.pro[1],
  f.mortality = get_mortality("Goat", "Female"),
  m.mortality = get_mortality("Goat", "Male")
)

sheep.parms = list(
  ages = ages,
  parturition = ARR(NetPro$mn.pro[2], parturition.Interval),
  # calculates annual reproduction rate
  parturition.Interval = parturition.Interval,
  part.age = part.age,
  MeanProlificacy = NetPro$mn.pro[2],
  sdProlificacy = NetPro$sd.pro[2],
  f.mortality = get_mortality("Sheep", "Female"),
  m.mortality = get_mortality("Sheep", "Male")
)

sheep.param.props = list(
  tcla = tcla,
  parms = sheep.parms,
  nbphase = nbphase,
  female.offtake = female.offtake
)

goat.param.props = list(
  tcla = tcla,
  parms = goat.parms,
  nbphase = nbphase,
  female.offtake = female.offtake
)

param.props = list(goat = goat.param.props, sheep = sheep.param.props)

#####################--Initialize Transition Matrix--######################
cat("Initializing transition matrices for each strategy...\n")
#-- create param table for each offtake strategy with varying fertility and mortality
set.seed(542)
param = lapply(param.props, function(p) {
  p$parms = vary.fert.mort(p$parms)
  lapply(offtake.mortality, function(o) {
    with(p , {
      build_param(
        tcla,
        parms,
        nbphase,
        female.offtake,
        correctionfec = TRUE,
        offtake = o
      )
    })
  })
})

#-- create simulation environment
cat("Creating simulation environment for each strategy...\n")
set.seed(600)
system.time({
  listpar = lapply(param.props, function(p) {
    make.listpar(
      param.props = p,
      nbcycle = nbcycle,
      offtake.mortality = offtake.mortality
      )
    })
})

#-- save stochastic environment to replicate results
# save(listpar, file = "../data/listpar.RData")
# load(file = "../data/listpar.RData")

cat("Calculating lambda and sex proportions...\n")
# -- calculate lambda and sex proportions using wrapper.repro helper fun
# -- puts all lambda and female proportions in lists
system.time({
  lambda.list = wrapper.repro(listpar, out = "lambda")
  sex.prop.list = wrapper.repro(listpar, out = "sex")
  #-- assign names to each list
})
# 
# cat("Running Calculation of lambda and sex proportions on multiple cores...\n")
#-- this is the parallelized version of the above code, which can be used to speed up calculations
# system.time({
# #-- use parallel computing to speed up bootstrap calculations
# cl <- parallel::makeCluster(n.cores, type = "PSOCK") #system('nproc'))
# parallel::clusterExport(cl, c("listpar", "tcla", "p0", "get_lambda", "wrapper.repro_parallel"))
# parallel::clusterEvalQ(cl, 
#              library(HerdDynamics),
#              library(boot),
#              library(dplyr))
# lambda.list = parallel::parLapply(cl, seq_along(listpar), 
#                                   wrapper.repro_parallel, 
#                                   listpar = listpar, 
#                                   out = "lambda", 
#                                   tcla = tcla,
#                                   p0=p0)
# sex.prop.list = parallel::parLapply(cl, seq_along(listpar), 
#                                   wrapper.repro_parallel, 
#                                   listpar = listpar, 
#                                   out = "sex", 
#                                   tcla = tcla,
#                                   p0=p0)
# })
# parallel::stopCluster(cl)
#-- assign names to each list
names(lambda.list) <- names(sex.prop.list) <- names(listpar)
#-- unlist each list
lambda.list = unlist(lambda.list, recursive = F)
sex.prop.list = unlist(sex.prop.list, recursive = F)

#####################--BOOTSTRAP LAMBDA--######################
cat("Bootstrapping lambda and sex proportions...\n")
#-- bootstrap resample lambda and store in data.frame
lambda.boot.df = do.call(
  rbind.data.frame, 
  lapply(lambda.list, FUN = boot.fun, n = n.boot) # n.boot set to 1000 in parameters file
  )

#-- bootstrap sex proportions and store in data.frame
sex.prop.boot.df = do.call(
  rbind.data.frame, 
  lapply(sex.prop.list, function(l) {
    s = unlist(l)
    boot.fun(s, n.boot)}
    )
  )

#-- create taxon and strategy columns
taxon.strat = str_split(rownames(lambda.boot.df), "\\.", simplify = T)

#-- bind columns from each bootstrapped results table
repro.boot.df = cbind.data.frame(
  lambda.boot.df %>%
    reframe(
      Taxon = taxon.strat[, 1],
      strategy = taxon.strat[, 2],
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

#-- check correlation between lambda and sex proportions
cat("Testing correlation between lambda and proportion.../n")
print(
  repro.boot.df %>% 
    filter(strategy != "Baseline") %>% 
      with(cor.test(Lambda, Proportion.Female, method = "pearson"))
  )

#-- use function HerdDynamics::getLambda to extract reproduction traits
repro = lapply(listpar, function(l) {
  lapply(l, function(x) {
    HerdDynamics::get_lambda(x[[1]], tcla = tcla, p0 = p0)
  })
})

#-- age classes vector for plotting and projections
ageClasses = c("", Payne_ages$ageClasses)

#-- use wrapper function to get initial herd structure in a data.frame for projection
xini.df = xini.to.data.frame(repro, ageClasses = ageClasses, ages = ages)
cat("Initial herd structure for projections:\n")
print(head(xini.df))

#####################--RUN PROJECTIONS--######################
cat("Running projections for each strategy...\n")
#-- project herd for each strategy and store results in a list
listpar = unlist(listpar, recursive = F)

results = lapply(listpar, function(l) {
  projectHerd2(listpar = l, p0 = p0)
})

#-- summarize total population through time
tot.pop.res = lapply(results, function(r) {
  pop_summary2(r, sex = FALSE, interval = "year")
})

#-- aggregate into data.frame
tot.pop.df = res.to.df(tot.pop.res)

#-- run simulation for all strategies, replicate results, and store in an array
cat("Running stochastic simulations for each strategy and replicating results...\n")

set.seed(1056)
cl <- parallel::makeCluster(n.cores, type = "PSOCK") # system('nproc'))
parallel::clusterEvalQ(cl, {
  library(boot)
  library(HerdDynamics)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

#-- export necessary functions and variables to cluster nodes
parallel::clusterExport(cl, 
                        varlist = c(
                          "generate_AdultMortalityRates",
                          "vary.fert.mort",
                          "offtake.mortality", 
                          "stochastic.rep", 
                          "projectHerd2", 
                          "pop_summary2", 
                          "make.listpar", 
                          "nbcycle", 
                          "nbrep",
                          "run_stochastic_rep",
                          "p0",
                          "ARR"
                        )
)

stochastic.sim.res = replicate(n = nbrep,
                                lapply(param.props, function(p) {
                                  lapply(offtake.mortality, function(o) {
                                    r = stochastic.rep(
                                      p0 = p0,
                                      param.props = p,
                                      offtake.mortality = o,
                                      nbcycle = nbcycle
                                    )
                                  })
                                }),
                                simplify = "array")

#-- create matrix, then df from results
sto.res.mat = apply(stochastic.sim.res, 1, FUN = list.to.mat, simplify = F)
# sto.res.mat = apply(stochastic.sim.res.nested, 1, FUN = list.to.mat, simplify = F)
sto.res.df = lapply(sto.res.mat, FUN = matrix.to.df)
parallel::stopCluster(cl)
#-- set taxon fields
sto.res.df$goat$taxon = "goat"
sto.res.df$sheep$taxon = "sheep"

#-- combine into one df
sto.res.df = do.call(rbind.data.frame, sto.res.df)
# sto.res.df$strategy = factor(sto.res.df$strategy, levels = c(names(offtake.mortality)))

#####################--OPTIMIZE CULLING RATES--######################
cat("Optimizing culling rates for each strategy...\n")
#-- test of optimization function for all years in the Energy strategy for goats, as an example.
optimize.res = lapply(listpar$goat.Energy, function(b) {
  optimize(f = get.off.adjust.poff,
           param.ref = b,
           interval = c(0, 5))$minimum
})

#-- phi
optimize.res[1:5]

#-- set threshold as mean of all lambda
all.lambda = unlist(lambda.list[-c(which(str_detect(names(lambda.list), "Baseline")))])
Lambda.threshold = quantile(all.lambda)
Lambda.threshold.low = Lambda.threshold[[2]]
Lambda.threshold.high = Lambda.threshold[[4]]

#-- reproject with optimization
cat("Reprojecting with optimized culling rates...\n")
new.results = lapply(listpar, function(l) {
  l = lapply(
    l,
    FUN = adjust.offtake,
    low.threshold = Lambda.threshold.low,
    high.threshold = Lambda.threshold.high,
    p0 = p0
  )
  projectHerd2(listpar = l, p0 = p0)
})

#-- send results to df
new.tot.pop.res = lapply(new.results, function(r) {
  pop_summary2(r, sex = FALSE, interval = "year")
})

# res2df(new.tot.pop.res) # this function is in the HerdDynamics package
new.pop.res.df = res.to.df(new.tot.pop.res)

#-- create column for offtake adjustment
tot.pop.df$offtake = "unadjusted"
new.pop.res.df$offtake = "adjusted"

#-- merge into a single data.frame
all.res.df = rbind.data.frame(tot.pop.df, new.pop.res.df)

#-- summarizes mean population size for each strategy, used for Figure 8
new.vprod = do.call(rbind.data.frame, lapply(new.results, function(l) {
  fprod(formula = ~ 1, l$vecprod)
}))
new.vprod$Taxon = str_split_i(rownames(new.vprod), "\\.", i = 1)
new.vprod$strategy = str_split_i(rownames(new.vprod), "\\.", i = 2)
new.vprod = new.vprod %>%
  mutate(
    strategy = factor(strategy, levels = rev(names(offtake.mortality)), labels = rev(
      c(
        names(offtake.mortality)[1:5],
        "Meat A",
        "Meat B",
        "Milk A",
        "Milk B",
        names(offtake.mortality)[10:16]
      )
    )),
    Taxon = factor(Taxon, levels = c("goat", "sheep"))
  )

#####################--LEVENE'S TEST--######################
#-- summarize multiplication rate for old and new results
m.df = lapply(results, function(r) {
  o = mmage.fm(formula = ~ cycle , vecprod = r$vecprod)
  return(data.frame(cycle = o$cycle, m = o$m))
})
#-- create df
m.df = res.to.df(m.df)

#-- new res
new.m.df = lapply(new.results, function(r) {
  o = mmage.fm(formula = ~ cycle , vecprod = r$vecprod)
  return(data.frame(cycle = o$cycle, m = o$m))
})
#-- create df
new.m.df = res.to.df(new.m.df)

#-- create column for offtake adjustment
m.df$offtake = "unadjusted"
new.m.df$offtake = "adjusted"

#-- merge into a single data.frame
all.m.df = rbind.data.frame(m.df, new.m.df)

#-- Levene's test for homogeneity of variance in multiplication rate between adjusted and unadjusted offtake
all.m.list = all.m.df %>%
  mutate(strategy = factor(
    strategy,
    levels = names(offtake.mortality),
    labels = c(
      names(offtake.mortality)[1:5],
      "Meat A",
      "Meat B",
      "Milk A",
      "Milk B",
      names(offtake.mortality)[10:16]
    )
  )) %>%
  group_by(strategy, Taxon) %>%
  group_split()

m.lev.res = lapply(all.m.list, function(l) {
   lev.res = l %>% with(car::leveneTest(m, offtake))
  data.frame(
    strategy = first(l$strategy),
    Taxon = first(l$Taxon),
    F.value = lev.res$`F value`[1],
    p = lev.res$`Pr(>F)`[1]
  )
})

lev.res.df = do.call(rbind.data.frame, m.lev.res)
#-- add lambda values from bootstrap results
lev.res.df$lambda[lev.res.df$Taxon == "goat"] = lambda.boot.df$t0[1:16]
lev.res.df$lambda[lev.res.df$Taxon == "sheep"] = lambda.boot.df$t0[17:32]

#####################--SENSITIVITY ANALYSIS--######################
#-- WARNING -- high computational demand. 
# source(paste0(r_path,"R/XX_SENSITIVITY.R")) 

#-- Sensitivity Analysis was run separately for each strategy and taxon, and results were saved as .rds files in the output folder. The code below reads in the sobol indices for each strategy and combines them into a single data.frame for plotting and analysis.
# sobol.files = list.files("./output/", pattern = "sobol_indices_*", full.names = T)[-1]
sobol.files = list.files("./output/", pattern = "*res_goat*|*res_sheep*", full.names = T)
sobol_res = lapply(sobol.files, readRDS)

si <- lapply(sobol_res, function(s) {
  s = unlist(s, recursive = F)
  do.call(rbind.data.frame, lapply(s, extract_indices)) %>%  
    mutate(
      res = str_split_i(rownames(.), "\\.", 2),
      taxon = str_split_i(rownames(.), "\\.", 3),
      strategy = str_split_i(rownames(.), "\\.", 4)
    ) 
})

sobol_indices = do.call(rbind.data.frame, si)


###########################################################################
#                                                                         #
#                               outputs                                   #
#                                                                         #                 
###########################################################################
if(!dir.exists(paste0(r_path,"output"))){
  dir.create(paste0(r_path,"output"))
}
outpath = paste0(r_path,"output/")
saveRDS(sobol_indices, file = paste0("output/","sobol_indices_df.rds"))
save(stochastic.sim.res, file = paste0(outpath,"stochastic_sim_res.RData"))
save(listpar, file =paste0(outpath, "listpar.RData"))
write.table(repro.boot.df, file = paste0(outpath,"reproboot.csv"), sep = ",")

#-- Figures and Tables
source(paste0(r_path,"R/02_PLOTS_TABLES.R"))
