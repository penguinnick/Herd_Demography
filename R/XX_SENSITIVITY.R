###########################################################################
#                                                                         #
#         Sensitivity analysis only                                       #
#                                                                         #   
# How to use:                                                             #                       
#         Execute code in HPC or limit calculations by                    #
#         changing parameters                                             #
#                                                                         #  
#                                                                         #
###########################################################################
 

#####################--LIBRARIES--######################
#-- uncomment to install necessary packages on HPC environment
# Add ggrepel and sensitivity together
# pkgs <- c("sensitivity", "ggrepel", "ggplot2", "dplyr", "parallel")

# lib_dir <- "~/R/x86_64-pc-linux-gnu-library/4.4"
# cran_repo <- "https://stat.ethz.ch/CRAN/"

# dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
# .libPaths(c(normalizePath(lib_dir), .libPaths()))

# to_install <- pkgs[!pkgs %in% installed.packages(lib = lib_dir)[,"Package"]]

# if (length(to_install)) {
#   cat("Installing", paste(to_install, collapse = ", "), "with dependencies...\n")
  
#   install.packages(to_install,
#                    lib = lib_dir,
#                    repos = cran_repo,
#                    dependencies = TRUE,
#                    quiet = TRUE,
#                    Ncpus = 8)
# }
# Verify and load
# lapply(pkgs, library, character.only = TRUE)

library(sensitivity)
library(parallel)
library(HerdDynamics)
library(dplyr)
library(mmage)
library(boot)
library(tidyr)
library(stringr)
library(pbapply)
library(parallelly)

#####################--SENSITIVITY ANALYSIS PARAMETERS--######################

np = 500 # number of parameter sets to run for sensitivity analysis
nboot.sens = 50
Lambda.threshold.low = 0.936
Lambda.threshold.high = 1.005
p0 = 150
female.offtake = 15

avail <- availableConnections()
free <- freeConnections()
n_cores <- min(120, free) 
n.cores = n_cores
# n.cores = detectCores() # use all cores except one for parallel processing 

#####################--ESSENTIAL FUNCTIONS--######################
#-- change this to your local path to the project folder
# r_path = "~/Herd_Demography_UBELIX/" # HPC directory
r_path = "./" # local directory

#-- load functions
source(paste0(r_path,"R/00_FUNS.R"))
# source("R/00_FUNS.R")

#####################--ESSENTIAL DATA--######################
load("data/listpar.RData")
listpar = unlist(listpar, recursive = FALSE)
listpar = listpar[c(3:5, 19:21)] # subset to just Paynes models for sensitivity analysis
#####################--SENSITIVITY ANALYSIS--######################
cat("Running sensitivity analysis for culling rate optimization...\n")
set.seed(123) # for reproducibility
#-- set parameter space
np = np
low.threshold.range = rnorm(mean = Lambda.threshold.low, sd = 0.01, n = np)
high.threshold.range = rnorm(mean = Lambda.threshold.high, sd = 0.01, n = np)
# p0.range = round(rnorm(mean=p0, sd=25, n=np))
female.offtake.range = rnorm(mean=female.offtake, sd=5, n=np)

#-- generate two examples of random number from parmeter distributions
X1 = cbind.data.frame(low.threshold.range, high.threshold.range, p0.range, female.offtake.range)

#-- repeat sampling
low.threshold.range = rnorm(mean = Lambda.threshold.low, sd = 0.01, n = np)
high.threshold.range = rnorm(mean = Lambda.threshold.high, sd = 0.01, n = np)
# p0.range = round(rnorm(mean=p0, sd=25, n=np))
female.offtake.range = rnorm(mean=female.offtake, sd=5, n=np)

X2 = cbind.data.frame(low.threshold.range, high.threshold.range, p0.range, female.offtake.range)

#-- Create a dataframe containing the range of parameter values to test
sens_adjust_offtake = sensitivity::sobol2007(model = NULL, X1, X2, nboot = nboot.sens)

#-- run the model of all parameter values. The parameter values are passed to the adjust.offtake function, which applies the optimization algorithm to adjust offtake rates based on the specified thresholds and initial herd size. The resulting adjusted parameters are then used to project the herd dynamics using projectHerd2(). The outputs of these projections are summarized to obtain mean population size, standard deviation of population size, time to extinction, and initial population size (p0). These outputs are returned as a data frame for each set of parameter values.
sensitivity_wrapper <- function(l, sens_adjust_offtake) {
  
  #-- Apply over all parameter combinations in sens_adjust_offtake$X
  tmp_res <- apply(sens_adjust_offtake$X, 1, function(params) {
    
    lt <- params[1]  # low.threshold
    ht <- params[2]  # high.threshold
    # p0s <- params[3] # p0
    fo <- params[4]  # female.offtake
    
    #-- get adjusted parameters based on optimization function with current sampled parameter values
    adjusted_params <- lapply(l, 
                              FUN = adjust.offtake,
                              low.threshold = lt,
                              high.threshold = ht,
                              p0 = p0s,
                              female.offtake = fo,
                              sensitivity.test = TRUE)
    
    #-- Run projection
    tmp <- projectHerd2(listpar = adjusted_params, p0 = p0s)
    
    #-- summarize results
    summary.df <- pop_summary2(tmp, sex = FALSE, interval = "year")
    
    #-- Return summary data.frame
    data.frame(
      mean.pop.size = mean(summary.df$pop),
      pop.sd = sd(summary.df$pop), 
      time.to.extinct = summary.df$time[which(summary.df$pop == 0)[1]],
      p0 = p0s,
      low.threshold = lt,
      high.threshold = ht,
      female.offtake = fo
    )
  })
  
  #-- Combine results
  do.call(rbind.data.frame, tmp_res)
}
cat("Setting up parallel processing for sensitivity analysis...\n")
system.time({
  #-- Set up parallel processing
  cl <- parallel::makeCluster(n.cores, type = "PSOCK")
  # cl <- makeCluster(detectCores() - 1) # leave one core free
  #-- load necessary libraries on each cluster node
  parallel::clusterEvalQ(cl, {
    library(sensitivity)
    library(boot)
    library(HerdDynamics)
    library(dplyr)
    library(tidyr)
    library(stringr)
  })
  
  #-- export necessary functions and variables to cluster nodes
  parallel::clusterExport(cl, varlist = c("adjust.offtake", "get.off.adjust.poff", "projectHerd2", "pop_summary2", "sens_adjust_offtake"))
  
  # print message reporting how many cores are in use
  cat("Using", detectCores() - 1, "cores for sensitivity analysis...\n")
  cat("running sensitivity analysis across all parameter sets...\n")
  
  #-- run sensitivity analysis in parallel across all sets of parameters in listpar (except Baseline)
  
  # Sensitivity_results <- parLapply(cl, listpar[-c(16,32)], sensitivity_wrapper, sens_adjust_offtake)
  #-- Paynes models only
  Sensitivity_results <- parallel::parLapply(cl, listpar, sensitivity_wrapper, sens_adjust_offtake)
  #-- same as line above but includes a progress bar
  # Sensitivity_results <-  pbapply::pblapply(
  #   listpar[c(3:5, 19:21)],
  #   cl = cl,
  #   FUN = sensitivity_wrapper, 
  #   sens_adjust_offtake
  # )
  
  #-- tell sensitivity object about results for each parameter set. This will allow us to calculate sensitivity indices for each output variable (mean population size, population standard deviation, time to extinction) with respect to the input parameters (low.threshold, high.threshold, initial population and female offtake rate).
  sens_adjust_offtake_sobol.pop.size = lapply(Sensitivity_results, function(s){
    sensitivity::tell(sens_adjust_offtake, s$mean.pop.size)
  })
  
  sens_adjust_offtake_sobol.pop.sd = lapply(Sensitivity_results, function(s){
    sensitivity::tell(sens_adjust_offtake, s$pop.sd)
  })
  # 
  # sens_adjust_offtake_sobol.extinct = lapply(Sensitivity_results, function(s){
  #   sensitivity::tell(sens_adjust_offtake, s$time.to.extinct)
  # })
  
  sobol_res = list(
    pop.size = sens_adjust_offtake_sobol.pop.size,
    pop.sd = sens_adjust_offtake_sobol.pop.sd
    # time.to.extinct = sens_adjust_offtake_sobol.extinct
  )
  #-- send to data.frame
  sobol_indices = lapply(sobol_res, function(s) {
    do.call(rbind.data.frame, lapply(s, extract_indices)) %>%  
      mutate(
        taxon = str_split_i(rownames(.), "\\.", 1),
        strategy = str_split_i(rownames(.), "\\.", 2)
      ) 
  })
  
  sobol_indices <- lapply(sobol_indices, function(s) {
    rownames(s) <- NULL
    return(s)
  })
})

#-- clean up cluster
stopCluster(cl)

# saveRDS(sobol_indices, file = paste0(r_path, "data/sobol_indices.rds"))
saveRDS(sobol_indices, file = "output/sobol_indices.rds")
saveRDS(sobol_res, file = "output/sobol_res.rds")