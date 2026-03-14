###########################################################################
#                                                                         #
#         PARAMETERS for JCAA Article                                     #
#                                                                         #
# How to use:                                                             #
#         Execute code to reproduce results and figures in the article.   #
#         Parameters can be changed here before running 02_ANALYSIS.R     #
#                                                                         #                
###########################################################################
# n.cores = 128 # number of cores to use for parallel processing. Set to 1 for no parallelization.
n.cores = parallel::detectCores() - 1 # use all cores except one for parallel processing
#####################--MORTALITY, FERTILITY, AND CULLING RATES--######################
#-- read in lambing and kidding rates
# pro.dat = read.csv("./data/prolificacy.csv")
pro.dat = read.csv(paste0(r_path,"data/prolificacy.csv"))
#-- summarize prolificacy rates from the literature
NetPro = pro.dat |> 
  dplyr::group_by(Taxon) |> 
  dplyr::reframe(mn.pro = mean(LitterSize),
          sd.pro = sd(LitterSize))

#-- read in parameters file for mortality, prolificacy, and parturition rates
# param.dat = read.csv("./data/parameters2.csv")
param.dat = read.csv(paste0(r_path,"data/parameters2.csv"))


#-- Ages, parturition interval, and age of first parturition set here
ages = HerdDynamics::Payne_ages$ages  
parturition.Interval = 300 # 300 days between births
part.age = 2

#-- This number is used to account for culling of females due to infertility. If Null, no female offtake assumed.
female.offtake = 15 # percentage of culling rate applied to females

#####################--PROJECTION RULES--######################
#-- set number of cycles and initial population size
nbcycle = 200 # YEARS
p0 = 150       # initial population size

#####################--BOOTSTRAP AND REPLICATION COUNTS--######################
#-- set number of bootstrap replicates
n.boot = 1000

#-- number of replications for confidence intervals
nbrep = 100 

#####################--SENSITIVITY ANALYSIS--######################
np = 10000 # number of parameter sets to run for sensitivity analysis
nboot.sens = 500
#####################--MISC--######################
#-- color palette for figures
cbPalette = rcartocolor::carto_pal(10, "Safe")
# #-- Facet structure for plots
design <- matrix(c(
  11, 12, 13, 14, 15,
  1, 2, 3, 4, 5,
  6, 7, 8, NA, NA,
  9, 10, 16, NA, NA
), nrow = 4, byrow = TRUE )

#-- labels for Fig 11
param_labels <- c(
  "female.offtake.range"  = "Female offtake rate",
  "high.threshold.range"  = expression(High~lambda~threshold),
  "low.threshold.range"   = expression(Low~lambda~threshold),
  "p0.range"              = "Initial population (p0)"
)

taxon_pal = c("#d95f02", "#1b9e77")
