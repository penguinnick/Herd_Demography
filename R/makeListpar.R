#' Function to create listpar for input into projectHerd.R
#' This function creates listpar, which is a list of "param" tables, corresponding to the number of timesteps for which a projection will be made
#' 
#' Function does this: 
#' given nbcycle, calculate nbstep as nbcycle*nbphase to simulate fertility and mortality for time t=1 through time T=nbcycle
#' creates a list of T parms with varied mortality and fertility values using vary.fert.mort function
#' uses build.param function to create a param table for the i-th strategy for for time t through time T
#' 
#' @param param.props a list containing tcla (df), nbphase (integer), female.offtake (integer), truncated (boolean), parms (list: ages, parurition, part.age, prolificacy, MeanProlificacy, sdProlificacy, f.mortality, m.mortality
#' @param nbcycle integer from which nbstep is determined (the number of timesteps). When nbphase = 1, nbstep = nbcycle
#' @param offtake.mortality a list of mortality rates associated with different harvest strategies

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
          # truncated = truncated,
          correctionfec = TRUE,
          nbphase = nbphase 
        )$param
      })
    })
  })
}


