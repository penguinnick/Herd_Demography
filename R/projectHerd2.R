#' Function for projecting herd growth. Modified version of fproj (mmage) that builds listpar in the function allowing for function to run on multiple strategies in a list using lapply
#' @param all.param A list containing param data.frame created with fvh2par function (see mmage package) and initial.herd calculated with getLambda function
#' @param nbcycle an integer specifying how many years to run projection on 
#' @param nbphase an integer specifying how many phases to run projection on. Set to 12 for monthly, set to 1 for yearly. 
#' Default is NULL, which assumes the same ages were used for males and females when building tcla.
#' @param p0 integer the initial population size
#' @return a list containing: lambda of herd, a dataframe with sex proportions of the herd, and a dataframe with initial herd traits, including reproductive value by sex/age class
#' @references Lesnoff, M. (*), 2015. mmage: A R package for age-structured population matrix models. CIRAD, Montpellier, France. http://livtools.cirad.fr.



projectHerd2 = function (listpar, p0, nbcycle=NULL, nbphase = 1,  vec=TRUE, VariableEnvironment = FALSE ){
  nbcycle = length(listpar)
  nbstep = nbcycle * nbphase
  if ( length( listpar ) != nbstep ){
    stop( "listpar not adjusted to nbphase. length of listpar must be equal to nbcycle*nbphase" )
  }
  # nbstep = length(listpar)
  
  # listpar = vector("list", length=nbstep)
  # with(all.param,{
    
    #-- extract initial tcla from input params
    # tcla.ini = param$param[ , 1:5 ]
    
    #-- this section creates variable environment if True
    #-- assign each strategy's matrix to a list position
    # if ( VariableEnvironment ) {
    #   for( i in 1:nbstep ){
    #     
    #       output = with(param, { build.param( 
    #         tcla = tcla.ini, 
    #         parms = parms, 
    #         offtake = offtake, 
    #         nbphase = nbphase, 
    #         female.offtake = female.offtake,
    #         Inf.Mortality = "auto",
    #         prolificacyRate = "auto", 
    #         correctionfec = correctionfec, 
    #         phi = phi, truncated = FALSE )
    #        })
    #       listpar[[i]] = output$param
    #   }
    # } else {
    #   for(i in 1:nbstep){
    #     listpar[[i]] = param$param
    #   }
    # }
    
    # source("../R/mmagezcal.R")
    cal <- HerdDynamics::mmage.zcal(nbphase, nbstep + 1) 
    # cal <- zcal(nbphase, nbstep + 1) # run zcal function (mmage)
    cal
    z <- as.data.frame(listpar[[1]]) # set z as a listpar item
    Lf <- length(z$sex[z$sex == "F"]) - 1  # get length of male and female age classes
    Lm <- length(z$sex[z$sex == "M"]) - 1
    tcla <- z[, 1:5]   # subset z
    tcla2 <- tcla[tcla$class > 0, ]  # subset again
    matx <- matrix(0, Lf + Lm, nbstep + 1) # build matrix for storing results
    mat = HerdDynamics::mmage.fmat(z, Lf, Lm)
    # mat = fmat(z, Lf, Lm) # added 5/30/24
    matx[, 1] = mmage::feig(mat$A)$v*p0  # added 5/30/24 # populate with initial herd numbers
    # matx[, 1] <- initial.herd$xini  # populate with initial herd numbers
    matoff <- matdea <- matb <- X <- matrix(0, Lf + Lm + 2, nbstep)  # build offtake matrix
    matpar <- matrix(0, Lf + Lm, nbstep)  # build parturition matrix
    
    # source("../R/mmagefmat.R")  # improved fmat function
    # calculate step-wise herd demography changes #-- consider cleaning up with apply funs?
    for (j in 1:nbstep) {
      u <- HerdDynamics::mmage.fmat(listpar[[j]], Lf, Lm)  
      # u <- fmat(listpar[[j]], Lf, Lm)
      matx[, j + 1] <- u$A %*% matx[, j]
      X[, j] <- u$FEC %*% matx[, j]
      matb[c(1, Lf + 2), j] <- X[c(1, Lf + 2), j]
      matdea[, j] <- u$DEA %*% X[, j]
      matoff[, j] <- u$OFF %*% X[, j]
      matpar[, j] <- u$PAR[1, ] * matx[, j]
    }
    matstru <- mmage::fnorm(matx)
    matini <- matx[, 1:(ncol(matx) - 1)]
    matend <- matx[, 2:ncol(matx)]
    if (is.vector(matini)) {
      matini <- as.matrix(matini)
      matend <- as.matrix(matend)
    }
    fun <- function(tab1, tab2) {
      u <- as.data.frame(tab2)
      names(u) <- paste("t", 0:(ncol(u) - 1), sep = "")
      u <- cbind(tab1, u)
      row.names(u) <- 1:nrow(u)
      u
    }
    matx <- fun(tcla2, matx)
    matstru <- fun(tcla2, matstru)
    X <- fun(tcla, X)
    matini <- fun(tcla2, matini)
    matend <- fun(tcla2, matend)
    matdea <- fun(tcla, matdea)
    matoff <- fun(tcla, matoff)
    matpar <- fun(tcla2, matpar)
    matb <- fun(tcla, matb)
    vecx <- vecprod <- NULL
    if (vec) {
      z <- mmage::fproj2vec(matx, matstru, matini, matend, matdea, 
                     matoff, matpar, matb)
      u <- z$vecx
      u <- merge(cal$cal, u, by = "tim")
      u <- u[, -match("timend", names(u))]
      vecx <- u
      head(vecx)
      u <- z$vecprod
      u <- merge(cal$cal, u, by = "tim")
      u <- u[, -match("timend", names(u))]
      vecprod <- u
      head(vecprod)
    }
    res <- list(matx = matx, matstru = matstru, X = X, matini = matini, 
                matend = matend, matdea = matdea, matoff = matoff, matpar = matpar, 
                matb = matb, vecx = vecx, vecprod = vecprod)
    res
  # })
  }
