###########################################################################
#                                                                         #
#         Essential Functions for Analysis in JCAA Article                #
#         Run this file first to load all functions not already           #
#         loaded by other scripts.                                        #
###########################################################################

#-- source functions 
funs_to_load = c("culling_multiplot2", "mortality_histogram_gg", "mmagefm", "projectHerd2")
lapply(funs_to_load, function(f) {
  fun.file = paste0(r_path, "R/", f, ".R")
  source(fun.file)
})

# Rlib.dir <- "../R/x86_64-pc-linux-gnu-library/4.4"

#-- "using" function for installing and loading packages 
using <- function(pkgs,
                  repos = getOption("repos"),
                  quietly = TRUE, 
                  # lib = .libPaths()[2],
                  ...) {
  # pkgs: character vector of package names
  if (is.character(pkgs) && length(pkgs) == 1L) {
    pkgs <- c(pkgs)
  }
  
  # which are not installed?
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  
  if (length(to_install) > 0L) {
    install.packages(to_install, repos = repos, ...) #, lib = c("~/R/x86_64-pc-linux-gnu-library/4.4"), libs_only = TRUE)
  }
  
  # load all (now-installed) packages
  invisible(
    sapply(
      pkgs,
      require,
      character.only = TRUE,
      quietly = TRUE
    )
  )
}

#-- "ARR" function to calculate annual reproduction rate, included in HerdDynamics package
ARR = function(mean_litter_size, parturition_interval) {
  (mean_litter_size * 365) / parturition_interval
}

#-- "get_mortality" function to get vector of age structured mortality rates for male and female goats and sheep
get_mortality <- function(taxon, sex) {
  param.dat %>%
    filter(Taxon == taxon, Sex == sex) %>%
    pull(Mortality)
}

#-- "set_prolificacy" function to update param.dat data.frame with randomized prolificacy rates based on mean and sd from literature review
set_prolificacy <- function(taxon, i) {
  idx <- with(param.dat, ParturitionRate > 0 & Taxon == taxon)
  
  param.dat$Prolificacy[idx] <<- generate_prolificacy_rates(
    meanPro = NetPro$mn.pro[i],
    sdPro   = NetPro$sd.pro[i],
    n       = 6
  )
}

#-- "generate_prolificacy_rates" function to generate vector of prolificacy rates for adult age classes based on mean and sd from literature review--
generate_AdultMortalityRates <- function(Mort, n = 6) {
  adult.mort = Mort[4:length(Mort)]
  sort(generate_prolificacy_rates(
    meanPro = mean(adult.mort),
    sdPro = sd(adult.mort),
    n = 6
  ),
  decreasing = FALSE)
}

#-- "vary.fert.mort" function to vary prolificacy and mortality
vary.fert.mort = function(parms, n = 6) {
  parms$prolificacy = generate_prolificacy_rates(
    meanPro = parms$MeanProlificacy,
    sdPro = parms$sdProlificacy,
    n = n
  )
  
  #-- vary parturition rate
  # parms$parturition = ARR(mean(parms$prolificacy), parms$parturition.Interval)
  parms$parturition = ARR(parms$prolificacy, parms$parturition.Interval)
  
  # generate_AdultMortalityRates <- function(Mort = f.mortality, n = 6) {
  #   adult.mort = Mort[4:length(Mort)]
  #   sort(generate_prolificacy_rates(
  #     meanPro = mean(adult.mort),
  #     sdPro = sd(adult.mort),
  #     n = 6
  #   ),
  #   decreasing = FALSE)
  # }
  
  #-- update male and female mortality rates
  parms$f.mortality = c(
    HerdDynamics::generate_infant_mortality_rates(parms$f.mortality),
    generate_AdultMortalityRates(Mort = parms$f.mortality)
  )
  
  parms$m.mortality = c(
    HerdDynamics::generate_infant_mortality_rates(parms$m.mortality),
    generate_AdultMortalityRates(Mort = parms$m.mortality)
  )
  return(parms)
}

#-- "make.listpar" functio to create list of parameters for each offtake strategy with varying fertility and mortality for sensitivity analysis
make.listpar = function(param.props, nbcycle, offtake.mortality) {
  with(param.props, {
    #-- get number of timesteps
    nbstep = nbcycle * nbphase
    #-- create list of parms
    parms.rep = replicate(n = nbstep,
                          vary.fert.mort(parms = parms),
                          simplify = FALSE)
    #-- create listpar
    lapply(offtake.mortality, function(o) {
      lapply(parms.rep, function(p) {
        build_param(
          tcla = tcla,
          parms = p,
          offtake = o,
          female.offtake = female.offtake,
          correctionfec = TRUE,
          nbphase = nbphase
        )$param
      })
    })
  })
}

#-- "wrapper.repro": a helper function to get demography information from every listpar table
wrapper.repro = function(listpar, out = c("lambda", "sex")) {
  out <- match.arg(out)
  lapply(seq_along(listpar), function(l) {
    p = listpar[[l]]
    lapply(p, function(s) {
      sapply(s, function(x) {
        r = get_lambda(x, tcla = tcla, p0 = p0) #$lambda
        if (out == "lambda") {
          return(r$lambda)
        } else {
          if (out == "sex") {
            return(r$sex.proportion[1, 2])
          }
        }
      })
    })
  })
}

#-- adapted wrapper.repro function for parallel computing
wrapper.repro_parallel = function(listpar, out = c("lambda", "sex"), tcla, p0) {
  out <- match.arg(out)
  p <- listpar[[l_index]]
  
  lapply(p, function(s) {
    sapply(s, function(x){
      r <- get_lambda(x, tcla = tcla, p0 = p0)
      if (out == "lambda") {
        return(r$lambda)
      } else {
        return(r$sex.proportion[1, 2])
      }
    })
  })
}

#-- "mean_function" calculates the mean of a given vector of data, used for bootstrapping
mean_function <- function(data, indices) {
  # This function will be applied to resampled data
  return(mean(data[indices]))
}

#-- "boot.fun" function to compute bootstrapped mean and 95% confidence intervals
boot.fun = function(s, n) {
  # Perform bootstrapping
  bootstrap_results = boot::boot(data = s,
                           statistic = mean_function,
                           R = n)
  # Obtain bootstrapped confidence interval
  boot_conf_interval = boot::boot.ci(bootstrap_results, type = "perc")
  #-- returns mean, lower and upper ci
  return(
    data.frame(
      t0 = boot_conf_interval$t0,
      low = boot_conf_interval$percent[4],
      up = boot_conf_interval$percent[5]
    )
  )
}


#-- "facet_factor_fun" function to create factors for plotting offtake strategies by author and strategy in culling_multiplot2 function
facet_factor_fun <- function(dat){
  dat %>% 
    mutate(Author = case_when(
      strategy %in% c("Security", "Energy") ~ "Redding (1981)",
      strategy %in% c("Meat", "Milk", "Wool") ~ "Payne (1973)",
      strategy %in% c("MeatA", "MeatB", "MilkA", "MilkB", "Fleece") ~ "Vigne and Helmer (2007)",
      strategy == "Baseline" ~ "No Offtake",
      .default = "Neolithic" )) %>%
    mutate(Author = factor(
      Author, 
      levels = c( 
        "Vigne and Helmer (2007)", 
        "Payne (1973)", 
        "Redding (1981)",
        "Neolithic", 
        "No Offtake"
        )
      ),
      strategy = factor(
        strategy, 
        levels = c(
          names(offtake.mortality)[6:10], 
          names(offtake.mortality)[3:5],
          names(offtake.mortality)[1:2],
          "Smilčić EN", "Smilčić MN", "Benkovac-Barice MN",
          "Islam Grčki MN", "Zemunik Donji MN", "Baseline")
        )
    )
}

#-- "xini.to.data.frame" wrapper function to get a dataframe with herd structure from repro object
xini.to.data.frame = function(repro, ageClasses, ages) {
  xini.list = lapply(repro, function(r) {
    # given repro, a list of reproductive params for n species,
    lapply(r, function(x) {
      with(x, {
        # summarize initial herds (total males and females for each age class)
        xini = initial.herd %>%
          group_by(class) %>%
          reframe(n = as.integer(sum(xini)))
      })
    })
  })
  xini.list = unlist(xini.list, recursive = F) # unlist species groupings
  taxon.strat = str_split(names(xini.list), "\\.") # create list of species and strategies
  # create Taxon, strategy, and AgeClass fields for each df in xini.list
  for (i in 1:length(xini.list)) {
    ts = taxon.strat[[i]]
    taxon = rep(unlist(ts)[1], 9)
    strat = rep(unlist(ts)[2], 9)
    xini.list[[i]]$Taxon = taxon
    xini.list[[i]]$strategy = strat
    xini.list[[i]]$AgeClass = ageClasses[-1]
    xini.list[[i]]$Age = as.character(ages)
  }
  # call all dataframes into a single df
  xini.df = do.call(rbind.data.frame, xini.list)
  # make strategy field a factor
  xini.df$strategy = factor(xini.df$strategy, levels = c(names(offtake.mortality)))
  xini.df$AgeClass = factor(xini.df$AgeClass, levels = ageClasses[-1])
  return(xini.df)
}


#-- "pop_summary2" function to summarize dynamics of the total population through time
pop_summary2 = function(x, sex = FALSE, interval = "year") {
  vecx = x$vecx
  if (interval == "month") {
    if (sex) {
      out = aggregate(x ~ tim + sex , data = vecx, FUN = sum)
      colnames(out) = c("time", "sex", "pop")
    } else {
      out = aggregate(x ~ tim, data = vecx, FUN = sum)
      colnames(out) = c("time", "pop")
    }
  }
  if (interval == "year") {
    if (sex) {
      out = aggregate(x ~ cycle + sex , data = vecx, FUN = sum)
      colnames(out) = c("time", "sex", "pop")
    } else {
      out = aggregate(x ~ cycle , data = vecx, FUN = sum)
      colnames(out) = c("time", "pop")
    }
  }
  return(out)
}

#-- "res.to.df" function to gather results into a single data.frame
res.to.df = function(tot.pop.res) {
  df = do.call(rbind.data.frame, args = c(tot.pop.res, make.row.names = FALSE))
  strats = names(tot.pop.res)
  df$strategy = unlist(lapply(strats, FUN = rep, (nrow(df) / length(strats))))
  ts = str_split(df$strategy, "\\.", simplify = T)
  df$Taxon = ts[, 1]
  df$strategy = ts[, 2]
  return(df)
}

#-- "stochastic.rep" function to replicate projections of herd dynamics with varying fertility and mortality parameters for sensitivity analysis
stochastic.rep = function(p0,
                          param.props,
                          offtake.mortality,
                          nbcycle) {
  #-- create new set of environmental parameters
  lp = make.listpar(
    param.props = param.props,
    nbcycle = nbcycle,
    offtake.mortality = list(offtake.mortality)
  )
  #-- reproject herd growth
  result = projectHerd2(listpar = unlist(lp, recursive = F) , p0 = p0)
  #-- summarize results
  pop_summary2(result, sex = FALSE, interval = "year")$pop
}

run_stochastic_rep <- function(o, p0, param.props, offtake.mortality, nbcycle) {
  stochastic.rep(
    p0 = p0,
    param.props = param.props,
    offtake.mortality = o,
    nbcycle = nbcycle
  )
}

#-- "list.to.mat" helper function to transform lists into matrix
list.to.mat = function(l) {
  n.col = length(l)
  l.mat = matrix(unlist(l), ncol = n.col)
  rownames(l.mat) = unique(names(unlist(l)))
  l.mat
}

#-- "matrix.to.df" helper function to create df from matrix output of replicated stochastic.rep function
matrix.to.df = function(mat) {
  df = as.data.frame(mat)
  df$time = as.numeric(str_replace_all(rownames(df), "[:alpha:]|[:punct:]", ""))
  df$strategy = str_replace_all(rownames(df), "[:digit:]", "")
  df %>% pivot_longer(cols = c(-time, -strategy))
}

#-- "summarize.pop.df" summarize dataframe helper, computes mean and confidence intervals for each strategy at each time step. Used for plotting
summarize.pop.df = function(df) {
  df %>%
    group_by(taxon, strategy, time) %>%
    summarise(
      mean = mean(value),
      low = quantile(value, 0.025),
      up = quantile(value, 0.975)
    )
}

### Optimization functions ###
#-- function to adjust offtake based on poff in param table
get.off.adjust.poff = function(phi, param.ref, m = 1) {
  u = param.ref
  zf = u$poff[u$sex == "F"]
  zm = u$poff[u$sex == "M"]
  u$poff = c(phi * zf, zm)
  Lf = length(u$sex[u$sex == "F"]) - 1
  Lm = length(u$sex[u$sex == "M"]) - 1
  param = u
  A = mmage::fmat(param, Lf, Lm)$A
  (mmage::feig(A)$lambda - m)^2
}

#-- function for adjusting p. female offtake
adjust.offtake = function(in.param,
                          low.threshold,
                          high.threshold,
                          female.offtake = female.offtake,
                          sensitivity.test = FALSE, # if TRUE, tests a range of female offtake rates for sensitivity test
                          p0) {
  #-- calculate phi.opt
  phi.opt = optimize(f = get.off.adjust.poff,
                     param.ref = in.param,
                     interval = c(0, 5))$minimum
  
  if(sensitivity.test) {
    #-- get unadjusted male offtake
    m.off = in.param$poff[in.param$sex == "M" & in.param$class > 0]
    
    #-- apply female offtake rate
    f.off = m.off * (female.offtake/100)
  } else {
    #-- get unadjusted female offtake
    f.off = in.param$poff[in.param$sex == "F" & in.param$class > 0]
  }
  
  
  #-- extract tcla
  tcla1 = in.param[, c(1:5)]
  
  #-- extract lclass
  # lclass1 = in.param$lclass[ in.param$sex=="F" & in.param$class > 0 ]
  
  #-- get lambda
  lambda = HerdDynamics::get_lambda(param = in.param, tcla = tcla1, p0 = p0)$lambda
  
  #-- compare lambda to threshold value
  if (lambda < low.threshold | lambda >= high.threshold) {
    #-- if lambda below threshold, adjust female offtake
    in.param$poff[in.param$sex == "F" &
                    in.param$class > 0] = f.off * phi.opt
  } else {
    in.param$poff[in.param$sex == "F" & in.param$class > 0] = f.off
  }
  
  if(sensitivity.test) {
    in.param$p0 = p0
  }
  return(in.param)
}

#-- Rewritten fprod code from mmage package
fprod = function (formula, vecprod, wnum = NULL, wden = NULL, digits = c(1,3)) {
  if (length(digits) == 1) 
    digits <- rep(digits, 2)
  f <- formula(deparse(formula))
  nam.var <- all.vars(f)
  nam.var
  z <- vecprod
  # z <- results$Energy$vecprod
  z[is.na(z)] <- 0
  if (is.null(wnum)) 
    wnum <- rep(1, nrow(z))
  else wnum <- z[, wnum]
  if (is.null(wden)) 
    wden <- rep(1, nrow(z))
  else wden <- z[, wden]
  nbcycle <- max(z$cycle)
  nbphase <- max(z$phase)
  nbcycle
  nbphase
  if ("cycle" %in% nam.var) {
    nbc <- 1
    z$xini <- ifelse(z$phase == 1, z$xini, 0)
    z$xend <- ifelse(z$phase == nbphase, z$xend, 0)
  }  else {
    nbc <- nbcycle
    z$xini <- ifelse(z$cycle == 1 & z$phase == 1, z$xini, 
                     0)
    z$xend <- ifelse(z$cycle == nbcycle & z$phase == nbphase, 
                     z$xend, 0)
  }
  z$delta <- wnum * z$delta
  z$dea <- wnum * z$dea
  z$off <- wnum * z$off
  z$xmean <- wden * z$xmean
  z$prod <- z$delta + z$off
  vecprod2 <- z
  head(vecprod2)
  newf <- formula(paste("cbind(xini, xend, xmean, delta, dea, off, prod) ~", 
                        f[2]))
  newf
  z <- vecprod2
  z <- aggregate(newf, data = z, FUN = sum)
  if (length(nam.var) > 1) 
    z <- z[do.call(order, z[, nam.var]), ]
  z$xmean <- z$xmean/(nbc * nbphase)
  u <- z$xmean * nbc
  z$rdelta <- z$delta/u
  z$rdea <- z$dea/u
  z$roff <- z$off/u
  z$rprod <- z$prod/u
  z$m <- (z$xend/z$xini)^(1/nbc)
  z$nbcycl <- rep(nbc, nrow(z))
  u <- c("xini", "xend", "xmean", "delta", "dea", "off", "prod")
  z[, match(u, names(z))] <- round(z[, match(u, names(z))], 
                                   digits = digits[1])
  u <- c("rdelta", "rdea", "roff", "rprod", "m")
  z[, match(u, names(z))] <- round(z[, match(u, names(z))], 
                                   digits = digits[2])
  z
}

#-- function for formatting Tables 4 and 5 
lambda.tables = function(tab) {
  tf = as_grouped_data(tab, groups = "Taxon")
  as_flextable(tf) %>%
    separate_header() %>%
    set_header_labels(
      name = "",
      MeatA = "Meat A",
      MeatB = "Meat B",
      MilkA = "Milk A",
      MilkB = "Milk B"
    ) %>%
    bold(
      j = 1,
      i = ~ !is.na(Taxon),
      bold = TRUE,
      part = "body"
    ) %>% # make Taxon bold
    bold(part = "header", bold = TRUE) %>%  # make headers bold
    colformat_double(i = c(2, 6), digits = 3) %>%
    colformat_double(i = c(3, 4, 7, 8), digits = 2) %>%
    autofit() %>%
    set_table_properties(layout = "autofit", align = "center") %>%
    paginate(
      group = "Taxon",
      # keep Sample.ID group across pages
      group_def = "rle"
    ) %>%
    theme_vanilla() %>%
    align(j = 1, align = "left")
}

#-- Function to calculate moving average
calculate_moving_average <- function(time_series,
                                     window_size,
                                     type = c("simple", "exponential"),
                                     alpha = 0.1) {
  type <- match.arg(type)
  
  if (type == "simple") {
    # Calculate simple moving average
    moving_average <- zoo::rollmean(time_series,
                                    k = window_size,
                                    fill = NA,
                                    align = "right")
  } else if (type == "exponential") {
    # Calculate exponential moving average
    moving_average <- zoo::rollapply(
      time_series,
      width = window_size,
      FUN = function(x) {
        n <- length(x)
        weights <- (1 - alpha)^(n:1)
        sum(weights * x) / sum(weights)
      },
      fill = NA,
      align = "right"
    )
  }
  
  return(moving_average)
}

#-- Function to detect peaks
detect_peaks <- function(x, span = 3) {
  z <- embed(x, span)
  s <- span %/% 2
  v <- max.col(z, ties.method = "first") == 1 + s
  c(rep(FALSE, s), v, rep(FALSE, s))
}

#-- Function to detect valleys
detect_valleys <- function(x, span = 3) {
  z <- embed(x, span)
  s <- span %/% 2
  v <- apply(z, 1, which.min) == 1 + s
  c(rep(FALSE, s), v, rep(FALSE, s))
}


plot.m = function(taxon = c("goat", "sheep")) {
  tax = match.arg(taxon)
  all.m.df %>%
    facet_factor_fun() %>% 
    filter(Taxon == tax) %>%
    group_by(Taxon, strategy, offtake) %>%
    reframe(
      cycle,
      Taxon,
      strategy,
      offtake,
      Author,
      m,
      m10 = calculate_moving_average(m, 10, type = "simple"),
      peak = detect_peaks(m10, span = 31),
      valley = detect_valleys(m10, span = 31)
    ) %>%
    filter(!is.na(m10)) %>%
    ggplot(aes(
      x = cycle,
      y = m10,
      group = offtake,
      colour = offtake
    )) +
    geom_line(lwd = 0.2) +
    geom_hline(
      yintercept = Lambda.threshold.low,
      linetype = "solid",
      lwd = 0.2,
      alpha = 0.8
    ) +
    geom_hline(
      yintercept = Lambda.threshold.high,
      linetype = "dashed",
      lwd = 0.2,
      alpha = 0.8
    ) +
    geom_point(
      data = . %>% filter(peak),
      aes(x = cycle, y = m10),
      shape = 24,
      fill = rcartocolor::carto_pal(name = "Safe")[12],
      alpha = 0.8,
      size = 0.8,
      show.legend = F
    ) +
    geom_point(
      data = . %>% filter(valley),
      aes(x = cycle, y = m10),
      shape = 25,
      fill = rcartocolor::carto_pal(name = "Safe")[12],
      alpha = 0.8,
      size = 0.8,
      show.legend = F
    ) +
    scale_color_manual(values = cbPalette[c(4, 6)], labels = c("m (optimized)", "λ (unadjusted)")) +
    ylab(paste0(expression("λ and m, ("), tax, ")")) + # "Lambda (sheep)") +
    xlab("Time (years)") +
    xlim(c(10, 200)) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    facet_manual(
      facets = vars(Author, strategy), 
      design = design, 
      strip = strip_nested(),
      scales = "free_y") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 8),
      axis.text.x = element_text(hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      line = element_line(linewidth = 0.2),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.3), 
      legend.title = element_blank(),
      strip.text = element_text(size = 6),
      strip.background = element_rect(fill = "white", color = "grey60"),
      panel.spacing = unit(0.2, "lines")
    ) 
}


#-- Function for plotting projections of herd dynamics under different offtake strategies, with facets for strategy
plot.projection = function(df, offtake.list) {
  df %>% 
    facet_factor_fun() %>%
    mutate(
      Taxon = factor(Taxon, levels = c("goat", "sheep"))) %>%
    ggplot() +
    geom_col(aes(time, pop, fill = Taxon), alpha = 0.8) +
    scale_fill_manual(
      breaks = c("goat", "sheep"), 
      values = taxon_pal, # c("#d95f02", "#1b9e77"),
      labels = c("Goat", "Sheep")
    ) +
    labs(x = "Time (years)", y = "Number of animals") +
    facet_manual(
      facets = vars(Author, strategy), 
      design = design, 
      strip = strip_nested(),
      scales = "free_y") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 8),
      axis.text.x = element_text(hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      line = element_line(linewidth = 0.2),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.3), 
      legend.title = element_blank(),
      strip.text = element_text(size = 6),
      strip.background = element_rect(fill = "white", color = "grey60"),
      panel.spacing = unit(0.2, "lines")
    ) 
}

#-- Extracts Sobol indices from sobol2007 object
extract_indices = function(sobol_obj){
  interactions = sobol_obj$T[,1] - sobol_obj$S[,1]
  data.frame(
    param = rownames(sobol_obj$S),
    main_effect = sobol_obj$S[,2], # extract bias-corrected estimator
    main_effect_low = sobol_obj$S[,4],
    main_effect_high = sobol_obj$S[,5],
    total_effect = sobol_obj$T[,1], # extract original estimator
    total_effect_low = sobol_obj$T[,4],
    total_effect_high = sobol_obj$T[,5],
    interactions = interactions
  )
}



