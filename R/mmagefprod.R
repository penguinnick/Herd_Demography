#' Rewritten fprod code from mmage

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
