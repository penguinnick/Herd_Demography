# rewritten fm function from mmage

mmage.fm = function (formula, vecprod, w = NULL, digits = c(1, 3)){
  
  # digits=c(1,3)
  # formula = ~1
  # vecprod = results$Wool$vecprod
  # w=NULL
  
  if (length(digits) == 1){
    digits <- rep(digits, 2)
  }
    
  f <- formula(deparse(formula))
  nam.var <- all.vars(f)
  nam.var
  z <- vecprod
  z[is.na(z)] <- 0
  if (is.null(w)){
    w <- rep(1, nrow(z))
  } else {
    w <- z[, w]
  }
  nbcycle <- max(z$cycle)
  nbphase <- max(z$phase)
  nbcycle
  nbphase
  z$xini <- w * z$xini
  z$xend <- w * z$xend
  z$xmean <- w * z$xmean
  if ("cycle" %in% nam.var) {
    nbc <- 1
    z$xini <- ifelse(z$phase == 1, z$xini, 0)
    z$xend <- ifelse(z$phase == nbphase, z$xend, 0)
  } else {
    nbc <- nbcycle
    z$xini <- ifelse(z$cycle == 1 & z$phase == 1, z$xini, 0)
    z$xend <- ifelse(z$cycle == nbcycle & z$phase == nbphase, z$xend, 0)
  }
  vecprod2 <- z
  newf <- formula(paste("cbind(xini, xend, xmean) ~", f[2]))
  newf
  z <- vecprod2
  # z <- aggregate(formula = newf, data = z, FUN = sum) # this fails!
  z <- aggregate(newf, data = z, FUN = sum)
  z$xmean <- z$xmean/(nbc * nbphase)
  if ("cycle" %in% nam.var) {
    z <- z[order(z$cycle), ]
    u <- by(data = z, INDICES = z[, "cycle"], function(x) mmage::fnorm(x$xini))
    z$stru.ini <- unlist(u)
    u <- by(data = z, INDICES = z[, "cycle"], function(x) mmage::fnorm(x$xend))
    z$stru.end <- unlist(u)
    u <- by(data = z, INDICES = z[, "cycle"], function(x) mmage::fnorm(x$xmean))
    z$stru.mean <- unlist(u)
  }  else {
    z$stru.ini <- mmage::fnorm(z$xini)
    z$stru.end <- mmage::fnorm(z$xend)
    z$stru.mean <- mmage::fnorm(z$xmean)
  }
  z$stru.ini <- 100 * z$stru.ini
  z$stru.end <- 100 * z$stru.end
  z$stru.mean <- 100 * z$stru.mean
  z$m <- (z$xend/z$xini)^(1/nbc)
  z$nbcycl <- rep(nbc, nrow(z))
  if (length(nam.var) > 1) {
    z <- z[do.call(order, z[, nam.var]), ]
  }
  z
  u <- c("xini", "xend", "xmean", "stru.ini", "stru.end", "stru.mean")
  z[, match(u, names(z))] <- round(z[, match(u, names(z))], 
                                   digits = digits[1])
  u <- c("m")
  z[, match(u, names(z))] <- round(z[, match(u, names(z))], 
                                   digits = digits[2])
  row.names(z) <- 1:nrow(z)
  z
}