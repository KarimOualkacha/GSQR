sqr1mcp <-
function(gamm, taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, 
                 pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
  #################################################################################
  # call Fortran core
  gamma <- 2 * gamma/((c/max(taux,1-taux)) * nobs)
  gamma <- as.double(gamma)
  fit <- .Fortran("sqr1mcp", gamm, taux, c, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                  nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output

  
  outlist <- getoutputSQR(fit, maxit, pmax, nvars, vnames,taux)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
  class(outlist) <- c("hsvm")
  outlist
}
