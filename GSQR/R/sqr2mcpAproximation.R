sqr2mcpAproximation <-
function(gamm, taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, 
                                pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
  #################################################################################
  # call Fortran core
  gamma <- 2 * gamma/(k * nobs)
  gamma <- as.double(gamma)
  fit <- .Fortran("asqr2mcp", gamm, taux, k, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
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
