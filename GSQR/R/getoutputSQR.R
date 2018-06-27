getoutputSQR <-
function(fit, maxit, pmax, nvars, vnames,taux) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = print(errmsg$msg, 
                                                                              call. = FALSE))
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(fit$beta[seq(nvars * nalam)], nvars, nalam, dimnames = list(vnames, 
                                                                               stepnames))
    df <- apply(abs(beta) > 0, 2, sum)
  } else {
    beta <- matrix(0, nvars, nalam, dimnames = list(vnames, stepnames))
    df <- rep(0, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}
