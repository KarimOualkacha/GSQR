Gsqr <-
function(x, y, group = NULL,method=c("GLasso","GMcp","GScad","AGMcp","AGScad"),check=c("f1","f2"), nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.05, 0.001), 
                 lambda = NULL, pf = sqrt(bs), dfmax = as.integer(max(group)) + 1, 
                 pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 0.0001, maxit = 3e+08, 
                 gamm = ifelse(((method == "GMcp")|(method == "AGMcp")), 3, 4) ,taux=0.5,delta=2,intercept=TRUE) {
  #################################################################################
  this.call <- match.call()
  check <- match.arg(check)
  c=k=delta
  #\tDesign matrix setup, error checking######################
  this.call <- match.call()
  if (!is.matrix(x)) 
    stop("x has to be a matrix")
  
  if (any(is.na(x))) 
    stop("Missing values in x not allowed!")
  y <- drop(y)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames)) 
    vnames <- paste("V", seq(nvars), sep = "")
  
  if (length(y) != nobs) 
    stop("x and y have different number of rows")
  
  if (!is.numeric(y)) 
    stop("The response y must be numeric. Factors must be converted to numeric")
  
  #################################################################################
  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else if (length(group) != nvars) 
    stop("group length does not match the number of predictors in x")
  
  bn <- as.integer(max(group))
  bs <- as.integer(as.numeric(table(group)))
  
  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) 
    stop("Groups must be consecutively numbered 1,2,3,...")
  
  ix <- rep(NA, bn)
  iy <- rep(NA, bn)
  j <- 1
  for (g in 1:bn) {
    ix[g] <- j
    iy[g] <- j + bs[g] - 1
    j <- j + bs[g]
  }
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)
  #################################################################################
  #  get upper bound
  gamma <- rep(NA, bn)
  for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
  ##for (g in 1:bn) gamma[g] <- 2*max(taux,1-taux)/c
  #################################################################################
  #parameter setup
  if (taux < 0) 
    stop("taux must be non-negtive")
  taux <- as.double(taux)
  if (c < 0) 
    stop("c must be non-negtive")
  c <- as.double(c)
  
  if (length(pf) != bn) 
    stop("The size of group-lasso penalty factor must be same as the number of groups")
  maxit <- as.integer(maxit)
  pf <- as.double(pf)
  gamm <- as.double(gamm)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  r <- y
  #################################################################################
  #lambda setup
  if (method!="SGLasso"){
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1) 
        stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)
    } else {
      #flmin=1 if user define lambda
      flmin <- as.double(1)
      if (any(lambda < 0)) 
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))
      nlam <- as.integer(length(lambda))
    }
  }else{
    if (is.null(lambda)) {
      ulam <- as.double(LambdaMax(x, y, index=group, alpha = alpha, min.frac = 0.05, nlam = 20, c=c, taux=taux,check=check))
      nlam <- as.integer(length(ulam))
    } else {
      if (any(lambda < 0)) 
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))
      nlam <- as.integer(length(lambda))
    }    
  }
  intr <- as.integer(intercept)
  #################################################################################
  # call Fortran core
  
  if (method == "GLasso") {
    fit <- switch(check, f1 = sqr1lasso(taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                        nlam, flmin, ulam, eps, maxit, vnames, group, intr), 
                  f2 = sqr2lasso(taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                 nlam, flmin, ulam, eps, maxit, vnames, group, intr))
  }
  if (method == "SGLasso") {
    fit <- switch(check, f1 = sqr1Slasso(alpha,taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                         nlam, ulam, eps, maxit, vnames, group, intr), 
                  f2 = sqr2Slasso(alpha,taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                  nlam, ulam, eps, maxit, vnames, group, intr))
  }
  if (method == "GMcp") {
    fit <- switch(check, f1 = sqr1mcp(gamm,taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                      nlam, flmin, ulam, eps, maxit, vnames, group, intr),
                  f2 = sqr2mcp(gamm,taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                               nlam, flmin, ulam, eps, maxit, vnames, group, intr))   
  }
  if (method == "AGMcp") {
    fit <- switch(check, f1 = sqr1mcpAproximation(gamm,taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                                  nlam, flmin, ulam, eps, maxit, vnames, group, intr),
                  f2 = sqr2mcpAproximation(gamm,taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                           nlam, flmin, ulam, eps, maxit, vnames, group, intr))   
  }
  if (method == "GScad") {
    fit <- switch(check, f1 = sqr1scad(gamm, taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                       nlam, flmin, ulam, eps, maxit, vnames, group, intr), 
                  f2 = sqr2scad(gamm, taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                nlam, flmin, ulam, eps, maxit, vnames, group, intr))
  }
  if (method == "AGScad") {
    fit <- switch(check, f1 = sqr1scadAproximation(gamm,taux, c, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                                   nlam, flmin, ulam, eps, maxit, vnames, group, intr),
                  f2 = sqr2scadAproximation(gamm,taux, k, bn, bs, ix, iy, gamma, nobs, nvars, x, y, pf, dfmax, pmax, 
                                            nlam, flmin, ulam, eps, maxit, vnames, group, intr)) 
  } 
  #################################################################################
  # output
  if (is.null(lambda)) 
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  class(fit) <- c("Gsqr", class(fit))
  fit
}
