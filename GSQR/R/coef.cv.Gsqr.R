coef.cv.Gsqr <-function(cv, s = c("lambda.1se", "lambda.min")) {

    cvm=cv$cv
    cvsd=cv$cv.error
    cvmin <- min(cvm)
    idmin <- cvm <= cvmin
    lambda.min <- max(cv$lambda[idmin])
    idmin <- match(lambda.min, cv$lambda)
    semin <- (cvm + cvsd)[idmin]
    idmin <- cvm <= semin
    lambda.1se <- max(cv$lambda[idmin])

  if (s=="lambda.1se"){
    l = which(cv$lambda == lambda.1se)
  }else{
    l = which(cv$lambda == lambda.min)
  }
  list(beta0=cv$finalfit$b0[l],beta=as.vector(t(cv$finalfit$beta)[l,]))
}



