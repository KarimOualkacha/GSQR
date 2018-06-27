cv.Gsqr <-
function (x, y, group, Kfold = 5,delta=2, taux=0.5, check = c("f1","f2"),gamm = ifelse(method == "GMcp", 3, 4),
  method = c("GLasso","SGLasso","GMcp","GScad","AGMcp","AGScad"), plot.it=FALSE, ...)
  {
    
    check <- match.arg(check)
    method <- match.arg(method)
    n=dim(x)[1]
    p=dim(x)[2]
    bs <- as.integer(as.numeric(table(group)))
    all.folds <- split(sample(1:length(y)), rep(1:Kfold, length = length(y)))
    object<-Gsqr(x=x, y=y, group=group,check=check,delta=delta, taux=taux, method=method ,gamm=gamm,...)
    finalfit <- object
    lambda <- object$lambda
    residmat <- matrix(0, length(lambda),Kfold)
    for (i in seq(Kfold)) {
      omit <- all.folds[[i]]
      object <- Gsqr(x=x[-omit, ],y=y[-omit],lambda=lambda, group=group,check=check,delta=delta, taux=taux, method=method,gamm=gamm,...)
      nbeta <- rbind(object$b0, object$beta)
      for (j in 1:length(lambda)){
        ypred <- cbind(1, x[omit, ]) %*% nbeta[,j]
        ypred <- y[omit]-ypred
        ypred[ypred<0] <- (1-taux)*abs(ypred[ypred<0])
        ypred[0<ypred] <- (taux)*abs(ypred[0<ypred])
        residmat[j, i] <- mean(ypred)
      }
    }
    cv <- apply(residmat, 1, median)
    cv.error <- sqrt(apply(residmat, 1, var)/(Kfold-1))
    
    if (plot.it) {
    matplot(lambda, cv, type = "l",main = paste("CV.Q-",method) ,ylim = range(cv, cv + cv.error, cv - cv.error))
    error.bars(lambda, cv + cv.error, cv - cv.error, width = 1/length(lambda),xlab = "lambda")
    }
    
    list(cv=cv,cv.error=cv.error,finalfit=finalfit,all.folds=all.folds,lambda = object$lambda)
}
