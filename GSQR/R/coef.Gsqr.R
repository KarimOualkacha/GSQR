coef.Gsqr <-
function(object, s = NULL, ...) {
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        if(length(s) == 1)
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
			nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
		} else
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
			nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
		}
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    return(nbeta)
}
