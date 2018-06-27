\name{Gsqr}
\alias{Gsqr}

\title{Fits the regularization paths for group penalized quantile regression}


\description{Fits a regularization path for the group penalized quantile regression at a sequence of regularization parameters lambda.}

\usage{
Gsqr(x, y, group = NULL, method = c("GLasso", "SGLasso", "GMcp", "GScad", "AGMcp", "AGScad"), check = c("f1", "f2"), nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.05, 0.001), lambda = NULL, pf = sqrt(bs), dfmax = as.integer(max(group)) + 1, pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08, gamm = ifelse(((method == "GMcp") | (method == "AGMcp")), 3, 4), taux = 0.5, delta=2, intercept = TRUE)
}


\arguments{

\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
\item{y}{response variable: vector of N observations.}
\item{group}{a vector of consecutive integers describing the grouping of the coefficients (see example below).}
\item{method}{a character string specifying the function penality to use, valid options are:
		\itemize{
		\item \code{"GLasso"}: the group Lasso penality,
		\item \code{"GScad"}: the group SCAD penality,
        \item \code{"GMcp"}: the group MCP penality,
		\item \code{"AGScad"}: the approximate group SCAD penality,
        \item \code{"AGMcp"}: the approximate group MCP penality,
		} Default is \code{"GLasso"}.}

\item{check}{a character string specifying the ckeck function approximation to use, valid options are:
		\itemize{
		\item \code{"f1"}: the check function 1,
		\item \code{"f2"}: the check function 2.
		} Default is "f1".}

\item{nlambda}{the number of \code{lambda} values - default is 100.}
\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 3e8. If models do not converge, consider increasing \code{maxit}.}

\item{lambda.factor}{The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined lambda sequence.} 

\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

\item{pf}{penalty factor, a vector in length of bn (bn is the total number of groups). Separate penalty weights can be applied to each group of \eqn{\beta}{beta's}s to allow}
\item{dfmax}{limit the maximum number of groups in the model. Useful for very large \code{bs} (group size), if a partial path is desired. Default is \code{bs+1}.}
\item{pmax}{limit the maximum number of groups ever to be nonzero. For example once a group enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,bs)}.}

\item{eps}{convergence termination tolerance. Defaults value is \code{1e-8}.}

\item{gamm}{ regularization parameter for the Group Scad or Group Mcp penalty of the coefficients}

\item{taux}{
   quantile to be targeted. Must be in (0,1).
}

\item{delta}{the parameter \eqn{\delta}{delta}  in the approximation check function . Default is 2.}

\item{intercept}{Whether to include intercept in the model. Default is TRUE.}	

}


\value{

		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix class}
		\item{lambda}{the actual sequence of \code{lambda} values used}
}


\author{Mohamed Ouhourane, Yi Yang, Karim Oualkacha
Maintainer: Mohamed OUHOURANE <mohamed.ouhourane@gmail.com>}
\references{
 Mohamed Ouhourane, Yi Yang, Karim Oulkacha: "Group Penalized Smooth Quantile Regression" 
}

\examples{
data(bardet)
group <- rep(1:20,each=5)
fit <- Gsqr(x=bardet$x,y=bardet$y,group=group,method="GLasso",check="f1",taux=0.5)
plot(fit)
}

