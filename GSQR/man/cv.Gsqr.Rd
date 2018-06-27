
\name{cv.Gsqr}
\alias{cv.Gsqr}

\title{Cross-validation for Gsqr}
\description{Does k-fold cross-validation for Gsqr.}


\usage{
cv.Gsqr(x, y, group, Kfold = 5,delta, taux = 0.5, check = c("f1", "f2"), gamm = ifelse(method == "GMcp", 3, 4), method = c("GLasso", "SGLasso", "GMcp", "GScad", "AGMcp", "AGScad"), plot.it = T,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
\item{y}{response variable: vector of N observations.}
\item{group}{a vector of consecutive integers describing the grouping of the coefficients (see example below).}
\item{Kfold}{number of folds - default is 5. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.}
\item{delta}{the parameter \eqn{\delta}{delta}  in the approximation check function . Default is 2.}
\item{taux}{quantile to be targeted. Must be in (0,1).}
\item{check}{a character string specifying the ckeck function approximation to use, valid options are:
		\itemize{
		\item \code{"f1"}: the check function 1,
		\item \code{"f2"}: the check function 2.
		} Default is \code{"f1"}.}

\item{gamm}{
regularization parameter for the Group Scad or Group Mcp penalty of the 
coefficients
}

\item{method}{a character string specifying the function penality to use, valid options are:
		\itemize{
		\item \code{"GLasso"}: the group Lasso penality,
		\item \code{"GScad"}: the group SCAD penality,
        \item \code{"GMcp"}: the group MCP penality,
		\item \code{"AGScad"}: the approximate group SCAD penality,
        \item \code{"AGMcp"}: the approximate group MCP penality,
		} Default is GLasso.}

        \item{plot.it}{plot cross-validation the error }
       
\item{\dots}{other arguments that can be passed to gglasso.}

}



\details{
The function runs Gsqr Kfold+1 times; the first to get the lambda sequence, and then the remainder to compute the fit with each of the folds omitted. The average error and standard deviation over the folds are computed.
}

\value{an object of class \code{cv.Gsqr} is returned, which is a
list with the ingredients of the cross-validation fit.

		\item{cv}{the mean cross-validated error - a vector of length}
		\item{cv.error}{estimate of standard error of \code{cv}.}
		\item{finalfit}{An object with S3 class \code{Gsqr.}}
		\item{all.folds}{ The folds used in cross-validation.}
		\item{lambda}{the values of \code{lambda} used in the fits.}
}

\author{Mohamed Ouhourane, Yi Yang, Karim Oualkacha
Maintainer: Mohamed OUHOURANE <mohamed.ouhourane@gmail.com>}
\references{
 Mohamed Ouhourane, Yi Yang, Karim Oulkacha: "Group Penalized Smooth Quantile Regression" 
}

\examples{
data(bardet)
group <- rep(1:20,each=5)
fit <- cv.Gsqr(x=bardet$x,y=bardet$y,group=group,method="GLasso",check="f1",taux=0.5)
}

