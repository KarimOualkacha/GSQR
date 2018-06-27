\name{coef.cv.Gsqr}
\alias{coef.cv.Gsqr}

\title{get coefficients or make coefficient predictions from a "cv.Gsqr" object.}
\description{
This function gets coefficients or makes coefficient predictions from a cross-validated \code{Gsqr} model,
using the stored \code{"Gsqr"} object, and the optimal value
chosen for \code{lambda}.
}

\arguments{
	\item{cv}{fitted \code{\link{cv.Gsqr}} object.}
	\item{s}{value(s) of the penalty parameter \code{lambda} at which
	predictions are required. Default is the value \code{s="lambda.1se"} stored
	on the CV \code{object}, it is the largest value of \code{lambda} such that error is
	within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be
	used, it is the optimal value of \code{lambda} that gives minimum
	cross validation error \code{cvm}.}
}



\details{This function makes it easier to use the results of cross-validation to get coefficients or make coefficient predictions.}


\author{Mohamed Ouhourane, Yi Yang, Karim Oualkacha
Maintainer: Mohamed OUHOURANE <mohamed.ouhourane@gmail.com>}
\references{
 Mohamed Ouhourane, Yi Yang, Karim Oulkacha: "Group Penalized Smooth Quantile Regression" 
}

\examples{
data(bardet)
group <- rep(1:20,each=5)
cv <- cv.Gsqr(x=bardet$x,y=bardet$y,group=group,method="GLasso",check="f1",taux=0.5)
coef.cv.Gsqr(cv,s = "lambda.1se")
}

