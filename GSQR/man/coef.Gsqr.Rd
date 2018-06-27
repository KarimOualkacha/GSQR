\name{coef.Gsqr}
\alias{coef.Gsqr}
\title{get coefficients or make coefficient predictions from an "Gsqr" object.}
\description{
Computes the coefficients at the requested values for \code{lambda} from a fitted \code{\link{Gsqr}} object.
}
\usage{
coef.Gsqr(object, s = NULL, ...)
}

\arguments{

		\item{object}{fitted \code{\link{Gsqr}} model object.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model.}
		\item{\dots}{not used. Other arguments to predict. } 


}


\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. 
}
\value{The coefficients at the requested values for \code{lambda}.}


\author{Mohamed Ouhourane, Yi Yang, Karim Oualkacha
Maintainer: Mohamed OUHOURANE <mohamed.ouhourane@gmail.com>}
\references{
 Mohamed Ouhourane, Yi Yang, Karim Oulkacha: "Group Penalized Smooth Quantile Regression" 
}

\examples{
data(bardet)
group <- rep(1:20,each=5)
fit <- Gsqr(x=bardet$x,y=bardet$y,group=group,method="GLasso",check="f1",taux=0.5)
coef(fit)
}
