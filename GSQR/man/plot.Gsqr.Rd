\name{plot.Gsqr}
\alias{plot.Gsqr}
\title{Plot solution paths from a "Gsqr" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{Gsqr}} object.
}
\usage{
\method{plot}{Gsqr}(x, group = FALSE, log.l = TRUE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{Gsqr}} model}
		\item{group}{what is on the Y-axis. Plot the norm of each group
		 if \code{TRUE}. Plot each coefficient if \code{FALSE}.}
		\item{log.l}{what is on the X-axis. Plot against the log-lambda
		sequence if \code{TRUE}. Plot against the lambda sequence if \code{FALSE}.}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
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



