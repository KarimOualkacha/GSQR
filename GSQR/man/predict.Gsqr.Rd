\name{predict.Gsqr}
\alias{predict.Gsqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Outputs predicted responses from an SGL model for new observations}
\description{
Outputs predicted response values for new user input observations at a specified s value}

\usage{
predict.Gsqr(object, newx, s = NULL)
}

\arguments{
\item{object}{fitted \code{"GSQR"} object}
\item{newx}{covariate matrix for new observations whose responses we wish to predict}
\item{s}{the index of the lambda value for the model with which we desire to predict}
}

\value{
prediction vector of news
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
predict(fit,bardet$x)
}
