% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kozak_vol.R
\name{kozak_vol}
\alias{kozak_vol}
\title{Estimate the total or partial volume of the tree, based on a fitted Kozak (2004) taper function.}
\usage{
kozak_vol(dbh, h, coef, p, hi, h0)
}
\arguments{
\item{dbh}{tree diameter at breast height, in centimeters.}

\item{h}{total tree height, in meters.}

\item{coef}{numerical vector containing eight coefficients of the Kozak taper equation.}

\item{p}{numerical value representing the first inflection point calculated in the segmented model of Max and Burkhart (1976).}

\item{hi}{final height of the tree section whose volume will be calculated, in meters. Default is the total tree height (h).}

\item{h0}{initial height of the tree section whose volume will be calculated, in meters. Default is 0 (ground height).}
}
\value{
a numeric value indicating the total or partial volume of the tree.
}
\description{
Estimates the total or partial volume of the tree from the diameter at breast height, total height, initial section height, final section height and coefficients of the Kozak (2004) taper equation.
}
\details{
the Kozak (2004) variable-form taper function is represented mathematically by the following expression

di ~ b0*(dbh^b1)\emph{(h^b2)}((1-(hi/h)^(1/4))/(1-(p^(1/3))))^(b3*(hi/h)^4+b4*(1/exp(dbh/h))+b5*((1-(hi/h)^(1/4))/(1-(p^(1/3))))^0.1+b6*(1/dbh)+b7*(h^(1-(hi/h)^(1/3)))+b8*((1-(hi/h)^(1/4))/(1-(p^(1/3)))))
}
\examples{

library(dplyr)
library(minpack.lm)
library(timbeR)

tree_scaling <- tree_scaling \%>\%
mutate(did = di/dbh,
       hih = hi/h)

kozak <- nlsLM(di ~ taper_kozak(dbh, h, hih, b0, b1, b2, b3, b4, b5, b6, b7, b8, p),
               start=list(b0=1.00,b1=.97,b2=.03,b3=.49,b4=-
                            0.87,b5=0.50,b6=3.88,b7=0.03,b8=-0.19, p = .1),
               data = tree_scaling,
               control = nls.lm.control(maxiter = 1000, maxfev = 2000)
)

coef_kozak <- coef(kozak)[-10]
p_kozak <- coef(kozak)[10]

h <- 20
dbh <- 25
di <- 5

kozak_vol(dbh, h, coef_kozak, p_kozak)

hi = 15
h0 = .2

kozak_vol(dbh, h, coef_kozak, p_kozak, hi, h0)

}
