% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poly5_logs.R
\name{poly5_logs}
\alias{poly5_logs}
\title{Simulate log extraction using a 5th degree polynomial that describes the taper of the tree.}
\usage{
poly5_logs(
  dbh,
  h,
  coef,
  assortments,
  stump_height,
  downgrade,
  broken,
  defect_height,
  eliminate,
  total_volume,
  only_vol
)
}
\arguments{
\item{dbh}{tree diameter at breast height, in centimeters.}

\item{h}{total tree height, in meters.}

\item{coef}{numerical vector containing six coefficients of the 5th degree polynomial function that describes the tree's taper.}

\item{assortments}{a data.frame with five columns and n rows, where n is the number of different wood assortments to be obtained from the tree stem. The first column must contain the names of the assortments, the second, numerical, contains the minimum diameters at the small end of the logs, in centimeters. The third column, numerical, contains the minimum lengths of the logs, in meters. The fourth column, numerical, contains the maximum lengths of the logs, in meters. The fifth column, numerical, contains the values in centimeters referring to the loss of wood due to cutting logs. The algorithm prioritizes the extraction of assortments along the stem in the order presented in the data.frame, starting from the first line, to the last.}

\item{stump_height}{tree cutting height, in meters. Default is 0.}

\item{downgrade}{if TRUE, the algorithm,from the defect_height onwards, simulates log extraction only for the last assortment in the assortments data.frame. Default is FALSE.}

\item{broken}{if TRUE, the algorithm will simulate the extraction of logs only up to the defect_height. Default is FALSE.}

\item{defect_height}{the height, in meters, from which the logs will be downgraded (if downgrade is TRUE) or log extraction simulation will be stopped (if broken is TRUE). Default is 0 for downgrade = TRUE (the whole tree is downgraded) and h * 0.5 for broken = TRUE (the tree is broken from half its original/estimated total height).}

\item{eliminate}{if TRUE, the algorithm does not get logs for any assortment present in the assortments table. All will be zero. Default is FALSE.}

\item{total_volume}{if TRUE, it adds an additional column to the results data.frame with the estimate of the total volume of the tree, from the ground height to h if broken argument is FALSE, or to defect_height if broken is TRUE. Default is FALSE.}

\item{only_vol}{if TRUE returns only volumes (does not return the number of logs). Default is FALSE.}
}
\value{
a list of two data.frames, the first (volumes) with the calculated volumes per assortment, and the second (logs) with the number of logs per assortment.
}
\description{
Simulate the extraction of logs from a tree from its measurements, taper function (5th degree polynomial), trunk quality characteristics and harvest parameters such as stump height and assortments.
}
\details{
when the \code{broken} and \code{downgrade} arguments are set to TRUE, the \code{defect_height} value is considered as the break height of the tree, and the entire tree is downgraded.
}
\examples{

library(dplyr)
library(minpack.lm)
library(timbeR)

tree_scaling <- tree_scaling \%>\%
mutate(did = di/dbh,
       hih = hi/h)

poli5 <- lm(did~hih+I(hih^2)+I(hih^3)+I(hih^4)+I(hih^5),tree_scaling)

coef_poli <- coef(poli5)

dbh <- 25
h <- 20

assortments <- data.frame(
  NAME = c('15-25','4-15'),
  SED = c(15,4),
  MINLENGTH = c(2.65,2),
  MAXLENGTH = c(2.65,4.2),
  LOSS = c(5,5)
)

poly5_logs(dbh, h, coef_poli, assortments)

}
