% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_and_remove.R
\name{select_and_remove}
\alias{select_and_remove}
\title{Remove unwanted data by selecting it}
\usage{
select_and_remove(data, xvar, yvar)
}
\arguments{
\item{data}{a data.frame.}

\item{xvar}{quoted name of the variable to be displayed in the x axis.}

\item{yvar}{quoted name of the variable to be displayed in the y axis.}
}
\value{
the data.frame given to the \code{data} argument, without the selected points.
}
\description{
Delete unwanted records from the dataset (e.g. outliers) by selecting them in a scatter plot.
}
\examples{

\dontrun{
library(dplyr)
library(timbeR)

tree_scaling <- tree_scaling \%>\%
mutate(did = di/dbh,
       hih = hi/h) \%>\%
select_and_remove(., 'hih', 'did')
}

}
