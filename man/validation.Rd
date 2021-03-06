% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{check_disjoint}
\alias{check_disjoint}
\alias{check_lt}
\title{Object validation for nascentrna package.}
\usage{
check_disjoint(gr)

check_lt(lesser, greater)
}
\arguments{
\item{gr}{GRanges object to validate.}

\item{lesser}{Numeric value that should be less than \code{greater}.}

\item{greater}{Numeric value used to validate \code{lesser}.}
}
\value{
\code{\link{check_disjoint}} invisibly returns TRUE if the input is
    disjoint.

\code{\link{check_lt}} invisibly returns TRUE if \code{lesser} is
    strictly less than \code{greater}
}
\description{
\code{\link{check_disjoint}} throws an error if the input has overlapping
regions.

\code{\link{check_lt}} throws an error if the \code{lesser} argument is equal
to or exceeds the \code{greater} argument.
}
