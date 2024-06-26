% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{data_helpers}
\alias{data_helpers}
\alias{number_of_missing}
\alias{is_dichotomous}
\title{Data summarization and classification methods}
\usage{
number_of_missing(x)

is_dichotomous(x)
}
\arguments{
\item{x}{A vector of any of the atomic types (see [\code{\link[base:vector]{base::vector()}}])}
}
\value{
Returns a single value determined by the individual functions
}
\description{
These related functions are intended to analyze a single data vector (e.g.
column from a dataset) and help predict its classification, or other relevant
attributes. These are simple yet opionated convenience functions.
}
\details{
The functions that are currently supported are:
\itemize{
\item \code{number_of_missing()} returns the number of missing values in a vector
\item \code{is_dichotomous()} returns TRUE if the vector is dichotomous, FALSE otherwise
}
}
