% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunFor.R
\name{optimal_size}
\alias{optimal_size}
\title{Selecting an optimal tree size}
\usage{
optimal_size(formula, data, npc = NULL, n_folds = 5)
}
\arguments{
\item{formula}{Formula of model fitted}

\item{data}{All the data}

\item{npc}{Number of PCs}

\item{n_folds}{Number of folds in cross-validation}
}
\value{
An integer for optimal tree size.
}
\description{
This function is to find a optimal tree size
}
\details{
This function is used to determine the optimal size of a tree fit.
See ?FunFor, the main function, for example.
}
\examples{
optimal_size(Y ~ ., data)
}
