% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossover.R
\name{crossover}
\alias{crossover}
\title{Single Crossover}
\usage{
crossover(parents)
}
\arguments{
\item{parents}{a 2 by C matrix of booleans, with each row representing a parent chromosome}
}
\description{
\code{crossover} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
}
