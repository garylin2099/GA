% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select}
\alias{select}
\title{Variable Selection using Genetic Algorithms}
\usage{
select(
  X,
  y,
  poolSize = round(2 * ncol(X)),
  regressionType = "gaussian",
  objectiveFunction = stats::AIC,
  oneParentRandom = FALSE,
  tournamentSelection = FALSE,
  groupNum = NULL,
  numCrossoverSplit = 1,
  mutationRate = 0.01,
  maxMutationRate = NULL,
  maxIter = 100,
  minIter = round(maxIter/10),
  diversityCutoff = max(0.05, 1/poolSize),
  objValPlot = TRUE,
  nCores = 1
)
}
\arguments{
\item{X}{design matrix, in a dataframe or matrix}

\item{y}{a vector of responses}

\item{poolSize}{the size of the chromosome pool for any generation, i.e. P. Default is 2C where C, or chromosome length, is the number of candidate covariates}

\item{regressionType}{a character string defining the distribution family. Default is gaussian.}

\item{objectiveFunction}{accepts a function object to specify objective function. If provided by user, the function must be able to take glm object as input and return a numeric scalar. Default is \code{\link{AIC}}.}

\item{oneParentRandom}{a logical applying to rank-based selection. If TRUE, then we select one parent with the probability being its fitness score and the other randomly, otherwise, select both parents according to their fitness score. Default is FALSE.}

\item{tournamentSelection}{a logical indicating whether to use tournament selection. If FALSE, rank-based selection is used. Default is FALSE.}

\item{groupNum}{number of groups to be partitioned in each round of the tournament selection. Value must be between 2 and poolSize / 2 (inclusive) to ensure meaningful tournaments}

\item{numCrossoverSplit}{number of crossover points.}

\item{mutationRate}{constant mutation rate. Value must be between 0 and 1. Default is 0.01.}

\item{maxMutationRate}{the mutation rate expected by the user when iteration limit is reached. If specified, then the mutation rate is linearly increasing in each iteration}

\item{maxIter}{maximum number of iterations of updating the generation. Default is 100.}

\item{minIter}{minimum number of iterations, specified to reduce the chance of pre-mature convergence. Default is one fifth of the maxIter}

\item{diversityCutoff}{the cutoff of diversity level for claiming convergence and stoping the iteration, where diversity level is defined as the number of unique chromosomes in a pool devided by the pool size. User can use a non-positive value to avoid convergence check. Default is the larger one between 0.05 and 1 / poolSize.}

\item{objValPlot}{a logical indicating whether to plot objective function values over generations. Default is true.}

\item{nCores}{number of cores to use in evaluating the objective function value for each chromosome. If value is greater than 1, the evaluation is executed in parallel. Default is 1.}
}
\description{
\code{select} is the main function of our genetic algorithms, it's a search heuristic that mimics evolution to find the best solution through variable selection.
A population of randomly generated candidate solutions is initialized, the candidate solutions are then evaluated, and their fitness value is calculated based on AIC scores or other user specific objective functions. The populations then goes through an iterative process of selection, crossover, mutation to produce the next generation until convergence.
}
\details{
1. Geof H. Givens, Jennifer A. Hoeting (2013) Combinatorial Optimization. Chapter 3 of Computational Statistics.
}
\examples{
set.seed(2019)
x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
x4 <- rnorm(100, 18, 1)
x5 <- rnorm(100, -7, 1)
x6 <- rnorm(100, -20, 1)
x7 <- rnorm(100, 15, 1)
x8 <- rnorm(100, -9, 1)
X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
y <- x1 + x3 + x4 - 2 * x6 - x8 + rnorm(100, 0, 0.1)
result1 <- select(X, y)
result2 <- select(X, y, tournamentSelection = TRUE, numCrossoverSplit = 2, maxMutationRate = 0.05)
}
