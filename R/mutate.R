#' Mutating Genes
#'
#' \code{crossover} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#' @param parents a 2 by C matrix of booleans, with each row representing a parent chromosome
#'
#' @export
currMutationRate <- function(mutationRate, maxMutationRate, iterCounter, maxIter) {
  return(mutationRate + (iterCounter - 1) / maxIter * (maxMutationRate - mutationRate))
}

#' Single Crossover
#'
#' \code{crossover} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#' @param parents a 2 by C matrix of booleans, with each row representing a parent chromosome
#'
#' @export
mutate <- function(childVector, mutationRate) {
  mutateGeneIndex <- runif(length(childVector)) <= mutationRate
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}
