#' Computing Current Mutation Rate
#'
#' \code{currMutationRate} computes the mutation rate in the current iteration of the genetic algorithm, assuming mutation rate is variable
#'
#' @param mutationRate constant mutation rate. Value must be between 0 and 1.
#' @param maxMutationRate the mutation rate expected by the user when iteration limit is reached. If specified, then the mutation rate is linearly increasing in each iteration.
#' @param iterCounter an integer specifying which iteration in the  the current update
#' @param maxIter maximum number of iterations of updating the generation.
#'
#' @export
currMutationRate <- function(mutationRate, maxMutationRate, iterCounter, maxIter) {
  return(mutationRate + iterCounter / maxIter * (maxMutationRate - mutationRate))
}

#' Mutating Genes
#'
#' \code{mutate} takes in a vector of genes and try mutating each of them independently with a mutation rate
#'
#' @param childVector a vector representing the genes of child chromosomes
#' @param mutationRate mutation rate
#'
#' @export
mutate <- function(childVector, mutationRate) {
  mutateGeneIndex <- runif(length(childVector)) <= mutationRate
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}
