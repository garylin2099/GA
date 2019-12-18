currMutationRate <- function(mutationRate, maxMutationRate, iterCounter, maxIter) {
  return(mutationRate + (iterCounter - 1) / maxIter * (maxMutationRate - mutationRate))
}

mutate <- function(childVector, mutationRate) {
  mutateGeneIndex <- runif(length(childVector)) <= mutationRate
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}
