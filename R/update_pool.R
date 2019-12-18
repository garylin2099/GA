
updatePool <- function(pool, fitness, tournamentSelection, groupNum, oneParentRandom, numCrossoverSplit,
                       mutationRate, maxMutationRate, iterCounter, maxIter) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  if (!tournamentSelection) {
    childVector <- rankSelect(pool, fitness, oneParentRandom, numCrossoverSplit)
  } else {
    childVector <- tournamentSelect(pool, fitness, numCrossoverSplit, groupNum)
  }

  if (is.null(maxMutationRate)) {
    childVector <- mutate(childVector, mutationRate)
  } else {
    childVector <- mutate(childVector, currMutationRate(mutationRate, maxMutationRate, iterCounter, maxIter))
  }

  return(matrix(childVector[1:(poolSize*chromoSize)], nrow = poolSize, byrow = TRUE))
}


