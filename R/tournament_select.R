#' Single Crossover
#'
#' \code{crossover} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#' @param parents a 2 by C matrix of booleans, with each row representing a parent chromosome
#'
#' @export
tournamentSelect <- function(pool, fitness, numCrossoverSplit, groupNum) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  winners <- matrix(rep(NA, poolSize * chromoSize), nrow = poolSize)
  winnerCounter <- 1
  while (winnerCounter <= poolSize) {
    groupSize <- floor(poolSize / groupNum)
    groupAssignment <- sample(1:poolSize, groupNum * groupSize)
    for (i in 1:groupNum) {
      indexInPool <- groupAssignment[((i - 1) * groupSize + 1):(i * groupSize)]
      maxIndexInPool <- indexInPool[which.max(fitness[indexInPool])]
      winners[winnerCounter,] <- pool[maxIndexInPool,]
      winnerCounter = winnerCounter + 1
      if (winnerCounter > poolSize) {
        break
      }
    }
  }
  childVector <- rep(NA, poolSize*chromoSize)
  childCount <- 0
  pairAssignment <- sample(1:poolSize, poolSize)
  while(childCount < poolSize) {
    parents <- winners[pairAssignment[(childCount+1):(childCount+2)],]
    if (numCrossoverSplit == 1) {
      childs <- crossover(parents)
    } else {
      childs <- multipleCrossover(parents, numCrossoverSplit)
    }
    childVector[(childCount*chromoSize) + 1:(2*chromoSize)] <- childs
    childCount = childCount + 2
  }
  return(childVector)
}
