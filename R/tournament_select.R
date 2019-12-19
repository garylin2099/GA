#' Selecting Parents by Tournament
#'
#' \code{tournamentSelect} selects parents by tournaments. The chromosomes in the parent generation is randomly paritioned into k disjoint subsets of equal size and perform a tournament amongst them. The individual with best objective function (e.g. least AIC score) in each group is chosen as a parent for the next generation. The function takes "groupNum" as one of the parameters, user may use the default value for which the "groupNum' is set to be the value that ensure each group size equals to 2, or user may define the number of subset they'd like to obtain. This process is repeated until the desired amount of population has been reached. Parents are then paired randomly for breeding.
#'
#' @param pool a matrix of booleans representing the chromosome pools
#' @param fitness a vector indicating the fitness score of each chromosome in the pool
#' @param numCrossoverSplit number of crossover points.
#' @param groupNum number of groups to be partitioned in each round of the tournament selection. Value must be between 2 and poolSize / 2 (inclusive) to ensure meaningful tournaments
#'
#'
#' @export

tournamentSelect <- function(pool, fitness, numCrossoverSplit, groupNum) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  # initialize a winner pool to contain winner from each group in each tournament
  winners <- matrix(rep(NA, poolSize * chromoSize), nrow = poolSize)
  winnerCounter <- 1
  # each while loop is a tournament, perform tournaments until the total number of winners generated equals poolSize
  while (winnerCounter <= poolSize) {
    groupSize <- floor(poolSize / groupNum)
    # group assignments using a random permutation
    groupAssignment <- sample(1:poolSize, groupNum * groupSize)
    # each for loop is a within-group competion and produce a winner for the winner pool
    for (i in 1:groupNum) {
      # get the original indices of group members from the generation pool, so as to get their fitness score
      indexInPool <- groupAssignment[((i - 1) * groupSize + 1):(i * groupSize)]
      # select chromosome with the highest fitness score as the winner, here we just borrow rank-based fitness score for use.
      # strictly speaking, tournament selection should be based on the rank of objective function values,
      # however, since rank-based fitness score is proportional to the rank, selection based on rank-based fitness here is also valid
      maxIndexInPool <- indexInPool[which.max(fitness[indexInPool])]
      winners[winnerCounter,] <- pool[maxIndexInPool,]
      winnerCounter = winnerCounter + 1
      if (winnerCounter > poolSize) {
        break
      }
    }
  }
  # after obtaining winners as parents, breeding childs
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
