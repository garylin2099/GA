#' Updating Chromosome Pool
#'
#' \code{updatePool} is used to update the chromosome pool in each iteration of the genetic algorithms. It returns a vector of length P*C, which is loaded by all genes of next generation chromosomes.
#'
#' @param pool a matrix of booleans representing the chromosome pools
#' @param fitness a vector indicating the fitness score of each chromosome in the pool
#' @param tournamentSelection a logical indicating whether to use tournament selection. If FALSE, rank-based selection is used. Default is FALSE.
#' @param groupNum number of groups to be partitioned in each round of the tournament selection. Value must be between 2 and poolSize / 2 (inclusive) to ensure meaningful tournaments
#' @param oneParentRandom a logical applying to rank-based selection. If TRUE, then we select one parent with the probability being its fitness score and the other randomly, otherwise, select both parents according to their fitness score. Default is FALSE.
#' @param numCrossoverSplit number of crossover points.
#' @param mutationRate constant mutation rate. Value must be between 0 and 1. Default is 0.01.
#' @param maxMutationRate the mutation rate expected by the user when iteration limit is reached. If specified, then the mutation rate is linearly increasing in each iteration
#' @param iterCounter an integer specifying which iteration in the  the current update
#' @param maxIter maximum number of iterations of updating the generation. Default is 100.
#'
#' @export

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


