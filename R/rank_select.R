#' Selecting Parents by Rank
#'
#' \code{rankSelect} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#'
#' @param pool a matrix of booleans representing the chromosome pools
#' @param fitness a vector indicating the fitness score of each chromosome in the pool
#' @param oneParentRandom a logical applying to rank-based selection. If TRUE, then we select one parent with the probability being its fitness score and the other randomly, otherwise, select both parents according to their fitness score.
#' @param numCrossoverSplit number of crossover points.
#'
#' @export

rankSelect <- function(pool, fitness, oneParentRandom, numCrossoverSplit) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  childVector <- rep(NA, poolSize*chromoSize)
  childCount <- 0
  while(childCount < poolSize) {
    if (oneParentRandom) {
      # select one parent with probability equal to its fitness and the other randomly
      parents <- rbind(pool[sample(1:poolSize, 1, prob = fitness),],
                       pool[sample(1:poolSize, 1),])
    } else {
      # select both parents with probability equal to their fitness
      parents <- pool[sample(1:poolSize, 2, prob = fitness),]
    }
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
