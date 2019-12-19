#' Selecting Parents by Rank
#'
#' \code{rankSelect} selects parents based on their fitness scores, which in turn depend on their rank in objective function values. In this operation, users have two choices in selecting parents as specified by the parameter "OneParenRandom". The default approach is to select each parent with probability equal to fitness, this approach is valid since the fitness score of all chromosome sums to one, indicating a selection probability. Another approach is to select a parent with probability proportional to fitness and to select the other parent completely at random.
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
