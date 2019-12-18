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
