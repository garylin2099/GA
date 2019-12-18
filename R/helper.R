
init <- function(poolSize, chromSize) {
  return(matrix(sample(c(TRUE, FALSE), chromSize * poolSize, replace = TRUE), nrow = poolSize, byrow = TRUE))
}

convergeCheck <- function(pool, currentIterCount, minIter, diversityCutoff) {
  if (currentIterCount >= minIter & nrow(unique(pool)) / nrow(pool) <= diversityCutoff) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


getMajorChromo <- function(pool) {
  uniqueChromo <- unique(pool)
  numEachChromo <- rep(NA, nrow(uniqueChromo))
  for (i in 1:nrow(uniqueChromo)) {
    numEachChromo[i] <- sum(apply(pool, 1, function(poolRow) {return(all(poolRow == uniqueChromo[i,]))}))
  }
  majorChromoIndex <- which(numEachChromo == max(numEachChromo))[1]
  return(uniqueChromo[majorChromoIndex,])
}

