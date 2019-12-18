getObjective <- function(X, y, pool, objectiveFunction, regressionType, nCores) {
  getObjVal <- function(chromo) {
    # safety check, if all genes of this chromosome are FALSE, then regress y on an constant
    if (sum(chromo) == 0) {
      model <- lm(y ~ 1)
    } else {
      model <- glm(cbind(y, X[chromo]), family = regressionType)
    }
    return(objectiveFunction(model))
  }
  if (nCores == 1) {
    return(apply(pool, 1, getObjVal))
  } else {
    return(future.apply::future_apply(pool, 1, getObjVal, future.scheduling = nCores))
  }

}

getFitness <- function(objectiveVal) {
  poolSize <- length(objectiveVal)
  r <- poolSize + 1 - rank(objectiveVal) # let small AIC have higher ranking value
  return(2*r/(poolSize*(poolSize+1)))
}
