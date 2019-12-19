#' Computing Objective Function Values
#'
#' \code{getObjective} returns a poolSize by 1 vector indicating the objective function value for each chromosome in the pool
#'
#' @param X design matrix, in a dataframe or matrix
#' @param y a vector of responses
#' @param pool a matrix of booleans representing the chromosome pools
#' @param objectiveFunction accepts a function object to specify objective function. If provided by user, the function must be able to take glm object as input and return a numeric scalar.
#' @param regressionType a character string defining the distribution family.
#' @param nCores number of cores to use in evaluating the objective function value for each chromosome. If value is greater than 1, the evaluation is executed in parallel.
#'
#' @export
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

#' Computing Fitness Scores
#'
#' \code{getFitness} computes the fitness score of chromosomes based on their objective function values
#'
#' @param objectiveVal a vector of objective function values
#'
#' @export
getFitness <- function(objectiveVal) {
  poolSize <- length(objectiveVal)
  r <- poolSize + 1 - rank(objectiveVal) # let small AICs have higher ranking value
  return(2*r/(poolSize*(poolSize+1)))
}
