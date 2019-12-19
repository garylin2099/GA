#' Variable Selection using Genetic Algorithms
#'
#' \code{select} is the main function of our genetic algorithms, it's a search heuristic that mimics evolution to find the best solution through variable selection.
#' A population of randomly generated candidate solutions is initialized, the candidate solutions are then evaluated, and their fitness value is calculated based on AIC scores or other user specific objective functions. The populations then goes through an iterative process of selection, crossover, mutation to produce the next generation until convergence.
#'
#'
#' 1. Geof H. Givens, Jennifer A. Hoeting (2013) Combinatorial Optimization (italicize). Chapter 3 of Computational Statistics (italicize).
#'
#' @param X design matrix, in a dataframe or matrix
#' @param y a vector of responses
#' @param poolSize the size of the chromosome pool for any generation, i.e. P. Default is 2C where C, or chromosome length, is the number of candidate covariates
#' @param regressionType a character string defining the distribution family. Default is gaussian.
#' @param objectiveFunction accepts a function object to specify objective function. If provided by user, the function must be able to take glm object as input and return a numeric scalar. Default is \code{\link{AIC}}.
#' @param oneParentRandom a logical applying to rank-based selection. If TRUE, then we select one parent with the probability being its fitness score and the other randomly, otherwise, select both parents according to their fitness score. Default is FALSE.
#' @param tournamentSelection a logical indicating whether to use tournament selection. If FALSE, rank-based selection is used. Default is FALSE.
#' @param groupNum number of groups to be partitioned in each round of the tournament selection. Value must be between 2 and poolSize / 2 (inclusive) to ensure meaningful tournaments
#' @param numCrossoverSplit number of crossover points.
#' @param mutationRate constant mutation rate. Value must be between 0 and 1. Default is 0.01.
#' @param maxMutationRate the mutation rate expected by the user when iteration limit is reached. If specified, then the mutation rate is linearly increasing in each iteration
#' @param maxIter maximum number of iterations of updating the generation. Default is 100.
#' @param minIter minimum number of iterations, specified to reduce the chance of pre-mature convergence. Default is one fifth of the maxIter
#' @param diversityCutoff the cutoff of diversity level for claiming convergence and stoping the iteration, where diversity level is defined as the number of unique chromosomes in a pool devided by the pool size. User can use a non-positive value to avoid convergence check. Default is the larger one between 0.05 and 1 / poolSize.
#' @param objValPlot a logical indicating whether to plot objective function values over generations. Default is true.
#' @param nCores number of cores to use in evaluating the objective function value for each chromosome. If value is greater than 1, the evaluation is executed in parallel. Default is 1.
#'
#'
#' @examples
#' set.seed(2019)
#' x1 <- rnorm(100, 2, 1)
#' x2 <- rnorm(100, 5, 1)
#' x3 <- rnorm(100, -1, 1)
#' x4 <- rnorm(100, 18, 1)
#' x5 <- rnorm(100, -7, 1)
#' x6 <- rnorm(100, -20, 1)
#' x7 <- rnorm(100, 15, 1)
#' x8 <- rnorm(100, -9, 1)
#' X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
#' y <- x1 + x3 + x4 - 2 * x6 - x8 + rnorm(100, 0, 0.1)
#' result1 <- select(X, y)
#' result2 <- select(X, y, tournamentSelection = TRUE, numCrossoverSplit = 2, maxMutationRate = 0.05)

#' @export
select <- function(X,
                   y,
                   poolSize = round(2 * ncol(X)),
                   regressionType = "gaussian",
                   objectiveFunction = stats::AIC,
                   oneParentRandom = FALSE,
                   tournamentSelection = FALSE,
                   groupNum = NULL,
                   numCrossoverSplit = 1, # need to be smaller than chromosome size
                   mutationRate = 0.01,
                   maxMutationRate = NULL,
                   maxIter = 100,
                   minIter = round(maxIter / 10),
                   diversityCutoff = max(0.05, 1 / poolSize),
                   objValPlot = TRUE,
                   nCores = 1) {
  # individual validity check
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  # main part of the validity checks
  checkInputValidity(X,y,poolSize,regressionType,objectiveFunction,oneParentRandom,
                     tournamentSelection,groupNum,numCrossoverSplit,mutationRate,
                     maxMutationRate,maxIter,minIter,diversityCutoff,nCores)
  # individual default value setup
  if (tournamentSelection && is.null(groupNum)) {
    groupNum <- floor(poolSize / 2) # by default partition pool into groups of two
  }

  # necessary setup including initialization
  allData <- cbind(y, X)
  chromoSize <- ncol(X)
  pool <- init(poolSize, chromoSize)
  objValEachIter <- data.frame(
    iter = rep(NA, poolSize * (maxIter + 1)),
    objectiveValue = rep(NA, poolSize * (maxIter + 1))
  )
  # select, crossover, and mutate to update pool (generation) in each iteration, keep updating it
  # until convergence or maximum number of iteration is hit
  for (i in 1:maxIter) {
    objVal <- getObjective(X, y, pool, objectiveFunction, regressionType, nCores)
    fitness <- getFitness(objVal)
    # a side step, record objective values of each generation for plotting
    objValEachIter[((i - 1) * poolSize + 1):(i * poolSize),] <- cbind(rep(i, poolSize), objVal)
    pool <- updatePool(pool, fitness, tournamentSelection, groupNum, oneParentRandom, numCrossoverSplit,
                       mutationRate, maxMutationRate, i, maxIter)
    if (convergeCheck(pool, i, minIter, diversityCutoff)) {
      objValConverge <- getObjective(X, y, pool, objectiveFunction, regressionType, nCores)
      objValEachIter[(i * poolSize + 1):((i + 1) * poolSize),] <- cbind(rep(i+1, poolSize), objValConverge)
      cat("number of iterations to achieve convergence is", i, "\n")
      break
    }
  }
  # plotting objective function values against number of iteration
  objValEachIter <- objValEachIter[!is.na(objValEachIter$iter),]
  if (objValPlot) {
    plotObjVal(objValEachIter$iter, objValEachIter$objectiveValue)
  }
  # generate final result
  majorChromo <- getMajorChromo(pool) # find the chromosome type that dominates the final pool
  resultModel <- glm(cbind(y, X[majorChromo]), family = regressionType)
  result <- list(
    "selectedVariables" = colnames(X)[majorChromo],
    "model" = resultModel,
    "ObjectiveValue" = objectiveFunction(resultModel),
    "convergence" = convergeCheck(pool, i, minIter, diversityCutoff),
    "iterationCount" = i,
    "objectiveValueRecord" = objValEachIter
  )
  cat("Variable selected are", colnames(X)[majorChromo], "\n")
  return(result)

}
