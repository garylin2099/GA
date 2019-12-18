#' Variable Selection using Genetic Algorithms
#'
#' select implements genetic algorithms for variable selection for GLMs by optimizing package or user specified objective functions such as AIC, BIC, and logloglikelihood.
#' Uses functions: \code{\link{generate_founders}}, \code{\link{evaluate_fitness}}, and \code{\link{create_next_generation}}.
#' Functions find optimal variables by using evolutationry biology concepts of natural selection, fitness, genetic crossover, and mutation. Founding generation of chromosomes is randomly generated and evaluated using an critieria such as AIC, BIC, or loglihood. Parents are selected by their fitness, and generate children chromosomes. As each generation breeds and produces new genreations, the algorithm moves towards the optimum.
#'
#' 1. Geof H. Givens, Jennifer A. Hoeting (2013) Combinatorial Optimization (italicize). Chapter 3 of Computational Statistics (italicize).
#'
#' @param Y vector of response variable
#' @param X a matrix or dataframe of predictor variables
#' @param family a character string describing the error distribution and link function to be used in the model. Default is gaussian.
#' @param objective_function function for computing objective. Default is \code{\link{AIC}}. User can specify custom function.
#' @param crossover_parents_function a function for crossover between mate pairs. User can specify custom function. Default is \code{\link{crossover_parents}}.
#' @param crossover_method a character string describing crossover method. Default is multi-point crossover. See \code{\link{crossover_parents}}.
#' @param pCrossover a numeric value for he probability of crossover for each mate pair.
#' @param start_chrom a numeric value for the  size of the popuation of chromosomes. Default is \code{choose(C, 2)} \eqn{\le 200}, where C is number of predictors.
#' @param mutation_rate a numeric value for rate of mutation. Default is \eqn{1 / (P \sqrt C)}, where P is number of chromosomes, and C is number of predictors.
#' @param converge a logical value indicating whether algorithm should attempt to converge or run for specified number of iterations. If \code{TRUE}, convergence will occur when highest ranked chromosomes is equal to mean of top 50\% in current and previous generation.
#' @param tol a numeric value indicating convergence tolerance. Default is 1e-4.
#' @param iter an integer specifying maximum number of generations algorithm will produce. Default is 100
#' @param minimize a logical value indicating whether optimize should be minimized (TRUE) or maximized (FALSE).
#' @param nCores an integer indicating number of parallel processes to run when evaluating fitness. Default is 1, or no paralleization. See \code{\link{evaluate_fitness}}.
#'
#'If user wants to use custom objective_function, they must use a function that is compatible with \code{\link{lm}} or \code{\link{glm}} fitted objects which output a numberic value of length 1.
#'
#' @examples
#'x1 <- rnorm(100, 2, 1)
#'x2 <- rnorm(100, 5, 1)
#'x3 <- rnorm(100, -1, 1)
#'x4 <- rnorm(100, 18, 1)
#'x5 <- rnorm(100, -7, 1)
#'x6 <- rnorm(100, -20, 1)
#'x7 <- rnorm(100, 15, 1)
#'x8 <- rnorm(100, -9, 1)
#'X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
#'y <- x1 + x3 + x4 - 2 * x6 - x8 + rnorm(100, 0, 0.1)
#'lm(cbind(y, X))
#'glm(cbind(y, X), family = "gaussian")

#'for (i in 1:1) {
#'  result <- select(X, y, maxMutationRate = 0.05, numCrossoverSplit = 3, tournamentSelection = FALSE)
#'}
#'result



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
                   nCores = 1) {
  # check if arguments are valid
  checkInputValidity(X,y,poolSize,regressionType,objectiveFunction,oneParentRandom,
                     tournamentSelection,groupNum,numCrossoverSplit,mutationRate,
                     maxMutationRate,maxIter,minIter,diversityCutoff,nCores)

  allData <- cbind(y, X)
  chromoSize <- ncol(X)
  pool <- init(poolSize, chromoSize)
  objValEachIter <- data.frame(
    iter = rep(NA, poolSize * maxIter),
    objectiveValue = rep(NA, poolSize * maxIter)
  ) # here
  for (i in 1:maxIter) {
    objVal <- getObjective(X, y, pool, objectiveFunction, regressionType, nCores)
    fitness <- getFitness(objVal)
    # record objective values (by default AIC) of each generation
    objValEachIter[((i - 1) * poolSize + 1):(i * poolSize),] <- cbind(rep(i, poolSize), objVal)
    pool <- updatePool(pool, fitness, tournamentSelection, groupNum, oneParentRandom, numCrossoverSplit, mutationRate, maxMutationRate, i, maxIter)
    if (convergeCheck(pool, i, minIter, diversityCutoff)) {
      cat("number of iterations to achieve converge is", i, "\n")
      break
    }
  }
  objValEachIter <- objValEachIter[!is.na(objValEachIter$iter),]
  objValue <- objValEachIter$objectiveValue
  plot(objValEachIter$iter, objValue + runif(length(objValue), -0.2 * min(abs(objValue)), 0.2 * min(abs(objValue))), # we add perturbation to help visualization
       main = "Objective Function Values (e.g. AIC) w Perturbation",
       xlab = "Generation i", ylab = "Objective Function Values",
       ylim = 1.5 * c(min(objValue), max(objValue)))
  majorChromo <- getMajorChromo(pool)
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
