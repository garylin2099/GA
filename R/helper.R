###################
# Helper  Functions
###################

checkInputValidity <- function(X,y,poolSize,regressionType,objectiveFunction,oneParentRandom,
                               tournamentSelection,groupNum,numCrossoverSplit,mutationRate,
                               maxMutationRate,maxIter,minIter,diversityCutoff,nCores) {
  if (!is.data.frame(X)) stop("X must be a dataframe or matrix")
  if (!is.vector(y)) {
    if (!is.data.frame(y) || ncol(y) != 1) {
      stop("y must be a vector or a dataframe of just one column")
    }
  }
  if (!is.numeric(poolSize)) stop("poolSize must be numeric")
  if (!is.character(regressionType)
      || !(regressionType %in% c("binomial","gaussian","Gamma","inverse.gaussian",
                                 "poisson","quasi","quasibinomial","quasipoisson"))) {
    stop("must provide valid ditribution family name")
  }
  if (!is.function(objectiveFunction)) stop("objectiveFunction must be a function")
  if (!is.logical(oneParentRandom)) stop("oneParentRandom argument must be a logical variable")
  if (oneParentRandom && tournamentSelection) stop("no need to specify oneParentRandom when using tournament selection")
  if (!is.logical(tournamentSelection)) stop("tournamentSelection argument must be a logical variable")
  if (!tournamentSelection) {
    if (!is.null(groupNum)) stop("no need to specify groupNum argument if tournament selection is not adopted")
  } else {
    if (!is.null(groupNum)) {
      if (!is.numeric(groupNum)) stop("groupNum must be numeric")
      if (groupNum < 2 || groupNum > poolSize / 2) stop("groupNum must be between 2 and poolSize / 2,
                                                        otherwise each group will contain all chromosomes or
                                                        only one chromosome, making tournaments meaningless")
    } else {
      if (poolSize < 2 * 2) {
        stop("poolSize too small to guarantee at least two tournament groups with two candidates per group")
      }
    }
  }
  if (!is.numeric(numCrossoverSplit)) stop("numCrossOverSplit must be numeric")
  if (numCrossoverSplit < 1 || numCrossoverSplit > ncol(X) - 1) stop("numCrossoverSplit must range
                                                                     from 1 to the number of features
                                                                     minus 1 (inclusively)")
  if (!is.numeric(mutationRate)) stop("mutationRate must be numeric")
  if (mutationRate < 0 || mutationRate > 1) stop("mutationRate must be between 0 and 1")
  if (!is.null(maxMutationRate)) {
    if (!is.numeric(maxMutationRate)) stop("maxMutationRate must be numeric")
    if (maxMutationRate < mutationRate) stop("maxMutationRate must be no less than mutationRate")
    if (maxMutationRate > 1) stop("maxMutationRate must be no more than 1")
  }
  if (!is.numeric(maxIter)) stop("maxIter must be numeric")
  if (maxIter <= 0) stop("maxIter must be positive")
  if (!is.numeric(minIter)) stop("minIter must be numeric")
  if (minIter <= 0) stop("minIter must be positive")
  if (maxIter < minIter) stop("maxIter must be no less than minIter")
  if (!is.numeric(diversityCutoff)) stop("diversityCutoff must be numeric")
  # user can use negative value or 0 to keep updating generation for maxIter iterations
  if (diversityCutoff > 1) stop("diversityCutoff must be no more than 1")
  if (!is.numeric(nCores)) stop("nCores must be numeric")
  if (nCores < 1) stop("nCores must have positive")
}

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

plotObjVal <- function(iter, objValue) {
  plot(iter, objValue + runif(length(objValue), -0.2 * min(abs(objValue)), 0.2 * min(abs(objValue))), # we add perturbation to help visualization
       main = "Objective Function Values (e.g. AIC) w Perturbation",
       xlab = "Generation i", ylab = "Objective Function Values",
       ylim = 1.5 * c(min(objValue), max(objValue)))
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

