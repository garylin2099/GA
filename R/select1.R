rm(list = ls())
set.seed(201912)

select <- function(X,
                   y,
                   poolSize = round(2 * ncol(X)),
                   objectiveFunction = stats::AIC,
                   mutationRate = 0.01,
                   mutationRateByEnd = NULL,
                   maxIter = 100,
                   minIter = maxIter / 5,
                   regressionType = "gaussian",
                   oneParentRandom = FALSE,
                   diversityCutoff = 1 / poolSize) {
  allData <- cbind(y, X)
  chromoSize <- ncol(X)
  pool <- init(poolSize, chromoSize)
  for (i in 1:maxIter) {
    objVal <- getObjective(X, y, pool, objectiveFunction, regressionType)
    fitness <- getFitness(objVal)
    pool <- updatePool(pool, fitness, oneParentRandom, mutationRate, i)
    if (convergeCheck(pool, i, minIter, diversityCutoff)) {
      cat("number of iterations to achieve converge is", i, "\n")
      break
    }
  }
  majorChromo <- getMajorChromo(pool)
  cat("Variable selected are", colnames(X)[majorChromo], "\n")
  cat(majorChromo, "\n")
  return(majorChromo)

}
init <- function(poolSize, chromSize) {
  return(matrix(sample(c(TRUE, FALSE), chromSize * poolSize, replace = TRUE), nrow = poolSize, byrow = TRUE))
}


getObjective <- function(X, y, pool, objectiveFunction, regressionType) {
  getObjVal <- function(chromo) {
    # safety check, if all genes of this chromosome are FALSE, then regress y on an constant
    if (sum(chromo) == 0) {
      model <- lm(y ~ 1)
    } else {
      model <- glm(cbind(y, X[chromo]), family = regressionType)
    }
    return(objectiveFunction(model))
  }
  return(apply(pool, 1, getObjVal))
}

getFitness <- function(objectiveVal) {
  poolSize <- length(objectiveVal)
  r <- poolSize + 1 - rank(objectiveVal) # let small AIC have higher ranking value
  return(2*r/(poolSize*(poolSize+1)))
}

updatePool <- function(pool, fitness, oneParentRandom, mutationRate, mutationRateByEnd, iterCounter) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  childVector <- rankSelect(pool, fitness, oneParentRandom)
  if (is.null(mutationRateByEnd)) {
    childVector <- mutate(childVector, mutationRate)
  } else {
    childVector <- mutate(childVector, currMutationRate(mutationRate, mutationRateByEnd, i))
  }

  return(matrix(childVector[1:(poolSize*chromoSize)], nrow = poolSize, byrow = TRUE))
}

rankSelect <- function(pool, fitness, oneParentRandom) {
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
    childs <- crossover(parents)
    childVector[(childCount*chromoSize) + 1:(2*chromoSize)] <- childs
    childCount = childCount + 2
  }
  return(childVector)
}

crossover <- function(parents) {
  chromoSize <- ncol(parents)
  splitPoint <- sample(2:chromoSize, 1)
  firstChild <- c(parents[1,1:(splitPoint-1)], parents[2,splitPoint:chromoSize])
  secondChild <- c(parents[2,1:(splitPoint-1)], parents[1,splitPoint:chromoSize])
  return(c(firstChild, secondChild))
}

mutate <- function(childVector, mutationRate) {
  mutateGeneIndex <- runif(length(childVector)) <= mutationRate
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}

currMutationRate <- function(mutationRate, mutationRateByEnd, iterCounter) {

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
  majorChromoIndex <- which(numEachChromo == max(numEachChromo))
  return(uniqueChromo[majorChromoIndex,])
}

#test data
x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
x4 <- rnorm(100, 18, 1)
x5 <- rnorm(100, -7, 1)
X <- data.frame(x1, x2, x3, x4, x5)
y <- 3 * x1 + x3 + 2 * x4 + rnorm(100, 0, 0.1)
lm(cbind(y, X))
glm(cbind(y, X), family = "gaussian")

for (i in 1:10) {
  select(X,y)
}







# testing code
pool <- init(round(1.5*ncol(X)), ncol(X))
pool
obj1 <- getObjective(X, y, pool, "gaussian")
obj2 <- getObjective(X, y, pool, "poisson")
fitness <- getFitness(obj1)

pool[sample(1:nrow(pool), 2, prob = fitness),]
rbind(pool[sample(1:nrow(pool), 1),],
      pool[sample(1:nrow(pool), 1),])

# test rank-based selection
oneParentRandom = FALSE
childrankSelect(pool, fitness1, oneParentRandom)

pool1 <- updatePool(pool, fitness)

