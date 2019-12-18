a="asd"
expect_output(str(a),"asd")

rm(list = ls())


init <- function(poolSize, chromSize) {
  return(matrix(sample(c(TRUE, FALSE), chromSize * poolSize, replace = TRUE), nrow = poolSize, byrow = TRUE))
}


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

crossover <- function(parents) {
  chromoSize <- ncol(parents)
  splitPoint <- sample(2:chromoSize, 1)
  firstChild <- c(parents[1,1:(splitPoint-1)], parents[2,splitPoint:chromoSize])
  secondChild <- c(parents[2,1:(splitPoint-1)], parents[1,splitPoint:chromoSize])
  return(c(firstChild, secondChild))
}

multipleCrossover <- function(parents, numSplit) {
  chromoSize <- ncol(parents)
  # a splitIndex of value i means the split point is before the ith gene
  splitIndex <- sample(2:chromoSize, numSplit)
  # manually append an index to facilitate assigning values to the last section of the child chromosome
  splitIndex <- c(sort(splitIndex), chromoSize + 1)
  firstChild <- rep(NA, chromoSize)
  secondChild <- rep(NA, chromoSize)
  previousSplitIndex <- 1
  for (i in 1:(numSplit + 1)) {
    firstChild[previousSplitIndex:(splitIndex[i]-1)] <- parents[2 - i %% 2, previousSplitIndex:(splitIndex[i]-1)]
    secondChild[previousSplitIndex:(splitIndex[i]-1)] <- parents[1 + i %% 2, previousSplitIndex:(splitIndex[i]-1)]
    previousSplitIndex <- splitIndex[i]
  }
  return(c(firstChild, secondChild))
}

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

tournamentSelect <- function(pool, fitness, numCrossoverSplit, groupNum) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  winners <- matrix(rep(NA, poolSize * chromoSize), nrow = poolSize)
  winnerCounter <- 1
  while (winnerCounter <= poolSize) {
    groupSize <- floor(poolSize / groupNum)
    groupAssignment <- sample(1:poolSize, groupNum * groupSize)
    for (i in 1:groupNum) {
      indexInPool <- groupAssignment[((i - 1) * groupSize + 1):(i * groupSize)]
      maxIndexInPool <- indexInPool[which.max(fitness[indexInPool])]
      winners[winnerCounter,] <- pool[maxIndexInPool,]
      winnerCounter = winnerCounter + 1
      if (winnerCounter > poolSize) {
        break
      }
    }
  }
  childVector <- rep(NA, poolSize*chromoSize)
  childCount <- 0
  pairAssignment <- sample(1:poolSize, poolSize)
  while(childCount < poolSize) {
    parents <- winners[pairAssignment[(childCount+1):(childCount+2)],]
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

currMutationRate <- function(mutationRate, maxMutationRate, iterCounter, maxIter) {
  return(mutationRate + (iterCounter - 1) / maxIter * (maxMutationRate - mutationRate))
}

mutate <- function(childVector, mutationRate) {
  mutateGeneIndex <- runif(length(childVector)) <= mutationRate
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}

updatePool <- function(pool, fitness, tournamentSelection, groupNum, oneParentRandom, numCrossoverSplit,
                       mutationRate, maxMutationRate, iterCounter, maxIter) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  if (!tournamentSelection) {
    childVector <- rankSelect(pool, fitness, oneParentRandom, numCrossoverSplit)
  } else {
    childVector <- tournamentSelect(pool, fitness, numCrossoverSplit, groupNum)
  }

  if (is.null(maxMutationRate)) {
    childVector <- mutate(childVector, mutationRate)
  } else {
    childVector <- mutate(childVector, currMutationRate(mutationRate, maxMutationRate, iterCounter, maxIter))
  }

  return(matrix(childVector[1:(poolSize*chromoSize)], nrow = poolSize, byrow = TRUE))
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

select <- function(X,
                   y,
                   poolSize = round(2 * ncol(X)),
                   objectiveFunction = stats::AIC,
                   tournamentSelection = FALSE,
                   groupNum = floor(poolSize / 3), # gurantees that each group has at least three chromosomes
                   numCrossoverSplit = 1, # need to be smaller than chromosome size
                   mutationRate = 0.01,
                   maxMutationRate = NULL,
                   maxIter = 100,
                   minIter = maxIter / 10,
                   regressionType = "gaussian",
                   oneParentRandom = FALSE,
                   diversityCutoff = max(0.05, 1 / poolSize),
                   nCores = 1) {
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
       main = "Objective Function Values (e.g AIC) of Each Generation",
       xlab = "Generation i", ylab = "Objective Function Values (e.g. AIC)",
       ylim = 1.5 * c(min(objValue), max(objValue)))
  majorChromo <- getMajorChromo(pool)
  resultModel <- glm(cbind(y, X[majorChromo]), family = regressionType)
  result <- list(
    "selectedVariablesIndex" = which(majorChromo),
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

#test data
x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
x4 <- rnorm(100, 18, 1)
x5 <- rnorm(100, -7, 1)
x6 <- rnorm(100, -20, 1)
x7 <- rnorm(100, 15, 1)
x8 <- rnorm(100, -9, 1)
X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
y <- x1 + x3 + x4 - 2 * x6 - x8 + rnorm(100, 0, 0.1)
lm(cbind(y, X))
glm(cbind(y, X), family = "gaussian")

for (i in 1:1) {
  result <- select(X, y, maxMutationRate = 0.05, numCrossoverSplit = 3, tournamentSelection = FALSE)
}
expect_equal(result$selectedVariablesIndex, c(1,3,4,6,8))
