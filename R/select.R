rm(list = ls())

select <- function(X, y, poolSize = round(1.5 * ncol(X)), mutateRate = 0.01, maxIter = 100, genGap = 1) {
  allData <- cbind(y, X)
  chromoSize <- ncol(X)
  pool <- init(poolSize, chromoSize)
  for (i in 1:maxIter) {
    objVal <- getObjective(allData, pool)
    fitness <- getFitness(objVal)
    pool <- updatePool(pool, fitness)
  }
  selectedIndex <- which(colSums(pool)==max(colSums(pool)))
  return(cat("Variable selected are", colnames(X)[selectedIndex]))

}
init <- function(poolSize, chromSize) {
  return(matrix(sample(c(TRUE, FALSE), chromSize * poolSize, replace = TRUE), nrow = poolSize, byrow = TRUE))
}


getObjective <- function(data, pool) {
  fitScore <- rep(NA, nrow(pool))
  getAIC <- function(chromo) {
    return(AIC(lm(data[c(TRUE, chromo)])))
  }
  return(apply(pool, 1, getAIC))
}

getFitness <- function(objectiveVal) {
  poolSize <- length(objectiveVal)
  r <- poolSize + 1 - rank(objectiveVal) # let small AIC have higher ranking value
  return(2*r/(poolSize*(poolSize+1)))
}

updatePool <- function(pool, fitness) {
  poolSize <- nrow(pool)
  chromoSize <- ncol(pool)
  childVector <- NULL
  childCount <- 0
  while(childCount < poolSize) {
    parents <- pool[sample(1:poolSize, 2, prob = fitness),]
    childs <- crossover(parents)
    childVector <- c(childVector, childs)
    childCount = childCount + 2
  }
  childVector <- mutate(childVector)
  return(matrix(childVector[1:(poolSize*chromoSize)], nrow = poolSize, byrow = TRUE))
}

crossover <- function(parents) {
  chromoSize <- ncol(parents)
  splitPoint <- sample(2:chromoSize, 1)
  firstChild <- c(parents[1,1:(splitPoint-1)], parents[2,splitPoint:chromoSize])
  secondChild <- c(parents[2,1:(splitPoint-1)], parents[1,splitPoint:chromoSize])
  return(c(firstChild, secondChild))
}

mutate <- function(childVector) {
  mutateGeneIndex <- runif(length(childVector)) <= 0.01
  childVector[mutateGeneIndex] <- !childVector[mutateGeneIndex]
  return(childVector)
}


#test data
x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
x4 <- rnorm(100, 18, 1)
x5 <- rnorm(100, -7, 1)
X <- data.frame(x1, x2, x3, x4, x5)
y <- 3 * x1 + x3 + 2 * x4 + rnorm(100, 0, 0.1)
lm(y ~ ., data = X)

select(X,y)









# # testing code
# poolSize <- round(1.5*ncol(X))
# chromoSize <- ncol(X)
# pool <- init(round(1.5*ncol(X)), ncol(X))
# pool
# obj <- getObjective(data, pool)
# fitness <- getFitness(obj)
# pool1 <- updatePool(pool, fitness)

