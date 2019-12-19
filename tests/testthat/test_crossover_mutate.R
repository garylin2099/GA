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

allData <- cbind(y, X)
chromoSize <- ncol(X)
poolSize=round(2 * ncol(X))
maxIter=100
pool <- GA::init(poolSize, chromoSize)
objValEachIter <- data.frame(
  iter = rep(NA, poolSize * maxIter),
  objectiveValue = rep(NA, poolSize * maxIter)
)
groupNum=NULL


objVal <- GA::getObjective(
  X, y, pool, objectiveFunction=stats::AIC, regressionType="gaussian", nCores=1
)
fitness <- GA::getFitness(objVal)
tournamentSelection = FALSE
numCrossoverSplit = 1
mutationRate = 0.01
maxMutationRate = NULL
maxIter = 100
minIter = round(maxIter / 10)
iterCounter=2



poolSize <- nrow(pool)
chromoSize <- ncol(pool)
childVector <- rep(NA, poolSize*chromoSize)
childCount <- 0
oneParentRandom=FALSE


parents <- rbind(pool[sample(1:poolSize, 1, prob = fitness),],
                     pool[sample(1:poolSize, 1),])


test_that('Crossover method 1 ',
          {test <- GA::crossover(parents)
          expect_equal(length(test),length(parents))
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_type(test, "logical")
          })

numCrossoverSplit = 3
test_that('Crossover method 2 ',
          {test <- GA::multipleCrossover(parents, numCrossoverSplit)
          expect_equal(length(test),length(parents))
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_type(test, "logical")
          })


groupNum <- floor(poolSize / 3)
childVector <- GA::tournamentSelect(pool, fitness, numCrossoverSplit, groupNum)

test_that('Mutate method 1 ',
          {test <- GA::mutate(childVector, mutationRate)
          expect_equal(length(test),chromoSize*poolSize )
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_type(test, "logical")

          })


maxMutationRate = 0.05
test_that('Mutate method 2 ',
          {test <- GA::mutate(childVector, currMutationRate(mutationRate, maxMutationRate, iterCounter, maxIter))
          expect_equal(length(test),chromoSize*poolSize )
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_type(test, "logical")
          })

