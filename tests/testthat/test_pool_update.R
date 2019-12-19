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


groupNum <- floor(poolSize / 3)
test_that('Tournament Selection',
          {test <- GA::tournamentSelect(pool, fitness, numCrossoverSplit, groupNum)
          expect_equal(length(test),chromoSize*poolSize )
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_type(test, "logical")

          })

test_that('Rank Selection',
          {test <- GA::rankSelect(pool, fitness, numCrossoverSplit, groupNum)
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "logical")
          expect_equal(length(test),chromoSize*poolSize )
          expect_type(test, "logical")

          })


test_that('Pool Update works for rank selection',
          {test <- GA::updatePool(pool, fitness,
                                  tournamentSelection, groupNum,
                                  oneParentRandom=FALSE, numCrossoverSplit,
                                  mutationRate, maxMutationRate, i, maxIter)
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "matrix")
          expect_equal(length(test),chromoSize*poolSize )
          expect_type(test, "logical")

          })

test_that('Pool Update works for rank selection with one parent random',
          {test <- GA::updatePool(pool, fitness,
                                  tournamentSelection, groupNum,
                                  oneParentRandom=TRUE, numCrossoverSplit,
                                  mutationRate, maxMutationRate, i, maxIter)
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "matrix")
          expect_equal(length(test),chromoSize*poolSize )
          expect_type(test, "logical")

          })

groupNum <- floor(poolSize / 3)

test_that('Pool Update works for tournament selection',
          {test <- GA::updatePool(pool, fitness,
                                  tournamentSelection=TRUE, groupNum,
                                  oneParentRandom=FALSE, numCrossoverSplit,
                                  mutationRate, maxMutationRate, i, maxIter)
          expect_true(any(test == 1 | test == 0))
          expect_is(test, "matrix")
          expect_equal(length(test),chromoSize*poolSize )
          expect_type(test, "logical")

          })


