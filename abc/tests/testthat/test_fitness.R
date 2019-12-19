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

#context('fitness_calculation')
#context('helper')

allData <- cbind(y, X)
chromoSize <- ncol(X)
poolSize=round(2 * ncol(X))
maxIter=100
pool <- GA::init(poolSize, chromoSize)
objValEachIter <- data.frame(
  iter = rep(NA, poolSize * maxIter),
  objectiveValue = rep(NA, poolSize * maxIter)
)


test_that('Objective function and Fitness score work',
          {test <- GA::getObjective(
                                X, y, pool, objectiveFunction=stats::AIC, regressionType="gaussian", nCores=1
                                    )
          expect_is(test, "numeric")
          expect_type(test, "double")
          expect_is(GA::getFitness(test), "numeric")
          expect_type(GA::getFitness(test), "double")

          })


test_that('Other objective functions work',
          {test <- GA::getObjective(
            X, y, pool, objectiveFunction=stats::logLik, regressionType="gaussian", nCores=1
          )
          expect_is(test, "numeric")
          expect_type(test, "double")
          })

test_that('Parallel Computation works',
          {test <- GA::getObjective(
            X, y, pool, objectiveFunction=stats::AIC, regressionType="gaussian", nCores=16
          )
          expect_is(test, "numeric")
          expect_type(test, "double")
          })
