y <- mtcars$mpg
x <- as.matrix(mtcars[, c(-1)])

x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
x4 <- rnorm(100, 18, 1)
x5 <- rnorm(100, -7, 1)
x6 <- rnorm(100, -20, 1)
x7 <- rnorm(100, 15, 1)
x8 <- rnorm(100, -9, 1)
X1 <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
y1 <- x1 + x3 + x4 - 2 * x6 - x8 + rnorm(100, 0, 0.1)


test_that('test for invalid input',
          {
            expect_error(GA::select(x, y, objective_function = "AIC"))
            expect_error(GA::select(x, y, oneParentRandom="hello world"))
            expect_error(GA::select(x, y, mutationRate=0.3,maxMutationRate = 0.1))
            expect_error(GA::select(x, y, numCrossoverSplit=30))
            expect_error(GA::select(x,y, "foo", regressionType = "gaussian"))
            expect_error(GA::select(x, regressionType = "gaussian"))
            expect_error(GA::select(x,y[-1]))
            expect_error(GA::select(x, y, regressionType = "chi-square"))
            expect_error(GA::select(x, y, nCores = -3))
            expect_error(GA::select(cbind(y, y), x))
            expect_error(GA::select( list(c(1,2,3,4),c(5,6,7,8),5)))

          })


test_that('GA works',
          {test <- GA::select(x, y, objValPlot = FALSE)
          expect_type(test, "list") #list
          expect_type(GA::select(x,y,objValPlot = FALSE)$model,
                      "list")

          })
test_that('GA works 2',
          {test <- select(X1, y1,
                          maxMutationRate = 0.02, numCrossoverSplit = 3,
                          tournamentSelection = FALSE, objValPlot = FALSE)
          expect_type(test, "list") #list
          expect_type(GA::select(X1,y1,objValPlot = FALSE)$model,
                      "list")
          expect_equal(GA::select(X1, y1,
                                  maxMutationRate = 0.02,
                                  numCrossoverSplit = 3,
                                  tournamentSelection = FALSE,
                                  objValPlot = FALSE)$convergence,TRUE)
          })

test_that('Diverge',
          {expect_equal(GA::select(x, y, maxIter = 100,
                                   mutationRate  = 0.8,
                                   objValPlot = FALSE)$convergence, FALSE)
            expect_equal(GA::select(x,y, maxIter = 100,
                                    mutationRate = 0.8,
                                    objValPlot = FALSE)$iterationCount, 100)
          })
