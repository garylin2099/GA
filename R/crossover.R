#' Single Crossover
#'
#' \code{crossover} accomplishes the process that parents crossover with a single split point in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#' @param parents a 2 by C matrix of booleans, with each row representing a parent chromosome
#'
#' @export
crossover <- function(parents) {
  chromoSize <- ncol(parents)
  splitPoint <- sample(2:chromoSize, 1)
  firstChild <- c(parents[1,1:(splitPoint-1)], parents[2,splitPoint:chromoSize])
  secondChild <- c(parents[2,1:(splitPoint-1)], parents[1,splitPoint:chromoSize])
  return(c(firstChild, secondChild))
}

#' Multiple Crossover
#'
#' \code{multipleCrossover} accomplishes the process that parents crossover with multiple split points in their chromosomes. This function returns a vector of 2*C representing the genes of two children.
#'
#' @param parents a 2 by C matrix of booleans, with each row representing a parent chromosome
#' @param numSplit number of split points in a crossover
#'
#' @export
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
