



#test data
x1 <- rnorm(100, 2, 1)
x2 <- rnorm(100, 5, 1)
x3 <- rnorm(100, -1, 1)
X <- data.frame(x1, x2, x3)
y <- 3 * x1 + x3 + rnorm(100, 0, 0.1)
data <- cbind(y, X)

# data <- matrix(c(y, x1, x2, x3), nrow = length(y))
# chrom1 <- c(TRUE, FALSE, TRUE)
# chrom2 <- c(TRUE, TRUE, FALSE)
# data <- data.frame(
#   y,x1,x2,x3
# )

lm1 <- lm(as.data.frame(data[,c(TRUE, chrom1)]))
lm2 <- lm(as.data.frame(data[,c(TRUE, chrom2)]))
AIC(lm1)
AIC(lm2)
