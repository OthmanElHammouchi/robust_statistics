huberRho <- function(x, b) {
  if (abs(x) > b) {
    res <- b * abs(x) - b**2 / 2
    return(res)
  } else {
    res <- x ** 2 / 2
  }
}

tukeyBisquare <- function(x, c) {
  if (abs(x) > c) {
    res <- c**2 / 6
    return(res)
  } else {
    res <- x**2 / 2 - x**4 / (2 * c**2) + x**6 / (6 * c**4)
    return(res)
  }
}

library(MASS)
mu <- c(2, 3)
Sigma <- matrix(c(2, 4, 2, 5), byrow = TRUE, nrow = 2)
test <- mvrnorm(1e3, mu, Sigma)
test <- as.data.frame(test)
names(test) <- c("x", "y")

cov.est <- cov(test)
eigen.res <- eigen(cov.est)
eigen.vectors <- as.data.frame(t(eigen.res$vectors))
names(eigen.vectors) <- c("x", "y")

ggplot() +
  geom_point(aes(x, y), data = test) +
  geom_segment(
    aes(x = mu[1], y = mu[2], xend = x + mu[1], yend = y + mu[2]),
    arrow = arrow(length = unit(0.5, "cm")),
    colour = "red",
    linewidth = 2,
    data = eigen.vectors
  )
