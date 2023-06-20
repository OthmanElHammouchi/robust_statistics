source(file.path("R", "common.r"))

## ----------------------------------------------------------------
##                  M-estimator implementation                   -
## ----------------------------------------------------------------

## Huber's rho for b = 1.5
huberPsi <- function(x) {
  b <- 1.5
  if (abs(x) > b) {
    return(b * sign(x))
  } else {
    return(x)
  }
}

## Tukey's psi for c = 4.68
tukeyPsi <- function(x) {
  c <- 4.68
  if (abs(x) > c) {
    return(0)
  } else {
    return(x * (1 - x**2 / c**2)**2)
  }
}

weight <- function(x, type) {
  if (x == 0) {
    return(1)
  } else {
    if (type == "huber") {
      return(huberPsi(x) / x)
    } else {
      return(tukeyPsi(x) / x)
    }
  }
}

## M-estimator for location
locMEst <- function(sample, type, tol = 1e-8) {
  sigma <- mad(sample) # scale estimate
  mu <- median(sample) # initial location estimate

  k <- 1
  diff <- Inf
  while (diff > tol * sigma) {
    mu.old <- mu # keep current estimate
    weights <- sapply(((sample - mu) / sigma), weight, type = type) # compute weights
    mu <- sum(weights * sample) / sum(weights) # update estimate
    diff <- abs(mu - mu.old) # compute change
    k <- k + 1
  }
  return(mu)
}

## Sensitivity curve
sensitivity <- function(x, type, sample) {
  sample.aug <- c(sample, x)
  n <- length(sample.aug)
  est <- locMEst(sample, type = type)
  est.aug <- locMEst(sample.aug, type = type)
  res <- (est.aug - est) / (1 / n)
  return(res)
}

## Theoretical influence functions
# Normalisation constant for Tukey computed using Python script tukey.py
b <- 1.5
huberIF <- sapply(z, huberPsi) / (pnorm(b) - pnorm(-b))

c <- 4.68
tukeyIF <- sapply(z, tukeyPsi) / 0.7573

## ----------------------------------------------------------------
##                Sensitivity curves for n = 100                 -
## ----------------------------------------------------------------

## Simulate symmetric standard normal sample
sample <- rnorm(100 / 2)
sample <- c(sample, -sample)

## Calculate sensitivity curves
z <- seq(-4, 4, by = 1e-2)
tukey <- sapply(z, sensitivity, type = "tukey", sample = sample)
huber <- sapply(z, sensitivity, type = "huber", sample = sample)

## Plot sensitivity curves
huber.df <- data.frame(z = z, huber = huber, huberIF = huberIF)
p <- ggplot(huber.df) +
  geom_line(aes(z, huber), colour = "red") +
  geom_line(aes(z, huberIF), colour = "black") +
  xlab("z") +
  ylab("Sensitivity") +
  theme(legend.box.spacing = unit(0, "mm"))
ggsave(file.path(plot.dir, "q1-sens_curves_huber_small.eps"),
  p,
  width = width / 2,
  height = height / 2.5,
  units = "mm"
)

tukey.df <- data.frame(z = z, tukey = tukey, tukeyIF = tukeyIF)
p <- ggplot(tukey.df) +
  geom_line(aes(z, tukey), colour = "red") +
  geom_line(aes(z, tukeyIF), colour = "black") +
  xlab("z") +
  ylab("Sensitivity") +
  theme(legend.box.spacing = unit(0, "mm"))
ggsave(file.path(plot.dir, "q1-sens_curves_tukey_small.eps"),
  p,
  width = width / 2,
  height = height / 2.5,
  units = "mm"
)

## ----------------------------------------------------------------
##                Sensitivity curves for n = 1e6                 -
## ----------------------------------------------------------------
## Simulate symmetric standard normal sample
sample <- rnorm(1e6 / 2)
sample <- c(sample, -sample)

## Compute sensitivity curves
# Warning: this code takes ~7 min. to run on an Intel i7-12700H
# processor with 20 virtual cores. For convenience, the result
# from the run has been saved and can be loaded by uncommenting
# and running the following line:
#
# res <- readRDS(file.path("results", "data", "q1_sens_data.RDS"))

z <- seq(-4, 4, by = 1e-2)
library(parallel)
library(doSNOW)
library(foreach)
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

pgb <- txtProgressBar(max = length(z))
opts <- list(progress = function(i) setTxtProgressBar(pgb, i))

res <- foreach(val = z, .combine = rbind, .options.snow = opts) %dopar% {
  tukey <- sensitivity(val, type = "tukey", sample = sample)
  huber <- sensitivity(val, type = "huber", sample = sample)
  c(tukey, huber)
}
stopCluster(cl)
res <- as.data.frame(res)
colnames(res) <- c("tukey", "huber")
rownames(res) <- NULL

## Plot sensitivity curves
huber.df <- data.frame(z = z, huber = res$huber, huberIF = huberIF)
p <- ggplot(huber.df) +
  geom_line(aes(z, huber), colour = "red") +
  geom_line(aes(z, huberIF), colour = "black") +
  xlab("z") +
  ylab("Sensitivity") +
  theme(legend.box.spacing = unit(0, "mm"))
ggsave(file.path(plot.dir, "q1-sens_curves_huber_large.eps"),
  p,
  width = width / 2,
  height = height / 2.5,
  units = "mm"
)

tukey.df <- data.frame(z = z, tukey = res$tukey, tukeyIF = tukeyIF)
p <- ggplot(tukey.df) +
  geom_line(aes(z, tukey), colour = "red") +
  geom_line(aes(z, tukeyIF), colour = "black") +
  xlab("z") +
  ylab("Sensitivity") +
  theme(legend.box.spacing = unit(0, "mm"))
ggsave(file.path(plot.dir, "q1-sens_curves_tukey_large.eps"),
  p,
  width = width / 2,
  height = height / 2.5,
  units = "mm"
)

# Sensitivity curves approach influence functions as n increases
