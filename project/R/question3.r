source(file.path("R", "common.r"))
load(file.path("data", "Topgear.rdata"))

## ---------------------------------------------------------------
##                              EDA                             -
## ---------------------------------------------------------------

## Data cleaning
rownames(data) <- NULL # reset row indices
data <- as.data.frame(lapply(data, function(col) {
  if (is.factor(col)) col <- droplevels(col)
  col
})) # drop unused factor levels

## Summary statistics for the response
# Discrepancy mean and median => possible outliers.
silent <- capture.output(stargazer(
  data[, "Price", drop = FALSE],
  out = file.path(tab.dir, "q3-summary.tex"),
  single.row = TRUE,
  header = FALSE,
  median = TRUE,
  digits = 2,
  float = FALSE
))

## Correlation plot
# Two identifiable feature groups through hclust.
# Interdependence especially high between physical power
# measurements => possible multicollinearity.
colfunc <- colorRampPalette(c("blue", "white", "red"))
eps(file = file.path(plot.dir, "q3-corplot.eps"), width = 0.04 * width, height = 0.04 * height / 2, pointsize = 8)
corrplot(cor(select(data, -c(Maker, Model))), order = "hclust", addrect = 2, col = colfunc(105), tl.col = "black")
dev.off()

## Cluster car makers
# Difficult to use 'Maker' as-is due to high number of levels
# => solution: cluster on mean price to obtain categories.
mean.price <- tapply(data$Price, data$Maker, mean)
diss <- dist(as.matrix(mean.price))
res <- fastkmed(diss, 3)
price.cat <- as.factor(res$cluster)
levels(price.cat) <- c("High-end", "Mid-range", "Entry-level")
price.cat <- cbind.data.frame(Maker = names(price.cat), Cat = price.cat)
data <- left_join(price.cat, data)

## Scatter plots of price against different predictors
# Clustering used to colour price categories.
# Results do not support the linear model: clear difference
# in slope between different groups, heteroscedasticity, extreme observations.
# The latter correspond to hypercar manufacturers Bugatti and Pagani
# and are therefore understandable.
vars <- c("Displacement", "BHP", "Torque", "Acceleration", "MPG")
hyper.cars <- c("Bugatti", "Pagani")
plt.list <- list()
for (i in seq_along(vars)) {
  p <- ggplot(data, aes(.data[[vars[i]]], Price)) +
    geom_point(aes(colour = Cat), size = 0.5) +
    geom_text_repel(
      aes(label = ifelse(Maker %in% hyper.cars, as.character(Maker), "")),
      size = (5 / 14) * 6
    ) +
    scale_colour_hue("Price category") +
    theme(
      legend.direction = "vertical",
      axis.title = element_text(size = 6)
    ) +
    xlab(vars[i]) +
    ylab("Price")

  plt.list <- c(plt.list, list(p))
}
p <- wrap_plots(plt.list, guides = "collect", nrow = 2) + guide_area()
ggsave(
  file.path(plot.dir, "q3-scatter-classic.eps"),
  p,
  height = width / 2,
  width = height,
  units = "mm"
)

## ---------------------------------------------------------------
##                  Classical linear regression                 -
## ---------------------------------------------------------------

## -----------------
##  Original data
## -----------------

## Fit linear model to the original data
model <- lm(Price ~ Displacement + BHP + Torque + Acceleration + MPG, data)

## QQ-plot
p <- ggplot(data.frame(resid = rstandard(model))) +
  geom_qq(aes(sample = resid), size = 0.5) +
  geom_qq_line(aes(sample = resid)) +
  ylab("Standardised residuals") +
  xlab("Theoretical quantiles")
ggsave(file.path(plot.dir, "q3-qq-classic.eps"), p, width = width / 2, height = height / 3, units = "mm")

## Standardised residuals vs. fitted values
plt.df <- data.frame(maker = data$Maker, resid = rstandard(model), fitted = predict(model))
p <- ggplot(plt.df, aes(fitted, resid)) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  geom_point(size = 1) +
  geom_text_repel(aes(label = ifelse(abs(resid) > 2.5, as.character(maker), "")), size = (5 / 14) * 6) +
  xlab("Fitted values") +
  ylab("Standardised residuals")
ggsave(file.path(plot.dir, "q3-resids-vs-fitted-classic.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

# Suspicions based on scatterplots confirmed: diagnostic plots
# clearly show that linear model is not suitable.
# Downward pattern in resids vs. fitted, fat tails in QQ-plot.

## ------------------------
##  Transformed response
## ------------------------

## Scatter plots with log-transformed response
# Linear model seems to fit much better now: variance is stabilised,
# previous outliers become part of line pattern. However: bad leverage
# outliers for Acceleration and MPG.
plt.list <- list()
for (i in seq_along(vars)) {
  p <- ggplot(data, aes(.data[[vars[i]]], log(Price))) +
    geom_point(aes(colour = Cat), size = 0.5) +
    geom_text_repel(
      aes(label = ifelse(Maker %in% hyper.cars, as.character(Maker), "")),
      size = (5 / 14) * 6
    ) +
    scale_colour_hue("Price category") +
    theme(
      legend.direction = "vertical",
      axis.title = element_text(size = 6)
    ) +
    xlab(vars[i]) +
    ylab("Price")

  if (vars[i] == "MPG") {
    p <- p + geom_text_repel(
      aes(label = ifelse(MPG > 100, paste(as.character(Maker), as.character(Model), sep = " "), "")),
      size = (5 / 14) * 6
    )
  }

  if (vars[i] == "Acceleration") {
    p <- p + geom_text_repel(
      aes(label = ifelse(Acceleration == 0, paste(as.character(Maker), as.character(Model), sep = " "), "")),
      size = (5 / 14) * 6
    )
  }

  plt.list <- c(plt.list, list(p))
}
p <- wrap_plots(plt.list, guides = "collect", nrow = 2) + guide_area()
ggsave(file.path(plot.dir, "q3-scatter-tf.eps"),
  p,
  height = width / 2,
  width = height,
  units = "mm"
)

## Fit linear model with transformed response
model.tf <- lm(log(Price) ~ Displacement + BHP + Torque + Acceleration + MPG, data)

## QQ-plot
p <- ggplot(data.frame(resid = rstandard(model.tf))) +
  geom_qq(aes(sample = resid), size = 0.5) +
  geom_qq_line(aes(sample = resid)) +
  ylab("Standardised residuals") +
  xlab("Theoretical quantiles")
ggsave(file.path(plot.dir, "q3-qq-tf.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

## Standardised residuals vs. fitted values
plt.df <- data.frame(maker = data$Maker, resid = rstandard(model.tf), fitted = predict(model.tf))
p <- ggplot(plt.df, aes(fitted, resid)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = ifelse(abs(resid) > 2.5, as.character(maker), "")), size = (5 / 14) * 6) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  xlab("Fitted values") +
  ylab("Standardised residuals")
ggsave(file.path(plot.dir, "q3-resids-vs-fitted-tf.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

# Residual plots show considerable improvement: tails are less pronounced,
# (though still more than normality), pattern is gone from resids vs. fitted.

## Investigate bad leverage outliers
# MPG outliers are electric-powered vehicles => MPG not meaningful.
# Acceleration outlier is clearly error (0 seconds to 62 mph).
data[data$MPG > 100, ]
data[data$Acceleration < 2, ]

## Neither classical nor robust fit flag outliers correctly
data[which(rstandard(model.tf) > 2.5), c("Maker", "Model")]
model.rob <- ltsReg(log(Price) ~ Displacement + BHP + Torque + Acceleration + MPG, data)
data[which(model.rob$residuals / model.rob$scale > 2.5), c("Maker", "Model")]

## ----------------------------------
##  Corrected and transformed data
## ----------------------------------

## Exclude problem points
data <- data[data$MPG <= 100, ]
data <- data[data$Acceleration != 0, ]

## Refit the linear model
model.tf <- lm(log(Price) ~ Displacement + BHP + Torque + Acceleration + MPG, data)
res <- coef(summary(model.tf))[, 1, drop = FALSE] # don't include p-values (fat tails)
silent <- capture.output(stargazer(
  res,
  out = file.path(tab.dir, "q3-reg-res-tf.tex"),
  single.row = TRUE,
  header = FALSE,
  median = TRUE,
  digits = 2,
  float = FALSE,
  flip = TRUE,
  report = "vc",
  omit.table.layout = "sn"
))

## QQ-plot
p <- ggplot(data.frame(resid = rstandard(model.tf))) +
  geom_qq(aes(sample = resid), size = 0.5) +
  geom_qq_line(aes(sample = resid)) +
  ylab("Standardised residuals") +
  xlab("Theoretical quantiles")
ggsave(file.path(plot.dir, "q3-qq-corrected.eps"),
  p,
  width = width / 2,
  height = height / 3.5,
  units = "mm"
)

## Standardised residuals vs. fitted values
plt.df <- data.frame(maker = data$Maker, resid = rstandard(model.tf), fitted = predict(model.tf))
p <- ggplot(plt.df, aes(fitted, resid)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = ifelse(abs(resid) > 2.5, as.character(maker), "")), size = (5 / 14) * 6) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  xlab("Fitted values") +
  ylab("Standardised residuals")
ggsave(file.path(plot.dir, "q3-resids-vs-fitted-corrected.eps"),
  p,
  width = width / 2,
  height = height / 3.5,
  units = "mm"
)

# Diagnostic plots broadly similar to previous model. Heavy tails in the
# distribution of the residuals persist, which is worrying. It suggests
# either that we should replace the Gaussian distribution with a more
# appropriate one, or that the model fails to capture some pattern.

## Classical outlier map
X <- model.matrix(model.tf)
X <- X[, 2:ncol(X)]
nvars <- ncol(X)
nobs <- nrow(X)

plt.df <- data.frame(maker = data$Maker,
  model = data$Model,
  resid = rstandard(model.tf),
  md = sqrt(mahalanobis(X, colMeans(X), cov(X)))
)
p <- ggplot(plt.df, aes(md, resid)) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  geom_vline(xintercept = sqrt(qchisq(0.975, nvars)), colour = "red") +
  geom_point(size = 0.5) +
  xlab("Mahalanobis distance") +
  ylab("Standardized LS residual") +
  geom_text_repel(
    aes(label = ifelse(abs(resid) > 2.5 | md > sqrt(qchisq(0.975, nvars)),
      paste(as.character(maker), as.character(model), sep = " "),
      "")
    ),
    size = (5 / 14) * 6
  )
ggsave(file.path(plot.dir, "q3-outlier-classic.eps"),
  p,
  width = height / 2,
  height = width / 1.5,
  units = "mm"
)

## ----------------------------------------------------------------
##                          Robust fit                           -
## ----------------------------------------------------------------

## Fit linear model using the robust LTS estimator
model.rob <- ltsReg(log(Price) ~ Displacement + BHP + Torque + Acceleration + MPG, data)
res <- coef(summary(model.rob))[, 1, drop = FALSE]
silent <- capture.output(stargazer(
  res,
  out = file.path(tab.dir, "q3-reg-res-rob.tex"),
  single.row = TRUE,
  header = FALSE,
  median = TRUE,
  digits = 2,
  float = FALSE,
  flip = TRUE,
  report = "vc",
  omit.table.layout = "sn"
))

## QQ-plot
p <- ggplot(data.frame(resid = model.rob$residuals / model.rob$scale)) +
  geom_qq(aes(sample = resid), size = 0.5) +
  geom_qq_line(aes(sample = resid)) +
  ylab("LTS residuals") +
  xlab("Theoretical quantiles")
ggsave(file.path(plot.dir, "q3-qq-rob.eps"),
  p,
  width = width / 2,
  height = height / 3.5,
  units = "mm"
)

## LTS residuals vs. fitted values
plt.df <- data.frame(
  maker = data$Maker,
  resid = model.rob$residuals / model.rob$scale,
  fitted = model.rob$fitted.values
)
p <- ggplot(plt.df, aes(fitted, resid)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = ifelse(abs(resid) > 2.5, as.character(maker), "")), size = (5 / 14) * 6) +
  geom_hline(yintercept = c(-2, 2), colour = "red") +
  xlab("Fitted values") +
  ylab("LTS residuals")
ggsave(file.path(plot.dir, "q3-resids-vs-fitted-rob.eps"),
  p,
  width = width / 2,
  height = height / 3.5,
  units = "mm"
)

## Robust outlier map
X <- model.rob$X
X <- X[, 2:ncol(X)]
nvars <- ncol(X)
nobs <- nrow(X)

plt.df <- data.frame(
  maker = data$Maker,
  model = data$Model,
  resid = model.rob$residuals / model.rob$scale,
  rd = model.rob$RD
)
p <- ggplot(plt.df, aes(rd, resid)) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  geom_vline(xintercept = sqrt(qchisq(0.975, nvars)), colour = "red") +
  geom_point(size = 0.5) +
  geom_text_repel(
    aes(label = ifelse(
      abs(resid) > 2.5 | rd > sqrt(qchisq(0.975, nvars)),
      paste(as.character(maker), as.character(model), sep = " "),
      "")
    ),
    size = (5 / 14) * 6
  ) +
  xlab("Robust distance computed by MCD") +
  ylab("Standardized LTS residual")
ggsave(file.path(plot.dir, "q3-outlier-rob.eps"),
  p,
  width = height / 2,
  height = width / 1.5,
  units = "mm"
)

# Robust outlier map shows much higher number of good leverage outliers
# than classical one. Again, this indicates either a fat-tailed response
# distribution or the presence of patterns not captured by the model.
# Based on the scatterplots, we try adding a group effect to see if
# this might help.

## ----------------------------------------------------------------
##                  Incorporating price category                 -
## ----------------------------------------------------------------

## Fit linear model including price category effect
model.tf <- lm(log(Price) ~ Displacement + BHP + Torque + Acceleration + MPG + Cat, data)
res <- coef(summary(model.tf))[, c(1, 4)]
colnames(res) <- c("Estimate", "p-value")
res.file <- file(file.path(tab.dir, "q3-reg-res-price-cat.tex"))
writeLines(capture.output(print(xtable(t(res), display = rep("g", 9)), floating = FALSE, booktabs = TRUE)), res.file)
close(res.file)

## QQ-plot
p <- ggplot(data.frame(resid = rstandard(model.tf))) +
  geom_qq(aes(sample = resid), size = 0.5) +
  geom_qq_line(aes(sample = resid)) +
  ylab("Standardised residuals") +
  xlab("Theoretical quantiles")
ggsave(file.path(plot.dir, "q3-qq-price-cat.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

## Standardised residuals vs. fitted values
plt.df <- data.frame(maker = data$Maker, resid = rstandard(model.tf), fitted = predict(model.tf))
p <- ggplot(plt.df, aes(fitted, resid)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = ifelse(abs(resid) > 2.5, as.character(maker), "")), size = (5 / 14) * 6) +
  geom_hline(yintercept = c(-2.5, 2.5), colour = "red") +
  xlab("Fitted values") +
  ylab("Standardised residuals")
ggsave(file.path(plot.dir, "q3-resids-vs-fitted-price-cat.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

# Diagnostics are considerably improved, to the point that
# the remaining tail observations appear to be insufficient to
# invalidate inference based on normality.

## Normality tests
# Both fail to reject normality, though Shapiro-Wilk only narrowly
shapiro.test(rstandard(model.tf))
ad.test(rstandard(model.tf))
