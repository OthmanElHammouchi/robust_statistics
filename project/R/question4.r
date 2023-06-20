source(file.path("R", "common.r"))
data <- warpbreaks

## ---------------------------------------------------------------
##                          Poisson GLM                         -
## ---------------------------------------------------------------

## Fit Poisson GLM
model <- glm(breaks ~ wool + tension, poisson, data)
silent <- capture.output(stargazer(
  model,
  out = file.path(tab.dir, "q4-glm-res.tex"),
  single.row = TRUE,
  header = FALSE,
  float = FALSE
))

## First few rows of design matrix
# Illustrates encoding of the categorical variables
silent <- capture.output(stargazer(
  head(model.matrix(model)),
  out = file.path(tab.dir, "q4-design-mat.tex"),
  single.row = TRUE,
  header = FALSE,
  float = FALSE
))

## Box plot of response for different tensions
# Inverse relation between number of breaks and tension level
p <- ggplot(data) +
  geom_boxplot(aes(tension, breaks)) +
  geom_point(aes(tension, mean(breaks)), colour = "red", size = 0.5) +
  xlab("Tension") +
  ylab("Breaks")
ggsave(file.path(plot.dir, "q4-boxplot.eps"), p, width = width / 2, height = height / 3, units = "mm")

## ---------------------------------------------------------------
##                      Dispersion analysis                     -
## ---------------------------------------------------------------

## Response vs. fitted value
# Many observations outside 95% confidence band, indicating
# overdispersion
fitted <- data.frame(fitted = predict(model, type = "response"))
conf.band <- with(fitted, cbind.data.frame(upper = qpois(0.975, fitted), lower = qpois(0.025, fitted)))

plot.df <- cbind.data.frame(data, fitted, conf.band)
p <- ggplot(plot.df, aes(x = fitted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "red") +
  geom_point(aes(y = breaks), size = 1) +
  xlab("Fitted value") +
  ylab("Breaks")
ggsave(file.path(plot.dir, "q4-resp-fit-pois.eps"), p, width = width / 2.2, height = height / 3, units = "mm")

## Overdispersion test
# Earlier suspicion is supported
disp.test <- dispersiontest(model)
pval <- signif(disp.test$p.value, 3)

## Fit quasi-Poisson model
model.quasi <- glm(breaks ~ wool + tension, quasipoisson, data)

## Response against fitted value
# Confidence bands look much better now
disp <- summary(model.quasi)$dispersion
fitted <- data.frame(fitted = predict(model.quasi, type = "response"))
conf.band <- with(fitted,
  cbind.data.frame(upper = fitted - 3 * sqrt(disp * fitted),
    lower = fitted + 3 * sqrt(disp * fitted))
)
plot.df <- cbind.data.frame(data, fitted, conf.band)
p <- ggplot(plot.df, aes(x = fitted)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "red") +
  geom_point(aes(y = breaks), size = 1) +
  xlab("Fitted value") +
  ylab("Breaks")
ggsave(file.path(plot.dir, "q4-resp-fit-quasi.eps"), p, width = width / 2.2, height = height / 3, units = "mm")

## ----------------------------------------------------------------
##                          Robust fit                           -
## ----------------------------------------------------------------

## Robust fit
# We could find no implementation of robust glm supporting quasi-
# Poisson, so we use Poisson. This means we can compare the point
# estimates for the coefficients, but inference is not possible
# because the model does not adjust for overdispersion.
model.rob <- glmrob(breaks ~ wool + tension, poisson, data)
res <- coef(summary(model.rob))[, 1, drop = FALSE]
silent <- capture.output(stargazer(
  res,
  out = file.path(tab.dir, "q4-glm-rob-res.tex"),
  single.row = TRUE,
  header = FALSE,
  median = TRUE,
  digits = 2,
  float = FALSE,
  flip = TRUE,
  report = "vc",
  omit.table.layout = "sn"
))
