source(file.path("R", "common.r"))
load(file.path("data", "foodData.rdata"))
data <- unique(foodData)

## Data cleaning and augmentation
labeled <- data
labeled$food <- rownames(data)
labeled <- left_join(labeled, food[, c("grp_id", "food")])
labeled <- unique(dplyr::select(left_join(labeled, food_group), -c(grp_id, food)))

##---------------------------------------------------------------
##                    Variable transformation                   -
##---------------------------------------------------------------

## QQ-plot of Copper_mcg variable
plt.df <- cbind.data.frame(
  name = rownames(foodData),
  setNames(qqnorm(foodData$Copper_mcg, plot.it = FALSE), c("theoretical", "sample"))
)
p <- ggplot(plt.df, aes(theoretical, sample)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = ifelse(sample > 6, name, "")), size = (5 / 14) * 6) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")
ggsave(file.path(plot.dir, "q2-qq-orig.eps"),
  p,
  width = width,
  height = height / 2,
  units = "mm"
)

## Robust power transforms of Copper_mcg
res.yj <- transfo(temp, type = "YJ")
res.yj <- transfo(data$Copper_mcg[data$Copper_mcg > 1e-3], type = "YJ")
res.bc <- transfo(data$Copper_mcg, type = "BC")
data.tf <- data.frame(yj = as.numeric(res.yj$Xt), bc = as.numeric(res.bc$Xt))

## QQ-plots of transformed variables
# Both transformations result in an ugly deformity in the lower tail
p <- ggplot(data.tf, aes(sample = yj)) +
  geom_qq(size = 0.5) +
  geom_qq_line() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

ggsave(file.path(plot.dir, "q2-qq-tf-yj.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

p <- ggplot(data.tf, aes(sample = bc)) +
  geom_qq(size = 0.5) +
  geom_qq_line() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

ggsave(file.path(plot.dir, "q2-qq-tf-bc.eps"),
  p,
  width = width / 2,
  height = height / 3,
  units = "mm"
)

## Sanity check
# No value of lambda on a reasonable grid can get rid of the lower deformity
lambda.grid <- seq(-4, 4, by = 1e-2)
nlams <- length(lambda.grid)
nobs <- length(foodData$Copper_mcg)
res.tf <- matrix(rep(0, 2 * nlams * nobs), ncol = 2)

pgb <- txtProgressBar(max = nlams)
for (i in seq_along(lambda.grid)) {
  setTxtProgressBar(pgb, i)
  lambda <- lambda.grid[i]
  res.tf[(1 + (i - 1) * nobs):(i * nobs), 1] <- BoxCox(foodData$Copper_mcg, lambda)
  res.tf[(1 + (i - 1) * nobs):(i * nobs), 2] <- lambda
}
res.tf <- setNames(as.data.frame(res.tf), c("tf", "lambda"))

ggplot(res.tf, aes(sample = tf)) +
  geom_qq() +
  geom_qq_line() +
  transition_time(lambda) +
  view_follow() +
  labs(title = "Lambda: {frame_time}")

##----------------------------------------------------------------
##                    Dimensionality reduction                   -
##----------------------------------------------------------------

## Apply robust Yeo-Johnson transform to all variables
res <- transfo(data, type = "YJ")
data.tf <- as.data.frame(res$Xt)

## Fit ROBPCA
model <- PcaHubert(data.tf, k = 15, kmax = ncol(data.tf), mcd = FALSE)

## Make screeplot and determine k
# 6 PCs seem sufficient
plt.df <- data.frame(idx = seq_len(length(model$eigenvalues)), eig = model$eigenvalues)
p <- ggplot(plt.df, aes(idx, eig)) +
  geom_col(colour = "black", fill = "white") +
  geom_line() +
  geom_point() +
  xlab("Number of components") +
  ylab("Variance")
ggsave(file.path(plot.dir, "q2-screeplot.eps"),
  p,
  width = height / 2,
  height = width / 2,
  units = "mm"
)

## Refit ROBPCA with chosen k
model <- PcaHubert(data.tf, k = 6, mcd = FALSE)

## Loadings matrix plot
# We can tentatively observe some patterns in the loadings matrix.
# There seem to be 2 distinct large variable groups onto which the 2 PCs load
# heavily, with the remainder only associated to a few.
# However, expert judgement is required to assign meanings to these.
sorted <- faSort(model$loadings, reflect = TRUE)$loadings # `reflect` => largest loadings positive (interpretability)
colnames(sorted) <- paste0("PC", seq_len(ncol(sorted)))
loadings <- melt(sorted)
colnames(loadings) <- c("orig", "pc", "loading")

cutoff <- sqrt(1 / ncol(data.tf))
colfunc <- colorRampPalette(c("blue", "white", "red"))
p <- ggplot(loadings, aes(pc, orig)) +
  geom_raster(aes(fill = loading)) +
  geom_text(aes(label = round(loading, 3), fontface = ifelse(abs(loading) > cutoff, 2, 1)), size = (5 / 14) * 8) +
  scale_fill_gradientn(colours = colfunc(105), limits = c(-1, 1)) +
  guides(fill = guide_colourbar(title = "Loading", title.vjust = 0.8)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 8, face = "bold")
  ) +
  scale_y_discrete(limits = rev)
ggsave(file.path(plot.dir, "q2-loadings.eps"),
  p,
  width = height / 2,
  height = width,
  units = "mm"
)

## Outlier plot
# Good and bad leverage points are dominated by carbohydrate-rich
# foods (baked products and breakfast cereals). This might aid in
# interpretation, but again expert judgement seems required.
plt.df <- data.frame(sd = model$sd, od = model$od, group = labeled$group)

incl.legend <- unique(as.character(sapply(list(
  plt.df[plt.df$sd > model$cutoff.sd & plt.df$od < model$cutoff.od, ], # good leverage
  plt.df[plt.df$od > model$cutoff.od & plt.df$sd < model$cutoff.sd, ], # orthogonal
  plt.df[plt.df$od > model$cutoff.od & plt.df$sd > model$cutoff.sd, ] # bad leverage
), function(df) {
  tmp <- arrange(count(group_by(as_tibble(df), group)), desc(n))
  tmp[!is.na(tmp$group), ][1:3, ]$group
}
)))

p <- ggplot(plt.df) +
  geom_point(aes(sd, od, colour = factor(group)), size = 0.5) +
  gghighlight(
    sd > model$cutoff.sd | od > model$cutoff.od,
    unhighlighted_params = list(colour = "black"),
    keep_scales = TRUE
  ) +
  geom_vline(aes(xintercept = model$cutoff.sd), colour = "red") +
  geom_hline(aes(yintercept = model$cutoff.od), colour = "red") +
  xlab("Score distance") +
  ylab("Orthogonal distance") +
  scale_colour_manual(
    "Food group",
    values = setNames(brewer.pal(length(incl.legend), "Set1"), incl.legend)
  ) +
  guides(colour = guide_legend(nrow = 2)) +
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.background = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.box.spacing = unit(0, "mm")
  )
ggsave(file.path(plot.dir, "q2-outlier-plot.eps"),
  p,
  width = height / 2,
  height = width / 2,
  units = "mm"
)
