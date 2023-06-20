# This file contains setup code common to all scripts
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(stargazer);
  library(datasets); library(robustbase); library(AER);
  library(corrplot); library(trafo); library(GGally);
  library(fitdistrplus); library(reshape2); library(cellWise);
  library(R.devices); library(fungible) ; library(DescTools);
  library(gganimate); library(rrcov); library(usdanutrients);
  library(dplyr); library(gghighlight); library(ggrepel);
  library(cluster); library(kmed); library(patchwork);
  library(latex2exp); library(nortest); library(xtable)
})

## Directories
plot.dir <- file.path("results", "plots")
tab.dir <- file.path("results", "tables")

## Page width and height in mm
width <- 430 * 0.3528
height <- 556 * 0.3528

## ggplot theme settings
theme_set(theme_bw())

theme_update(
  legend.position = "top",
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  strip.text.x = element_text(size = 8)
)
