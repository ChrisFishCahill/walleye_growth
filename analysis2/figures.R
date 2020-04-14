#---------------------
# Load packages

library(ggplot2)
theme_set(theme_light())
library(TMB)
library(future)
library(tidyr)
library(INLA)
library(purrr)
library(furrr)
library(dplyr)
library(TMBhelper)
library(arm)
library(gtable)
library(grid)
library(gridExtra)

#---------------------
# Read in results and make Table 1 (misery):

cv_h_block <- readRDS("analysis2/cv_h_block.rds")

unique(cv_h_block$convergence)
MyTable <- cv_h_block %>%
  group_by(sig_varies_fitted, fit_interaction) %>%
  summarize(h_block_score = sum(cv_score) / n_distinct(cv_fold))

cv_lolo <- readRDS("analysis2/cv_lolo.rds")

unique(cv_lolo$convergence)
MyTable <- cv_lolo %>%
  group_by(sig_varies_fitted, fit_interaction) %>%
  summarize(cv_lolo_score = sum(cv_score) / n_distinct(cv_fold)) %>%
  dplyr::left_join(MyTable, .x, by = c("sig_varies_fitted", "fit_interaction"))

MyTable$`REML AIC` <- NA
MyTable$`ML AIC` <- NA

# Add in the AIC values:
ml <- readRDS(file = "analysis2/ML_fits.rds")
reml <- readRDS(file = "analysis2/REML_fits.rds")

remls <- map_dbl(.x = reml, .f = function(.x) {
  return(.x$opt$AIC)
})[1:nrow(MyTable)] # get rid of random slopes model
mls <- map_dbl(.x = ml, .f = function(.x) {
  return(.x$opt$AIC)
})[1:nrow(MyTable)]

which_order <- c(
  "ar1 st reduced", "ar1 st full", "both reduced",
  "both full", "by lake reduced", "by lake full",
  "by time reduced", "by time full"
)

MyTable$`REML AIC` <- remls[order(factor(names(remls), levels = which_order))]
MyTable$`ML AIC` <- mls[order(factor(names(mls), levels = which_order))]

MyTable$fit_interaction <- ifelse(MyTable$fit_interaction == "TRUE", "Yes", "No")

MyTable$`REML AIC` <- round(MyTable$`REML AIC`, 1)
MyTable$`REML AIC` <- MyTable$`REML AIC` - min(MyTable$`REML AIC`)
MyTable$`REML AIC` <- sprintf("%0.1f", MyTable$`REML AIC`)

MyTable$`ML AIC` <- round(MyTable$`ML AIC`, 1)
MyTable$`ML AIC` <- MyTable$`ML AIC` - min(MyTable$`ML AIC`)
MyTable$`ML AIC` <- sprintf("%0.1f", MyTable$`ML AIC`)

MyTable$cv_lolo_score <- round(MyTable$cv_lolo_score, 2)
MyTable$cv_lolo_score <- sprintf("%0.2f", MyTable$cv_lolo_score)
MyTable$h_block_score <- round(MyTable$h_block_score, 2)
MyTable$h_block_score <- sprintf("%0.2f", MyTable$h_block_score)

colnames(MyTable) <- c(
  "bold(Model)", "bold(Interaction~term)",
  "bold(LOLO~CV)", "bold(H~block~CV)", paste0("bold(Delta", "~REML~AIC)"),
  paste0("bold(Delta", "~ML~AIC)")
)

g <- tableGrob(MyTable,
  rows = NULL,
  theme = ttheme_minimal(
    parse = FALSE,
    colhead = list(fg_params = list(fontsize = 12, fontface = 3, parse = TRUE))
  )
)

g <- gtable_add_rows(g, heights = unit(0.15, "npc"), pos = 0)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 1, b = 1, l = 1, r = 6
)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 2, b = 2, l = 1, r = 6
)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 10, b = 10, l = 1, r = 6
)

#ggsave("analysis2/table_1.jpeg", g, width = 7, height = 4.5, dpi=600, units="in")

#---------------------
