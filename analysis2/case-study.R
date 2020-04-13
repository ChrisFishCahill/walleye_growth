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
plan(multisession, workers = future::availableCores() / 2)
TMB::compile("analysis2/vb_alta.cpp")

#---------------------
# Read in the data and add in the block variable for cross validation

data <- readRDS("analysis2/vB_analysis_august_2019_cahill.rds")
data[, c("X_TTM_c", "Y_TTM_c")] <- data[, c("X_TTM_c", "Y_TTM_c")] / 1000 # put distance in km
Loc <- unique(data[, c("X_TTM_c", "Y_TTM_c")])

# Add in the block variable to data for cross validation
Loc$Block <- NA
Loc[which(Loc$X_TTM_c < 400 & Loc$Y_TTM_c > 6400), "Block"] <- 1
Loc[which(Loc$X_TTM_c < 410 & Loc$Y_TTM_c < 6200 & is.na(Loc$Block)), "Block"] <- 2
Loc[which(Loc$X_TTM_c < 450 & Loc$Y_TTM_c < 6000 & is.na(Loc$Block)), "Block"] <- 2
Loc[which(Loc$X_TTM_c < 500 & Loc$Y_TTM_c > 6100 & is.na(Loc$Block)), "Block"] <- 3
Loc[which(Loc$X_TTM_c < 660 & Loc$Y_TTM_c > 6000 & is.na(Loc$Block)), "Block"] <- 4
Loc[which(Loc$X_TTM_c < 625 & Loc$Y_TTM_c < 6000 & is.na(Loc$Block)), "Block"] <- 5
Loc[which(Loc$Y_TTM_c < 5800 & is.na(Loc$Block)), "Block"] <- 6
Loc[which(Loc$X_TTM_c > 660 & Loc$Y_TTM_c > 5800 & is.na(Loc$Block)), "Block"] <- 7

# ggplot(Loc, aes(x = X_TTM_c, y = Y_TTM_c, color = as.factor(Block))) +
#   geom_point() +
#   theme_bw() +
#   scale_colour_manual(values = RColorBrewer::brewer.pal(7, "Dark2"))

data <- dplyr::left_join(data, Loc, by = c("X_TTM_c", "Y_TTM_c"))

#---------------------
# Build mesh and inputs for INLA approach:

# mesh = inla.mesh.2d(loc=Loc, max.edge=c(62,1000)) #Better mesh, but slower
mesh <- inla.mesh.create(Loc[, c("X_TTM_c", "Y_TTM_c")], refine = TRUE, extend = -0.5, cutoff = 0.01) # faster mesh

# png(file="Mesh.png",width=9.50,height=7.00,units="in",res=600)
# plot(mesh)
# points(Loc, col = "Steelblue", pch = 1)
# dev.off()

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

#---------------------
# set up fit function for mapping, purrr, furrr

get_fit <- function(Linf = 55, T0 = -1, SigO = 1.0, sd = 0.3,
                    rho = 0.9, kappa = 0.9, ln_sd_linf = 0, ln_sd_tzero = 0,
                    ln_b_sex = 0, ln_sd_omega_lake = 0, ln_sd_omega_time = 0,
                    sig_varies_fitted = c("fixed", "by lake", "by time", "both", "ar1 st"),
                    silent = TRUE, partition_i = NULL,
                    REML = TRUE, ...) {
  sig_varies_fitted <- match.arg(sig_varies_fitted)
  cat(
    crayon::green(
      clisymbols::symbol$tick
    ),
    fitted = "model fitted = ", sig_varies_fitted,
    sep = ""
  )
  # estimate on all data unless partition declared
  if (is.null(partition_i)) {
    partition_i <- rep(0, nrow(data))
  }
  data <- list(
    Nobs = nrow(data),
    length_i = data$TL,
    age_i = data$Age,
    lake_i = data$Lake - 1L,
    time_i = data$Year - 1L,
    Nlakes = length(unique(data$Lake)),
    sex_i = data$SexCode,
    X_ij_omega = model.matrix(~ 1 + data$wallEffDen.Std + data$compEffDen.Std +
      data$GDD.Std + data$wallEffDen.Std:data$compEffDen.Std),
    spdeMatrices = spdeMatrices,
    predTF_i = partition_i
  )

  parameters <- list(
    ln_global_linf = log(Linf),
    ln_sd_linf = ln_sd_linf,
    global_tzero = T0,
    ln_sd_tzero = ln_sd_tzero,
    ln_b_sex = ln_b_sex,
    b_j_omega = rep(0, ncol(data$X_ij_omega)),
    ln_sd_omega_lake = ln_sd_omega_lake,
    eps_omega_lake = rep(0, data$Nlakes),
    ln_sd_omega_time = ln_sd_omega_time,
    eps_omega_time = rep(0, length(unique(data$time_i))),
    eps_linf = rep(0, data$Nlakes),
    eps_t0 = rep(0, data$Nlakes),
    eps_omega_st = matrix(0, nrow = mesh$n, ncol = length(unique(data$time_i))),
    log_sd = log(sd),
    ln_kappa = log(kappa),
    ln_tau_O = log(SigO),
    rho_unscaled = qlogis((rho + 1) / 2)
  )
  map <- list()
  if (sig_varies_fitted %in% c("fixed", "by time", "ar1 st")) {
    map <- c(map, list(
      eps_omega_lake = as.factor(rep(NA, data$Nlakes)),
      ln_sd_omega_lake = factor(NA)
    ))
  }
  if (sig_varies_fitted %in% c("fixed", "by lake", "ar1 st")) {
    map <- c(map, list(
      eps_omega_time = as.factor(rep(NA, length(unique(data$time_i)))),
      ln_sd_omega_time = factor(NA)
    ))
  }
  if (sig_varies_fitted != "ar1 st") {
    map <- c(map, list(
      eps_omega_st = as.factor(matrix(NA, nrow = mesh$n, ncol = length(unique(data$time_i)))),
      ln_kappa = factor(NA),
      rho_unscaled = factor(NA),
      ln_tau_O = factor(NA)
    ))
  }

  random <- c("eps_linf", "eps_t0", "eps_omega_lake", "eps_omega_time", "eps_omega_st")
  if (REML == TRUE) {
    random <- union(random, c(
      "b_j_omega",
      "ln_global_linf",
      "global_tzero",
      "ln_b_sex"
    ))
  }
  if (!"vb_alta" %in% names(getLoadedDLLs())) {
    cat(crayon::blue(clisymbols::symbol$star), "Loading DLL\n")
    dyn.load(dynlib("analysis2/vb_alta"))
  }

  obj <- TMB::MakeADFun(data, parameters,
    DLL = "vb_alta",
    random = random,
    map = map,
    silent = silent
  )
  opt <- TMBhelper::fit_tmb(
    obj = obj,
    control = list(eval.max = 1000, iter.max = 1000),
    getsd = F, newtonsteps = 1 # newtonsteps required for ML model convergence
  )
  rep <- TMB::sdreport(obj, bias.correct = TRUE)
  final_gradient <- obj$gr(opt$par)
  if (any(abs(final_gradient) > 0.001) || rep$pdHess == FALSE) {
    opt$Convergence_check <- "Model did not converge: check results"
    cat(
      crayon::red(
        clisymbols::symbol$cross
      ),
      sig_varies_fitted, "did not converge: check results",
      sep = " "
    )
  }
  list(opt = opt, rep = rep, obj = obj)
}

# Testing:
# test <- get_fit(sig_varies_fitted = "by lake", silent = F, REML = F)
# out <- purrr::pmap(tofit, get_fit, silent = F) %>%
#                    setNames( c("by lake", "by time", "both", "ar1 st"))

#---------------------
# Fit models using restricted maximum likelihood & maximum likelihood

tofit <- dplyr::tibble(
  sig_varies_fitted = c("by lake", "by time", "both", "ar1 st")
)

system.time({ # ~3 minutes
  out <- furrr::future_pmap(tofit, get_fit,
    silent = TRUE,
    REML = TRUE
  ) %>%
    setNames(c("by lake", "by time", "both", "ar1 st"))
})
saveRDS(out, file = "analysis2/REML_fits.rds")
reml <- readRDS(file = "analysis2/REML_fits.rds")

system.time({ # ~10 minutes
  out <- furrr::future_pmap(tofit, get_fit,
    silent = TRUE,
    REML = FALSE
  ) %>%
    setNames(c("by lake", "by time", "both", "ar1 st"))
})
saveRDS(out, file = "analysis2/ML_fits.rds")
ml <- readRDS(file = "analysis2/ML_fits.rds")

#---------------------
# cross-validation routines

run_cv_experiment <- function(which_experiment = c("h block", "lolo"),
                              cv_fold = cv_fold,
                              sig_varies_fitted = sig_varies_fitted) {
  which_experiment <- match.arg(which_experiment)
  if (which_experiment == "h block") {
    partition_i <- ifelse(data$Block != cv_fold, 0, 1)
  } else {
    partition_i <- ifelse(data$Lake == cv_fold, 1, 0)
  }
  out <- get_fit(sig_varies_fitted = sig_varies_fitted, partition_i = partition_i, silent = TRUE)
  tibble::tibble(
    sig_varies_fitted = sig_varies_fitted,
    pred_jnll = out$obj$report()$pred_jnll,
    cv_score = out$obj$report()$pred_jnll / sum(partition_i),
    cv_fold = cv_fold,
    convergence = out$opt$Convergence_check
  )
}

#---------------------
# h-block (spatial) cross validation

tofit <- tidyr::expand_grid(
  cv_fold = seq_len(length(unique(data$Block))),
  sig_varies_fitted = c("by lake", "by time", "both", "ar1 st")
)

system.time({ # 9 minutes
  out <- furrr::future_pmap_dfr(tofit, run_cv_experiment, which_experiment = "h block")
})

saveRDS(out, file = "analysis2/cv_h_block.rds")
cv_h_block <- readRDS("analysis2/cv_h_block.rds")

unique(cv_h_block$convergence)
cv_h_block %>%
  group_by(sig_varies_fitted) %>%
  summarize(h_block_score = sum(cv_score) / n_distinct(cv_fold))

#---------------------
# leave one lake out cross validation

tofit <- tidyr::expand_grid(
  cv_fold = seq_len(length(unique(data$Lake))),
  sig_varies_fitted = c("by lake", "by time", "both", "ar1 st")
)

system.time({ # 113 minutes
  out <- furrr::future_pmap_dfr(tofit, run_cv_experiment, which_experiment = "lolo")
})

saveRDS(out, file = "analysis2/cv_lolo.rds")
cv_lolo <- readRDS("analysis2/cv_lolo.rds")

unique(cv_lolo$convergence)
cv_lolo %>%
  group_by(sig_varies_fitted) %>%
  summarize(cv_lolo_score = sum(cv_score) / n_distinct(cv_fold))

#---------------------
