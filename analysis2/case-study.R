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

data <- readRDS("analysis2/vB_analysis_august_2019_cahill.rds")
TMB::compile("analysis2/vb_alta.cpp")

Loc <- unique(data[, c("X_TTM_c", "Y_TTM_c")]) / 1000 # Put distance in kms

# mesh = inla.mesh.2d(loc=Loc, max.edge=c(62,1000)) #Better mesh, but slower
mesh <- inla.mesh.create(Loc, refine = TRUE, extend = -0.5, cutoff = 0.01) # faster mesh

# png(file="Mesh.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(Loc, col = "Steelblue", pch = 1)

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

get_fit <- function(Linf = 55, T0 = -1, SigO = 1.0, cv = 0.05, omega_global = 15,
                    rho = 0, kappa = 10.0,
                    sig_varies_fitted = c("fixed", "by lake", "by time", "both", "ar1 st"),
                    silent = TRUE, Partition_i = rep(0, nrow(data)),
                    REML = 1L) {
  sig_varies_fitted <- match.arg(sig_varies_fitted)
  cat(
    crayon::green(
      clisymbols::symbol$tick
    ),
    fitted = "model fitted = ", sig_varies_fitted,
    sep = ""
  )

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
    predTF_i = Partition_i
  )

  parameters <- list(
    ln_global_linf = log(Linf),
    ln_sd_linf = 0,
    global_tzero = T0,
    ln_sd_tzero = 0,
    ln_b_sex = 0,
    b_j_omega = rep(0, ncol(data$X_ij_omega)),
    ln_sd_omega_lake = 0,
    eps_omega_lake = rep(0, data$Nlakes),
    ln_sd_omega_time = 0,
    eps_omega_time = rep(0, length(unique(data$time_i))),
    eps_linf = rep(0, data$Nlakes),
    eps_t0 = rep(0, data$Nlakes),
    eps_omega_st = matrix(0, nrow = mesh$n, ncol = length(unique(data$time_i))),
    ln_cv = log(cv),
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

  opt <- nlminb(obj$par, obj$fn, obj$gr,
    eval.max = 1000, iter.max = 1000
  )
  rep <- TMB::sdreport(obj, bias.correct=TRUE)
  final_gradient <- obj$gr(opt$par)
  if (any(abs(final_gradient) > 0.01) || rep$pdHess == FALSE) {
    browser()
    opt$convergence <- 1L
  }
  AIC <- (2 * length(opt$par) - 2 * (-opt$objective))
  list(opt = opt, rep = rep, AIC = AIC, obj = obj)
}

# Testing:
by_lake <- get_fit(sig_varies_fitted = "by lake", cv = 0.15, SigO = 0.5, silent = F)
by_time <- get_fit(sig_varies_fitted = "by time", cv = 0.15, SigO = 2.0, silent = F)
both <- get_fit(sig_varies_fitted = "both", cv = 0.5, SigO = 1.38, silent = F)
ar1 <- get_fit(sig_varies_fitted = "ar1 st", SigO = 1.38, cv = 0.05, kappa = 0.05371844, silent = F)

by_lake$rep$value
ar1$obj$report()
ar1$obj$report()

#TODO: why is adreport broken
#some alain zuur stuff here and make a pretty table
#Model, Number of parameters, Objective, AIC

str(ar1$obj$report())

ar1$AIC
by_lake$AIC
by_time$AIC
both$AIC
