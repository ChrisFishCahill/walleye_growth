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

Loc <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] ) / 1000 #Put distance in kms

#mesh = inla.mesh.2d(loc=loc_xy, max.edge=c(62,1000)) #Better mesh, but slower
mesh = inla.mesh.create( Loc, refine=TRUE, extend=-0.5, cutoff=0.01 ) #faster mesh

mesh$n

#png(file="Mesh.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(Loc, col="Steelblue", pch=1)

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

rho_sd_prior = 2
rho_mean_prior = 0
tau_O_mean_prior = 0
tau_O_sd_prior = 3

data <- list(
  Nobs = nrow(data),
  length_i = data$TL,
  age_i = data$Age,
  lake_i = data$Lake - 1L,
  time_i = data$Year - 1L,
  Nlakes = length(unique(data$Lake)),
  spdeMatrices = spdeMatrices,
  rho_sd_prior = rho_sd_prior,
  rho_mean_prior = rho_mean_prior,
  tau_O_mean_prior = tau_O_mean_prior,
  tau_O_sd_prior = tau_O_sd_prior
)

parameters <- list(
  ln_global_linf = log(55),
  ln_sd_linf = 0,
  global_tzero = -1,
  ln_sd_tzero = 0,
  ln_global_omega = log(15),
  ln_sd_omega_lake = 0,
  ln_sd_omega_time = 0,
  eps_omega_lake = rep(0, data$Nlakes),
  eps_omega_time = rep(0, length(unique(data$time_i))),
  eps_linf = rep(0, data$Nlakes),
  eps_t0 = rep(0, data$Nlakes),
  eps_omega_st = matrix(0, nrow = mesh$n, ncol = length(unique(data$time_i))),
  ln_cv = log(0.069),
  ln_kappa = log(3),
  ln_tau_O = 0,
  rho_unscaled = 3.4
)

map <- list(
  eps_omega_time = as.factor(rep(NA, length(unique(data$time_i)))),
  ln_sd_omega_time = factor(NA),
    eps_omega_lake = as.factor(rep(NA, data$Nlakes)),
    ln_sd_omega_lake = factor(NA)
  )

dyn.load(dynlib("analysis2/vb_alta"))

obj <- TMB::MakeADFun(data, parameters,
                      DLL = "vb_alta",
                      random = c("eps_omega_lake", "eps_omega_time", "eps_omega_st", "eps_linf", "eps_t0"),
                      map = map
)

opt = nlminb(obj$par, obj$fn, obj$gr,
             eval.max = 1000, iter.max = 1000)
sdrep = sdreport( obj )
final_gradient = obj$gr( opt$par )
if( any(abs(final_gradient)>0.01) || sdrep$pdHess==FALSE ) stop("Not converged")

rangeIndex = which(row.names(summary(sdrep,"report"))=="range")
range = summary(sdrep,"report")[rangeIndex,]



rep <- sdreport(obj, bias.correct=T)



