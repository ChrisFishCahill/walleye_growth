#TODO--is there a smarter way to pass the mesh from get_sim_data to fit_sim?
#TODO--convergence checks?
#TODO--fixed estimation does not work
library(ggplot2)
theme_set(theme_light())
library(TMB)
library(future)
library(tidyr)
library(INLA)
library(purrr)
source("sim2/INLA_helpers.R")

plan(multisession, workers = future::availableCores()/2)

get_sim_data <- function(Nyears = 10, Nlakes = 12, Nfish = 25,
                         Linf = 55, T0 = -1, SigO = 0.8, cv = 0.2, omega_global = 14,
                         rho = 0.5, kappa = 0.5,
                         sig_varies = c("fixed", "by lake", "by time", "both", "ar1")) {
  sig_varies <- match.arg(sig_varies)
  omega_dev_st <- matrix(0, nrow = Nlakes, ncol = Nyears) # omega_dev_st all zero unless "ar1" selected
  Loc <- cbind("x" = runif(Nlakes, min = 0, max = 10), "y" = runif(Nlakes, min = 0, max = 10))
  if (sig_varies == "fixed") {
    omega_dev_lake <- rnorm(Nlakes, 0, 0)
    omega_dev_time <- rnorm(Nyears, 0, 0)
  } else if (sig_varies == "by lake") {
    omega_dev_lake <- rnorm(Nlakes, 0, SigO)
    omega_dev_time <- rnorm(Nyears, 0, 0)
  } else if (sig_varies == "by time") {
    omega_dev_lake <- rnorm(Nlakes, 0, 0)
    omega_dev_time <- rnorm(Nyears, 0, SigO)
  } else if (sig_varies == "both") {
    omega_dev_lake <- rnorm(Nlakes, 0, SigO)
    omega_dev_time <- rnorm(Nyears, 0, SigO)
  } else if (sig_varies == "ar1") {
    # set lake and year devs off:
    omega_dev_lake <- rnorm(Nlakes, 0, 0)
    omega_dev_time <- rnorm(Nyears, 0, 0)
    # simulate space-time devs a la INLA/GMRFlib:
    mesh <- inla.mesh.create(Loc, refine = TRUE, extend = -0.5, cutoff = 0.01)
    omega_dev_k <- rspde(Loc,
      range = sqrt(8) / kappa,
      sigma = SigO, n = Nyears, mesh = mesh,
      return.attributes = TRUE, seed = sample.int(1e6, 1)
    )
    omega_dev_st <- omega_dev_k[1:Nlakes, 1:Nyears]
    for (j in 2:Nyears) {
      omega_dev_st[, j] <- rho * omega_dev_st[, j - 1] + sqrt(1 - rho^2) * omega_dev_k[, j]
    }
  }

  to_sim <- tidyr::expand_grid(lake = 1:Nlakes, year = 1:Nyears)
  purrr::map2_dfr(to_sim$lake, to_sim$year, function(lake, year) {
    ages <- sample(0:25, Nfish, replace = T)
    eta_it <- exp(log(omega_global) + omega_dev_time[year] + omega_dev_lake[lake] + omega_dev_st[lake, year])
    lpreds <- Linf * (1 - exp(-(eta_it / Linf) * (ages - T0)))
    which_x <- Loc[lake, 1]
    which_y <- Loc[lake, 2]
    which_omega_dev_st <- omega_dev_st[lake, year]
    y_i <- rlnorm(Nfish, log(lpreds), cv)
    tibble::tibble(y_i, ages,
      lake = lake, year = year,
      linf = Linf, t0 = T0, omega_global = omega_global,
      rho = rho, kappa = kappa, SigO = SigO,
      x = rep(which_x, Nfish), y = rep(which_y, Nfish),
      omega_dev_st = rep(which_omega_dev_st, Nfish)
    )
  })
}

out <- purrr::map_dfr(seq_len(5), function(x) {
  get_sim_data(Nlakes = 5, sig_varies = "by lake")
}, .id = "sim_iter")
ggplot(out, aes(ages, y_i)) +
  facet_grid(sim_iter ~ lake) +
  geom_point(alpha = 0.2)

out <- get_sim_data(Nlakes = 5, Nyears = 7, Nfish = 100, cv = 0.01, sig_varies = "by time")
ggplot(out, aes(ages, y_i)) +
  facet_grid(year ~ lake) +
  geom_point(alpha = 0.2)

out <- get_sim_data(Nlakes = 5, Nyears = 7, Nfish = 100, cv = 0.01, sig_varies = "by lake")
ggplot(out, aes(ages, y_i)) +
  facet_grid(year ~ lake) +
  geom_point(alpha = 0.2)

out <- get_sim_data(Nlakes = 5, Nyears = 7, Nfish = 100, cv = 0.01, sig_varies = "both")
ggplot(out, aes(ages, y_i)) +
  facet_grid(year ~ lake) +
  geom_point(alpha = 0.2)

out <- get_sim_data(Nlakes = 50, Nyears = 7, Nfish = 100, cv = 0.01, rho = 0.5, kappa = 0.5, sig_varies = "ar1")

ggplot(out, aes(x, y, col = omega_dev_st)) +
  geom_point() +
  facet_wrap(~year) +
  scale_color_gradient2()

out <- get_sim_data(Nlakes = 7, Nyears = 5, Nfish = 100, cv = 0.01, rho = 0.5, kappa = 0.5, sig_varies = "ar1")
ggplot(out, aes(ages, y_i)) +
  geom_point() +
  facet_wrap(~lake)

out <- purrr::map_dfr(seq_len(5), function(x) {
  get_sim_data(Nlakes = 7, sig_varies = "ar1")
}, .id = "sim_iter")
ggplot(out, aes(ages, y_i)) +
  facet_grid(sim_iter ~ lake) +
  geom_point(alpha = 0.2)

TMB::compile("sim2/vb_cyoa.cpp")

fit_sim <- function(Nyears = 10, Nlakes = 12, Nfish = 20,
                    Linf = 55, T0 = -1, SigO = 0.8, cv = 0.2, omega_global = 14,
                    rho = 0.5, kappa = 0.5,
                    sig_varies = c("fixed", "by lake", "by time", "both", "ar1"),
                    sig_varies_fitted = c("fixed", "by lake", "by time", "both", "ar1"),
                    iter = NA) {
  # browser()
  #sig_varies <- "ar1"
  #sig_varies_fitted <- "fixed"
  #iter=1
  sig_varies <- match.arg(sig_varies)

  sim_dat <- get_sim_data(
    Nyears = Nyears, Nlakes = Nlakes, Nfish = Nfish,
    Linf = Linf, T0 = T0, SigO = SigO, cv = cv, omega_global = omega_global,
    rho = rho, kappa = kappa,
    sig_varies = sig_varies
  )

  Loc <- unique(sim_dat[, c("x", "y")])

  mesh <- inla.mesh.create(Loc, refine = TRUE, extend = -0.5, cutoff = 0.01)
  spde <- inla.spde2.matern(mesh, alpha = 2)
  spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

  data <- list(
    Nobs = nrow(sim_dat),
    length_i = sim_dat$y_i,
    age_i = sim_dat$ages,
    lake_i = sim_dat$lake - 1L,
    time_i = sim_dat$year - 1L,
    Nlakes = length(unique(sim_dat$lake)),
    spdeMatrices = spdeMatrices
  )
  parameters <- list(
    ln_global_linf = log(sim_dat$linf[1]),
    ln_sd_linf = 0,
    global_tzero = sim_dat$t0[1],
    ln_sd_tzero = 0,
    ln_global_omega = log(sim_dat$omega_global[1]),
    ln_sd_omega_lake = 0,
    ln_sd_omega_time = 0,
    eps_omega_lake = rep(0, data$Nlakes),
    eps_omega_time = rep(0, length(unique(sim_dat$year))),
    eps_linf = rep(0, data$Nlakes),
    eps_t0 = rep(0, data$Nlakes),
    eps_omega_st = matrix(0, nrow = mesh$n, ncol = Nyears),
    ln_cv = 0,
    ln_kappa = log(kappa),
    ln_tau_O = 0,
    rho_unscaled = 2 * plogis(rho) - 1
  )
  map <- list(
    ln_sd_tzero = factor(NA),
    ln_sd_linf = factor(NA),
    eps_linf = as.factor(rep(NA, data$Nlakes)),
    eps_t0 = as.factor(rep(NA, data$Nlakes))
  )
  if (sig_varies_fitted %in% c("fixed", "by lake")) {
    map <- c(map, list(
      eps_omega_time = as.factor(rep(NA, length(unique(sim_dat$year)))),
      ln_sd_omega_time = factor(NA)
    ))
  }
  if (sig_varies_fitted %in% c("fixed", "by time")) {
    map <- c(map, list(
      eps_omega_lake = as.factor(rep(NA, length(unique(sim_dat$lake)))),
      ln_sd_omega_lake = factor(NA)
    ))
  }
  if (sig_varies_fitted != "ar1") {
    map <- c(map, list(
      eps_omega_st = as.factor(matrix(NA, nrow = mesh$n, ncol = Nyears)),
      ln_kappa = factor(NA),
      rho_unscaled = factor(NA)
    ))
  }

  #sink(tempfile())
  dyn.load(dynlib("sim2/vb_cyoa"))
  obj <- TMB::MakeADFun(data, parameters,
    DLL = "vb_cyoa",
    random = c("eps_omega_lake", "eps_omega_time", "eps_omega_st", "eps_linf", "eps_t0"),
    map = map,
    silent = TRUE
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  dyn.unload(dynlib("sim2/vb_cyoa"))
  #sink()
  tibble::tibble(
    sig_varies = sig_varies,
    sig_varies_fitted = sig_varies_fitted,
    ln_global_omega = opt$par[["ln_global_omega"]],
    true_ln_global_omega = log(sim_dat$omega_global[1]),
    iter = iter
  )
}

totest <- tidyr::expand_grid(
  iter = seq_len(10L),
  sig_varies = c("by lake", "by time", "both", "ar1"),
  sig_varies_fitted = c("by lake", "by time", "both", "ar1")
)
nrow(totest)

set.seed(1234)
out <- purrr::pmap_dfr(totest, fit_sim) # for testing

system.time({
  out <- furrr::future_pmap_dfr(totest, fit_sim)
})

saveRDS(out, file = "sim2/sim2.rds")
out <- readRDS("sim2/sim2.rds")
out %>%
  dplyr::mutate(sig_varies = paste0("Sim = ", sig_varies)) %>%
  dplyr::mutate(sig_varies_fitted = paste0("Fitted = ", sig_varies_fitted)) %>%
  ggplot(aes(ln_global_omega)) +
  geom_histogram(bins = 25) +
  geom_vline(xintercept = out[["true_ln_global_omega"]][1]) +
  facet_grid(sig_varies_fitted ~ sig_varies) +
  xlab(expression(omega)) +
  ylab("Count")
ggsave("sim2/hist-sim.pdf", width = 7, height = 5)

out %>%
  dplyr::mutate(Matched = ifelse(sig_varies_fitted == sig_varies, TRUE, FALSE)) %>%
  ggplot(aes(sig_varies, exp(ln_global_omega), colour = sig_varies_fitted, fill = Matched)) +
  geom_boxplot() +
  geom_hline(yintercept = exp(out[["true_ln_global_omega"]][1])) +
  xlab(expression(Simulated ~ omega ~ variation)) +
  labs(colour = expression(Fitted ~ omega ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(omega[0]))
ggsave("sim2/boxplot-sim.pdf", width = 7, height = 5)
