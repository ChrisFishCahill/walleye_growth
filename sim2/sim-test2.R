library(ggplot2)
theme_set(theme_light())
library(TMB)
library(future)
library(tidyr)
library(INLA)
library(purrr)
library(furrr)
library(dplyr)
library(fs)
devtools::install_github("seananderson/ggsidekick")

# helper functions from http://www.r-inla.org/spde-book
source("sim2/INLA_helpers.R")

plan(multisession, workers = future::availableCores() / 2)

get_sim_data <- function(Nyears = 10, Nlakes = 15, Nfish = 20,
                         Linf = 55, T0 = -1, SigO = 0.8, cv = 0.2, omega_global = 14,
                         rho = 0.5, kappa = 0.5,
                         sig_varies = c("fixed", "by lake", "by time", "both", "ar1 st")) {
  sig_varies <- match.arg(sig_varies)
  omega_dev_st <- matrix(0, nrow = Nlakes, ncol = Nyears)
  Loc <- cbind("x" = runif(Nlakes, min = 0, max = 10), "y" = runif(Nlakes, min = 0, max = 10))
  mesh <- INLA::inla.mesh.create(Loc, refine = TRUE, extend = -0.5, cutoff = 0.01)
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
  } else if (sig_varies == "ar1 st") {
    omega_dev_lake <- rnorm(Nlakes, 0, 0)
    omega_dev_time <- rnorm(Nyears, 0, 0)
    # simulate space-time devs a la INLA/GMRFlib:
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
  out <- purrr::map2_dfr(to_sim$lake, to_sim$year, function(lake, year) {
    ages <- sample(0:25, Nfish, replace = TRUE)
    eta_it <- exp(log(omega_global) + omega_dev_time[year] +
      omega_dev_lake[lake] + omega_dev_st[lake, year])
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
  list(dat = out, mesh = mesh)
}

# out <- purrr::map_dfr(seq_len(5), function(x) {
#   get_sim_data(Nlakes = 5, sig_varies = "by lake")$dat
# }, .id = "sim_iter")
# ggplot(out, aes(ages, y_i)) +
#   facet_grid(sim_iter ~ lake) +
#   geom_point(alpha = 0.2)
#
# out <- get_sim_data(
#   Nlakes = 5, Nyears = 7,
#   Nfish = 100, cv = 0.01, sig_varies = "by time"
# )$dat
# ggplot(out, aes(ages, y_i)) +
#   facet_grid(year ~ lake) +
#   geom_point(alpha = 0.2)
#
# out <- get_sim_data(
#   Nlakes = 5, Nyears = 7,
#   Nfish = 100, cv = 0.01, sig_varies = "by lake"
# )$dat
# ggplot(out, aes(ages, y_i)) +
#   facet_grid(year ~ lake) +
#   geom_point(alpha = 0.2)
#
# out <- get_sim_data(
#   Nlakes = 5, Nyears = 7,
#   Nfish = 100, cv = 0.01, sig_varies = "both"
# )$dat
# ggplot(out, aes(ages, y_i)) +
#   facet_grid(year ~ lake) +
#   geom_point(alpha = 0.2)
#
# out <- get_sim_data(
#   Nlakes = 50, Nyears = 7,
#   Nfish = 100, cv = 0.01, rho = 0.5, kappa = 0.5,
#   sig_varies = "ar1 st"
# )$dat
#
# ggplot(out, aes(x, y, col = omega_dev_st)) +
#   geom_point() +
#   facet_wrap(~year) +
#   scale_color_gradient2()
#
# out <- get_sim_data(
#   Nlakes = 7, Nyears = 5,
#   Nfish = 100, cv = 0.01, rho = 0.5, kappa = 0.5,
#   sig_varies = "ar1 st"
# )$dat
# ggplot(out, aes(ages, y_i)) +
#   geom_point() +
#   facet_wrap(~lake)
#
# out <- purrr::map_dfr(seq_len(5), function(x) {
#   get_sim_data(Nlakes = 7, sig_varies = "ar1 st")$dat
# }, .id = "sim_iter")
# ggplot(out, aes(ages, y_i)) +
#   facet_grid(sim_iter ~ lake) +
#   geom_point(alpha = 0.2)

TMB::compile("sim2/vb_cyoa.cpp")

fit_sim <- function(Nyears = 10, Nlakes = 15, Nfish = 20,
                    Linf = 55, T0 = -1, SigO = 0.8, cv = 0.2, omega_global = 14,
                    rho = 0.5, kappa = 0.5,
                    sig_varies = c("fixed", "by lake", "by time", "both", "ar1 st"),
                    sig_varies_fitted = c("fixed", "by lake", "by time", "both", "ar1 st"),
                    iter = NA, silent = TRUE,
                    rho_sd_prior = 2, rho_mean_prior = 0,
                    tau_O_mean_prior = 0, tau_O_sd_prior = 3) {
  sig_varies <- match.arg(sig_varies)
  cat(
    crayon::green(
      clisymbols::symbol$tick
    ),
    " sim = ", sig_varies, "; fitted = ", sig_varies_fitted, "; iter = ", iter, "\n",
    sep = ""
  )

  sim <- get_sim_data(
    Nyears = Nyears, Nlakes = Nlakes, Nfish = Nfish,
    Linf = Linf, T0 = T0, SigO = SigO, cv = cv, omega_global = omega_global,
    rho = rho, kappa = kappa,
    sig_varies = sig_varies
  )
  sim_dat <- sim$dat
  mesh <- sim$mesh
  spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
  spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

  data <- list(
    Nobs = nrow(sim_dat),
    length_i = sim_dat$y_i,
    age_i = sim_dat$ages,
    lake_i = sim_dat$lake - 1L,
    time_i = sim_dat$year - 1L,
    Nlakes = length(unique(sim_dat$lake)),
    spdeMatrices = spdeMatrices,
    rho_sd_prior = rho_sd_prior,
    rho_mean_prior = rho_mean_prior,
    tau_O_mean_prior = tau_O_mean_prior,
    tau_O_sd_prior = tau_O_sd_prior
  )
  parameters <- list(
    ln_global_linf = log(sim_dat$linf[1]) + rnorm(1, 0, 0.1),
    ln_sd_linf = 0,
    global_tzero = sim_dat$t0[1] + rnorm(1, 0, 0.1),
    ln_sd_tzero = 0,
    ln_global_omega = log(sim_dat$omega_global[1]) + rnorm(1, 0, 0.1),
    ln_sd_omega_lake = 0,
    ln_sd_omega_time = 0,
    eps_omega_lake = rep(0, data$Nlakes),
    eps_omega_time = rep(0, length(unique(sim_dat$year))),
    eps_linf = rep(0, data$Nlakes),
    eps_t0 = rep(0, data$Nlakes),
    eps_omega_st = matrix(0, nrow = mesh$n, ncol = Nyears),
    ln_cv = log(cv) + rnorm(1, 0, 0.1),
    ln_kappa = log(kappa) + rnorm(1, 0, 0.1),
    ln_tau_O = log(SigO) + rnorm(1, 0, 0.1), # get close
    rho_unscaled = qlogis((0 + 1) / 2) + rnorm(1, 0, 0.05)
  )
  map <- list(
    ln_sd_tzero = factor(NA),
    ln_sd_linf = factor(NA),
    eps_linf = as.factor(rep(NA, data$Nlakes)),
    eps_t0 = as.factor(rep(NA, data$Nlakes))
  )
  if (sig_varies_fitted %in% c("fixed", "by lake", "ar1 st")) {
    map <- c(map, list(
      eps_omega_time = as.factor(rep(NA, length(unique(sim_dat$year)))),
      ln_sd_omega_time = factor(NA)
    ))
  }
  if (sig_varies_fitted %in% c("fixed", "by time", "ar1 st")) {
    map <- c(map, list(
      eps_omega_lake = as.factor(rep(NA, length(unique(sim_dat$lake)))),
      ln_sd_omega_lake = factor(NA)
    ))
  }
  if (sig_varies_fitted != "ar1 st") {
    map <- c(map, list(
      eps_omega_st = as.factor(matrix(NA, nrow = mesh$n, ncol = Nyears)),
      ln_kappa = factor(NA),
      rho_unscaled = factor(NA),
      ln_tau_O = factor(NA)
    ))
  }

  if (!"vb_cyoa" %in% names(getLoadedDLLs())) {
    cat(crayon::blue(clisymbols::symbol$star), "Loading DLL\n")
    dyn.load(dynlib("sim2/vb_cyoa"))
  }
  obj <- TMB::MakeADFun(data, parameters,
    DLL = "vb_cyoa",
    random = c("eps_omega_lake", "eps_omega_time", "eps_omega_st", "eps_linf", "eps_t0"),
    map = map,
    silent = silent
  )
  opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr, eval.max = 1000, iter.max = 500),
    error = function(e) list(par = list(ln_global_omega = NA, convergence = 1))
  )

  if (sig_varies_fitted == "ar1 st") {
    rho_hat <- 2 * plogis(opt$par[["rho_unscaled"]]) - 1
    ln_global_omega_hat <- opt$par[["ln_global_omega"]]
    # If fails to converge, or omega is stupid high, or rho bumps against the boundary:
    if (opt$convergence != 0 || ln_global_omega_hat > 15 || round(rho_hat, 2) == 1) {
      map <- map$rho_unscaled <- factor(NA)
      # If rho is stuck, fix it and re-estimate:
      if (round(rho_hat, 2) == 1) {
        parameters$rho_unscaled <- qlogis((0.99 + 1) / 2)
      } else {
        parameters$rho_unscaled <- qlogis((0 + 1) / 2)
      }
      obj <- TMB::MakeADFun(data, parameters,
        DLL = "vb_cyoa",
        random = c("eps_omega_lake", "eps_omega_time", "eps_omega_st", "eps_linf", "eps_t0"),
        map = map,
        silent = silent
      )
      opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr, eval.max = 1000, iter.max = 500),
        error = function(e) list(par = list(ln_global_omega = NA, convergence = 1))
      )
    }
  }
  ln_global_omega_hat <- opt$par[["ln_global_omega"]]
  ln_global_linf_hat <- opt$par[["ln_global_linf"]]
  if (is.na(ln_global_omega_hat) || opt$convergence != 0) {
    opt$par[["ln_global_omega"]] <- NA
    opt$convergence <- 1
  } else if (!is.na(ln_global_omega_hat) && ln_global_omega_hat > ln_global_linf_hat) {
    opt$par[["ln_global_omega"]] <- NA
    opt$convergence <- 1
  }

  # dyn.unload(dynlib("sim2/vb_cyoa"))
  tibble::tibble(
    sig_varies = sig_varies,
    sig_varies_fitted = sig_varies_fitted,
    ln_global_omega = opt$par[["ln_global_omega"]],
    true_ln_global_omega = log(sim_dat$omega_global[1]),
    iter = iter, convergence = opt$convergence
  )
}

# totest <- dplyr::tibble(
#   iter = seq_len(1L),
#   sig_varies = c("both"),
#   sig_varies_fitted = c("ar1 st")
# )
# set.seed(1)
# out <- purrr::pmap_dfr(totest, fit_sim, silent = F) # testing
# out <- furrr::future_pmap_dfr(totest, fit_sim,
# .options = future_options(seed = 123L) #for testing parallel

# Visualize the priors (penalties) for ar1 st model:
tau_O_mean_prior <- 0
tau_O_sd_prior <- 3
rho_sd_prior <- 2
rho_mean_prior <- 0

.n <- 200
df <- data.frame(
  parameter = c(rep("tau_O", .n), rep("rho", .n)),
  value = c(seq(0, 10, length.out = .n), seq(-1, 1, length.out = .n)),
  density = c(
    dnorm(seq(0, 10, length.out = .n), mean = tau_O_mean_prior, sd = tau_O_sd_prior),
    dnorm(seq(-1, 1, length.out = .n), mean = rho_mean_prior, sd = rho_sd_prior)
  )
)
df %>%
  ggplot(aes(value, density)) +
  geom_line() +
  facet_wrap(~parameter, scales = "free") +
  scale_y_continuous(limits = c(0, NA))

pars_to_sim1 <- tidyr::expand_grid(
  rho = c(0.1, 0.5, 0.9),
  kappa = c(0.57), # 50% range on 10x10 grid
  SigO = c(0.5)
)
pars_to_sim2 <- tidyr::expand_grid(
  rho = c(0.5),
  kappa = c(2.7, 0.57, 0.28), # 10%, 50%, and 90% range on 10x10 grid
  SigO = c(0.5)
)
pars_to_sim3 <- tidyr::expand_grid(
  rho = c(0.5),
  kappa = c(0.57), # 50% range on 10x10 grid
  SigO = c(0.2, 0.5, 0.8)
)
pars_to_sim <- bind_rows(pars_to_sim1, pars_to_sim2) %>%
  bind_rows(pars_to_sim3)

set.seed(9473772)
pars_to_sim <- pars_to_sim %>%
  tibble::add_column(
    seed = sample.int(1e6, nrow(pars_to_sim)),
    sim = seq_len(nrow(pars_to_sim))
  )

run_sim_experiment <- function(rho = rho, kappa = kappa,
                               SigO = SigO, seed = seed,
                               sim = sim) {
  f <- paste0("sim2/sim_", sim, ".rds")
  if (file.exists(f)) {
    return(NULL)
  } else {
    out <- furrr::future_pmap_dfr(totest, fit_sim,
      rho = rho, kappa = kappa, SigO = SigO,
      .options = furrr::future_options(seed = seed)
    )
    saveRDS(out, file = f)
  }
}

totest <- tidyr::expand_grid(
  iter = seq_len(300L),
  sig_varies = c("by lake", "by time", "both", "ar1 st"),
  sig_varies_fitted = c("by lake", "by time", "both", "ar1 st")
)

if (FALSE) {
  system.time({
    # 4.73 hrs on Cahill's desktop
    purrr::pwalk(pars_to_sim, run_sim_experiment)
  })
}

#-------------------------------------------------------------
# Map through directory, load sim.rds files, make dope plots
#-------------------------------------------------------------

paths <- fs::dir_ls("sim2/", regexp = "^sim2/sim_[0-9]+.rds$")
out <- purrr::map_dfr(paths, readRDS, .id = "sim")

out <- out %>%
  rowwise() %>%
  mutate(sim = strsplit(sim, "/")[[1]][2]) %>%
  mutate(sim = strsplit(sim, ".r")[[1]][1])

buggered <-
  out %>%
  group_by(sim) %>%
  filter(convergence == 1L) %>%
  distinct(iter) %>%
  mutate(failed_iter = iter)

# percent of simulations that failed to converge:
100 * (table(buggered$sim) / 300)

out <- left_join(out, buggered)
summary(out$failed_iter)

types <- tibble(sim = unique(out$sim))
types$param <- rep(c("rho", "spatial~range", "sigma[O]"), each = 3)
types$levels <- rep(c("Low", "Medium", "High"), 3)
out <- left_join(out, types, by = "sim")

out$levels <- factor(out$levels, levels = c("Low", "Medium", "High"))
out$sig_varies_fitted <- factor(out$sig_varies_fitted, levels = c("by lake", "by time", "both", "ar1 st"))
out$sig_varies <- factor(out$sig_varies, levels = c("by lake", "by time", "both", "ar1 st"))

true_omega <- exp(out$true_ln_global_omega)[1]
g <- out %>%
  # out[sample(seq_len(nrow(out)), 6e3), ] %>%
  dplyr::filter(is.na(failed_iter)) %>%
  dplyr::mutate(Matched = sig_varies_fitted == sig_varies) %>%
  ggplot(aes(sig_varies, exp(ln_global_omega),
    colour = sig_varies_fitted, fill = Matched
  )) +
  geom_boxplot(outlier.size = 0.6) +
  facet_grid(vars(param), vars(levels), labeller = label_parsed) +
  geom_hline(yintercept = true_omega, lty = 1, alpha = 0.6) +
  xlab(expression(Simulated ~ omega[0] ~ variation)) +
  labs(colour = expression(Fitted ~ omega[0] ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(Estimated~omega[0])) +
  coord_cartesian(ylim = c(true_omega * 0.5, true_omega * 2)) +
  ggsidekick::theme_sleek()
#ggsave("sim2/simulations-facet-grid.pdf", width = 10, height = 6)

true_omega <- exp(out$true_ln_global_omega)[1]
g <-
  out %>%
  # out[sample(seq_len(nrow(out)), 6e3), ] %>%
  dplyr::filter(is.na(failed_iter)) %>%
  dplyr::mutate(Matched = sig_varies_fitted == sig_varies) %>%
  ggplot(aes(sig_varies, exp(ln_global_omega),
    colour = sig_varies_fitted, fill = Matched
  )) +
  geom_violin(draw_quantiles = c(0.5), adjust = 0.8) +
  facet_grid(vars(param), vars(levels), labeller = label_parsed) +
  geom_hline(yintercept = true_omega, lty = 1, alpha = 0.6) +
  xlab(expression(Simulated ~ omega[0] ~ variation)) +
  labs(colour = expression(Fitted ~ omega[0] ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(Estimated~omega[0])) +
  coord_cartesian(ylim = c(true_omega * 0.5, true_omega * 2)) +
  ggsidekick::theme_sleek()
#ggsave("sim2/simulations-facet-grid-violin.pdf", width = 10, height = 6)

true_omega <- exp(out$true_ln_global_omega)[1]
plot_dat <-
  out %>%
  dplyr::filter(is.na(failed_iter)) %>%
  dplyr::mutate(Matched = sig_varies_fitted == sig_varies) %>%
  group_by(param, levels, Matched, sig_varies, sig_varies_fitted) %>%
  summarise(
    lwr = quantile(ln_global_omega, 0.1),
    med = quantile(ln_global_omega, 0.50),
    upr = quantile(ln_global_omega, 0.9),
    lwr2 = quantile(ln_global_omega, 0.25),
    upr2 = quantile(ln_global_omega, 0.75)
    )

dodge_width <- 0.7
g <- plot_dat %>% ggplot(aes(x = sig_varies, y = exp(med),
  colour = sig_varies_fitted, shape = Matched)) +
  geom_hline(yintercept = true_omega, lty = 1, alpha = 0.6) +
  geom_linerange(aes(ymin = exp(lwr), ymax = exp(upr)), position = position_dodge(width = dodge_width), lwd = 0.4) +
  geom_linerange(aes(ymin = exp(lwr2), ymax = exp(upr2)), position = position_dodge(width = dodge_width), lwd = 0.8) +
  geom_point(position = position_dodge(width = dodge_width), fill = "white", size = 2) +
  facet_grid(vars(param), vars(levels), labeller = label_parsed) +
  xlab(expression(Simulated ~ omega[0] ~ variation)) +
  labs(colour = expression(Fitted ~ omega[0] ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(Estimated~omega[0])) +
  coord_cartesian(ylim = c(true_omega * 0.66, true_omega * 1.5)) +
  ggsidekick::theme_sleek() +
  scale_shape_manual(values = c(19, 21)) +
  scale_y_log10(breaks = seq(8, 36, 2))
#ggsave("sim2/simulations-facet-grid-pointrange.pdf", width = 8, height = 5)

#------------------------------------------------
#Examine bias using MRE:
#Median(E(1:100) - T(1:100) / T(1:100))

plot_dat <-
  out %>%
  dplyr::filter(is.na(failed_iter)) %>%
  dplyr::mutate(Matched = sig_varies_fitted == sig_varies) %>%
  group_by(param, levels, Matched, sig_varies, sig_varies_fitted) %>%
  mutate(MRE = 100*(exp(ln_global_omega) - true_omega ) / true_omega)

plot_dat <- plot_dat %>%
  group_by(param, levels, Matched, sig_varies, sig_varies_fitted) %>%
  summarise(
    lwr = quantile(MRE, 0.1),
    med = quantile(MRE, 0.50),
    upr = quantile(MRE, 0.9),
    lwr2 = quantile(MRE, 0.25),
    upr2 = quantile(MRE, 0.75)
    )

g <- plot_dat %>% ggplot(aes(x = sig_varies, y = med,
  colour = sig_varies_fitted, shape = Matched)) +
  geom_hline(yintercept = 0, lty = 1, alpha = 0.6) +
  geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = dodge_width), lwd = 0.4) +
  geom_linerange(aes(ymin = lwr2, ymax = upr2), position = position_dodge(width = dodge_width), lwd = 0.8) +
  geom_point(position = position_dodge(width = dodge_width), fill = "white", size = 2) +
  facet_grid(vars(param), vars(levels), labeller = label_parsed) +
  xlab(expression(Simulated ~ omega[0] ~ variation)) +
  labs(colour = expression(Fitted ~ omega[0] ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(Median ~ relative ~ error ~ '%')) +
  ggsidekick::theme_sleek() +
  scale_shape_manual(values = c(19, 21))
g

#ggsave("sim2/simulations-facet-grid-pointrange-mre.pdf", width = 8, height = 5)

#------------------------------------------------
#Examine accuracy using MARE:
#Median(abs(E(1:100) - T(1:100)) / T(1:100))

plot_dat <-
  out %>%
  dplyr::filter(is.na(failed_iter)) %>%
  dplyr::mutate(Matched = sig_varies_fitted == sig_varies) %>%
  group_by(param, levels, Matched, sig_varies, sig_varies_fitted) %>%
  mutate(MARE = 100*abs((exp(ln_global_omega) - true_omega ) / true_omega))

plot_dat <- plot_dat %>%
  group_by(param, levels, Matched, sig_varies, sig_varies_fitted) %>%
  summarise(
    lwr = quantile(MARE, 0.1),
    med = quantile(MARE, 0.50),
    upr = quantile(MARE, 0.9),
    lwr2 = quantile(MARE, 0.25),
    upr2 = quantile(MARE, 0.75)
    )

g <- plot_dat %>% ggplot(aes(x = sig_varies, y = med,
  colour = sig_varies_fitted, shape = Matched)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = dodge_width), lwd = 0.4) +
  geom_linerange(aes(ymin = lwr2, ymax = upr2), position = position_dodge(width = dodge_width), lwd = 0.8) +
  geom_point(position = position_dodge(width = dodge_width), fill = "white", size = 2) +
  facet_grid(vars(param), vars(levels), labeller = label_parsed) +
  xlab(expression(Simulated ~ omega[0] ~ variation)) +
  labs(colour = expression(Fitted ~ omega[0] ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(Median ~ absolute ~ relative ~ error ~ '%')) +
  ggsidekick::theme_sleek() +
  scale_shape_manual(values = c(19, 21))
g

#ggsave("sim2/simulations-facet-grid-pointrange-mare.pdf", width = 8, height = 5)

# Extra plotting code to plot each panel as separate .pdf:
#
# out = out %>% mutate(sim_name = case_when(
#   sim == "sim_1" ~ "rho == 0.1",
#   sim == "sim_2" ~ "rho == 0.5",
#   sim == "sim_3" ~ "rho == 0.9",
#   sim == "sim_4" ~ "kappa == 2.7",
#   sim == "sim_5" ~ "kappa == 0.57",
#   sim == "sim_6" ~ "kappa == 0.28",
#   sim == "sim_7" ~ "sigma[O] == 0.2",
#   sim == "sim_8" ~ "sigma[O] == 0.5",
#   sim == "sim_9" ~ "sigma[O] == 0.8"
# ))


# plots <- out %>%
#   dplyr::filter(is.na(failed_iter)) %>%
#   dplyr::mutate(Matched = ifelse(sig_varies_fitted == sig_varies, TRUE, FALSE)) %>%
#   group_by(sim_name) %>%
#   nest() %>%
#   mutate(plot = map2(data, sim_name, ~ ggplot(data = .x, aes(sig_varies, exp(ln_global_omega),
#     colour = sig_varies_fitted, fill = Matched
#   )) +
#     geom_boxplot() +
#     geom_hline(yintercept = exp(out[["true_ln_global_omega"]][1])) +
#     xlab(expression(Simulated ~ omega ~ variation)) +
#     labs(colour = expression(Fitted ~ omega ~ variation)) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_fill_manual(values = c("white", "grey60")) +
#     ylab(expression(omega[0])) +
#     ylim(exp(min(na.exclude(out$ln_global_omega))),
#       exp(max(na.exclude(out$ln_global_omega)))) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     annotate("text", x="by lake", y=47, label=sim_name, parse=T,
#       size=7)
#     ))
#
# plots$plot[[9]]
#
# pdf("sim2/simulations.pdf", width=7, height=5)
# plots$plot
# dev.off()
