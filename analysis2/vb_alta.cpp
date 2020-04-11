#include <TMB.hpp>
/*
* Cahill & Anderson 11 April 2020
*/

template <class Type>
bool isNA(Type x)
{
  return R_IsNA(asDouble(x));
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres;
  logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;

  // Data
  DATA_INTEGER(Nobs);                  // Observations in the dataset
  DATA_VECTOR(length_i);               // length fish i
  DATA_VECTOR(age_i);                  // age of fish i
  DATA_IVECTOR(lake_i);                // grouping factor for lake
  DATA_IVECTOR(time_i);                // grouping factor for time
  DATA_INTEGER(Nlakes);                // Number of lakes
  DATA_IVECTOR(sex_i);                 // factor for sex
  DATA_MATRIX(X_ij_omega);             // covariate matrix for omega
  DATA_STRUCT(spdeMatrices, spde_t);   // SPDE objects
  DATA_VECTOR(predTF_i);               // Indicator for CV fold

  DATA_SCALAR(rho_sd_prior);           // rho SD prior
  DATA_SCALAR(rho_mean_prior);         // rho mean prior
  DATA_SCALAR(tau_O_sd_prior);         // tau_O_ SD prior
  DATA_SCALAR(tau_O_mean_prior);       // tau_O_ mean prior

  // Parameters
  PARAMETER(ln_global_linf);
  PARAMETER(ln_sd_linf);

  PARAMETER(global_tzero);
  PARAMETER(ln_sd_tzero);

  PARAMETER(ln_b_sex);
  PARAMETER_VECTOR(b_j_omega);

  PARAMETER(ln_sd_omega_lake);
  PARAMETER(ln_sd_omega_time);

  // Random coefficients
  PARAMETER_VECTOR(eps_omega_lake);
  PARAMETER_VECTOR(eps_omega_time);
  PARAMETER_VECTOR(eps_linf);
  PARAMETER_VECTOR(eps_t0);
  PARAMETER_ARRAY(eps_omega_st);

  // Likelihood noise term
  PARAMETER(ln_cv);

  // Spatial-temporal terms
  PARAMETER(ln_kappa);     //matern kappa
  PARAMETER(ln_tau_O);     //spatial noise
  PARAMETER(rho_unscaled); //autocorrelation term

  // back-transform rho
  Type rho = Type(2)*invlogit(rho_unscaled)-Type(1);

  vector<Type> length_pred(Nobs);
  Type jnll = 0;

  // Probability of random coefficients:
  for (int l = 0; l < Nlakes; l++) {
    jnll -= dnorm(eps_linf(l), Type(0.0), exp(ln_sd_linf), true);
    jnll -= dnorm(eps_t0(l), Type(0.0), exp(ln_sd_tzero), true);
    jnll -= dnorm(eps_omega_lake(l), Type(0.0), exp(ln_sd_omega_lake), true);
  }
  for (int t = 0; t < eps_omega_time.size(); t++) {
     jnll -= dnorm(eps_omega_time(t), Type(0.0), exp(ln_sd_omega_time), true);
  }

  // Derived quantities
  Type Range = sqrt(8) / exp(ln_kappa);
  Type SigO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));

  SparseMatrix<Type> Q = R_inla::Q_spde(spdeMatrices, exp(ln_kappa));
  jnll += SCALE(SEPARABLE(AR1(rho), GMRF(Q)), 1.0 / exp(ln_tau_O))(eps_omega_st);

  // penalties:
  jnll -= dnorm(rho, rho_mean_prior, rho_sd_prior, true);
  jnll -= dnorm(ln_tau_O, tau_O_mean_prior, tau_O_sd_prior, true);

  // calculate fixed effects for omega:
  vector<Type> eta_fixed_i = X_ij_omega * b_j_omega;

  for (int i = 0; i < Nobs; i++) {
    Type omega = exp(eta_fixed_i(i) +      // fixed effects
      eps_omega_lake(lake_i(i)) +          // std lake ran eff
      eps_omega_time(time_i(i)) +          // std time ran eff
      eps_omega_st(lake_i(i), time_i(i))); // ar1 space time

    Type linf = exp(ln_global_linf +       // intercept
      ln_b_sex*sex_i(i) +                  // sex effect
      eps_linf(lake_i(i)));                // std ran eff

    Type t0 = global_tzero +      // intercept
      eps_t0(lake_i(i));          // std ran eff

    length_pred(i) = linf * (1 - exp(-(omega / linf) * (age_i(i) - t0)));

    if (!isNA(length_i(i)))
      jnll -= dlnorm(length_i(i), log(length_pred(i)), exp(ln_cv), true);
  }
  REPORT(Range);
  REPORT(SigO);
  REPORT(rho);
  return jnll;
}

