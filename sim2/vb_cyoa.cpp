#include <TMB.hpp>
/*
 * Choose your own correlation adventure for the von Bertalanffy
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
  DATA_INTEGER(Nobs);     // Observations in the dataset
  DATA_VECTOR(length_i);  // length fish i
  DATA_VECTOR(age_i);     // age of fish i
  DATA_IVECTOR(lake_i);   // grouping factor for lake
  DATA_IVECTOR(time_i);   // grouping factor for time
  DATA_INTEGER(Nlakes);   // Number of lakes

  // SPDE objects
  DATA_STRUCT(spdeMatrices, spde_t);

  // Parameters
  PARAMETER(ln_global_linf);
  PARAMETER(ln_sd_linf);

  PARAMETER(global_tzero);
  PARAMETER(ln_sd_tzero);

  PARAMETER(ln_global_omega);
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
  PARAMETER(ln_kappa);       //matern kappa
  PARAMETER(ln_tau_O);       //spatial noise
  PARAMETER( rho_unscaled ); //autocorrelation term

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
  SparseMatrix<Type> Q = R_inla::Q_spde(spdeMatrices, exp(ln_kappa));
  jnll += SCALE(SEPARABLE(AR1(rho), GMRF(Q)), 1.0 / exp(ln_tau_O))(eps_omega_st);

  // Derived quantities
  Type Range = sqrt(8) / exp(ln_kappa);
  Type SigO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));

  for (int i = 0; i < Nobs; i++) {
    Type omega = exp(ln_global_omega +               // intercept
                     eps_omega_lake(lake_i(i)) +     // std lake ran eff
                     eps_omega_time(time_i(i)) +     // std time ran eff
                     eps_omega_st(lake_i(i), time_i(i)));  // ar1 space time

    Type linf = exp(ln_global_linf +       // intercept
                    eps_linf(lake_i(i)));  // std ran eff

    Type t0 = global_tzero +      // intercept
              eps_t0(lake_i(i));  // std ran eff

    length_pred(i) = linf * (1 - exp(-(omega / linf) * (age_i(i) - t0)));

    if (!isNA(length_i(i)))
      jnll -= dlnorm(length_i(i), log(length_pred(i)), exp(ln_cv), true);
  }
  return jnll;
  ADREPORT(Range)
  ADREPORT(SigO)
}
