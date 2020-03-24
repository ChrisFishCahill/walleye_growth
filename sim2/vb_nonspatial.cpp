#include <TMB.hpp>
/*
 *      Hierarchical von Bertalanffy with density dependent growth regression
 *                     Galluci and Quinn (1979) Parameterization
 *
 *                                Model structure:
 *    length_pred_i = linfinity*(1-exp(-( omega/linfinity) * (age_i (i) - t0 )))
 *
 *Chris Cahill 2019
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
  // Data
  DATA_INTEGER(Nobs);     // Observations in the dataset
  DATA_VECTOR(length_i);  // length fish i
  DATA_VECTOR(age_i);     // age of fish i
  DATA_IVECTOR(lake_i);   // grouping factor for lake
  DATA_IVECTOR(time_i);   // grouping factor for time
  DATA_INTEGER(Nlakes);   // Number of lakes

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

  // Likelihood noise term
  PARAMETER(ln_cv);

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

  for (int i = 0; i < Nobs; i++) {
    Type omega = exp(ln_global_omega +            // intercept
                     eps_omega_lake(lake_i(i)) +  // std ran eff
                     eps_omega_time(time_i(i)));  // std ran eff

    Type linf = exp(ln_global_linf +       // intercept
                    eps_linf(lake_i(i)));  // std ran eff

    Type t0 = global_tzero +      // intercept
              eps_t0(lake_i(i));  // std ran eff

    length_pred(i) = linf * (1 - exp(-(omega / linf) * (age_i(i) - t0)));

    if (!isNA(length_i(i)))
      jnll -= dlnorm(length_i(i), log(length_pred(i)), exp(ln_cv), true);
  }
  return jnll;
}
