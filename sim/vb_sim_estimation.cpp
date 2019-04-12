#include <TMB.hpp>
/*
*      Spatial-temporal von Bertalanffy with density dependent growth regression
*                     Galluci and Quinn (1979) Parameterization
*
*                                Model structure:
*    length_pred_i = linfinity*(1-exp(-( omega/linfinity) * (age_i (i) - t0 )))
*
*                                     Where:
* linfinity = linfinity_global + sex_effect_i + eps_linf
* t0  = t0_global + eps_t0
* omega   = omega_global + eta_fixed_i + eps_omega_st
*
*                       Probability of Random Coefficients:
* epsilon_linf ~ N(0, ln_sd_linf)
* epsilon_t0 ~ N(0, ln_sd_tzero)
* eps_omega_st = rho*eps_omega_s,t-1 + u_st
* where u_st ~ N(0, SIMGA) --> SIGMA is estimated as per INLA approach
*
*                                 Likelihoods:
* if(CTL == 1 == Normal)
* if(CTL == 2 == Lognormal)
* if(CTL == 3 == Gamma)
*
*Chris Cahill and Sean Anderson March 2019
*/
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;

  // Data
  DATA_INTEGER(Nobs);          //Observations in the dataset
  DATA_VECTOR(length_i);       //length fish i
  DATA_VECTOR(age_i);          //age of fish i
  DATA_IVECTOR(lake_i);        //grouping factor for lake
  DATA_MATRIX(X_ij_omega);     //covariate matrix for omega
  DATA_INTEGER(Nlakes);        //Number of lakes

  // SPDE objects
  DATA_STRUCT( spdeMatrices, spde_t );

  DATA_IVECTOR( s_i );         //Random effect index for location(s)
  DATA_IVECTOR( t_i );         //Random effect index for year(t)
  DATA_INTEGER( n_t );         //number of years
  DATA_IVECTOR( t_prev_i );    //Random effect index for year(t-1)

  DATA_INTEGER( CTL );         //Control for likelihood
  DATA_INTEGER(ar1);

  //Parameters
  PARAMETER(ln_global_omega);

  PARAMETER(ln_global_linf);
  PARAMETER(ln_sd_linf);

  PARAMETER(global_tzero);
  PARAMETER(ln_sd_tzero);
  PARAMETER_VECTOR(b_j_omega);

  //Random coefficients
  PARAMETER_ARRAY(eps_omega_st);
  PARAMETER_VECTOR(eps_linf);
  PARAMETER_VECTOR(eps_t0);

  //Likelihood noise term
  PARAMETER(ln_cv);

  //spatial-temporal terms
  PARAMETER( ln_kappa );       //matern kappa
  PARAMETER( ln_tau_O );       //spatial noise
  PARAMETER( rho_unscaled );   //autocorrelation term

  //calculate the fixed effects:
  vector<Type> eta_fixed_i = X_ij_omega * b_j_omega;

  vector<Type> length_pred(Nobs);
  vector<Type> sigma (Nobs);
  vector<Type> eps_i (Nobs);

  Type rho = Type(2)*invlogit(rho_unscaled)-Type(1);

  // Objective function
  Type jnll = 0;

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));

  // Probability of random coefficients for linf in lake l:
  for(int l=0; l<Nlakes; l++){
    jnll -= dnorm(eps_linf(l), Type(0.0), exp(ln_sd_linf), true);
  }

  // Probability of random coefficients for t0 in lake l:
  for(int l=0; l<Nlakes; l++){
    jnll -= dnorm(eps_t0(l), Type(0.0), exp(ln_sd_tzero), true);
  }

  // Probability of spatial-temporal random coefficients
  SparseMatrix<Type> Q = R_inla::Q_spde(spdeMatrices, exp(ln_kappa));

  if (ar1 == 0) {
    for(int t=0; t < n_t; t++)
      jnll += SCALE(GMRF(Q), 1.0 / exp(ln_tau_O))(eps_omega_st.col(t));
  } else {
    jnll += SCALE(SEPARABLE(AR1(rho), GMRF(Q)), 1.0 / exp(ln_tau_O))(eps_omega_st);
  }

    for(int i=0; i<Nobs; i++){
    Type omega = 0;
    Type linf  = 0;
    Type t0    = 0;
    Type sigma = 0;

    // if (t_i(i) == Type(0)){
    //   eps_i(i) = eps_omega_st( s_i(i), t_i(i) )*sqrt(1-rho*rho);
    // } else {
    //   eps_i(i) = rho * eps_omega_st((s_i(i), t_prev_i(i))) + //previous raneff at location s*rho
    //     eps_omega_st((s_i(i), t_i(i)));
    // }
    eps_i(i) = eps_omega_st( s_i(i), t_i(i) );


    omega = exp(ln_global_omega) +                  //intercept
            eta_fixed_i(i) +                        //fixed effects
            eps_i(i);                               //AR-1 space-time

    linf  = exp(ln_global_linf) +                   //intercept
            eps_linf(lake_i(i));                    //std ran eff

    t0    = global_tzero +                          //intercpet
            eps_t0(lake_i(i));                      //std ran eff

    length_pred(i) = linf*(1-exp(-(omega/linf) * (age_i (i) - t0 )));
    sigma    = exp(ln_cv)*length_pred(i);

    //CYO likelihood
    if(CTL == 1){
    //Normal
      if( !isNA(length_i(i)) ) jnll -= dnorm( length_i(i), length_pred(i), sigma, true );
    }
    if(CTL == 2){
    //Lognormal
      if( !isNA(length_i(i)) ) jnll -= dlnorm(length_i(i), log(length_pred(i)) - pow(sigma, 2)/2, sigma, true );
    }
    if(CTL == 3){
    //Gamma
      if( !isNA(length_i(i)) ) jnll -= dgamma( length_i(i), 1/pow(exp(ln_cv),2), length_pred(i)*pow(exp(ln_cv),2), true ); ;
    }
  }

  // Reporting
  REPORT(Range);
  REPORT(SigmaO);
  REPORT(rho);
  REPORT(length_pred);
  REPORT(eps_i);

  return jnll;
}
