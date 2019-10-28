#include <TMB.hpp>
/*
*      Hierarchical von Bertalanffy with density dependent growth regression
*                     Galluci and Quinn (1979) Parameterization
*
*                                Model structure:
*    length_pred_i = linfinity*(1-exp(-( omega/linfinity) * (age_i (i) - t0 )))
*
*                                     Where:
* linfinity = linfinity_global + sex_effect_i + eps_linf
* t0  = t0_global + eps_t0
* omega   = omega_global + eta_fixed_i + eps_omega
*
*                       Probability of Random Coefficients:
* epsilon_linf ~ N(0, ln_sd_linf)
* epsilon_t0 ~ N(0, ln_sd_tzero)
* eps_omega ~ N(0, ln_sd_omega)
*
*                                 Likelihoods:
* if(CTL == 1 == Normal)
* if(CTL == 2 == Lognormal)
* if(CTL == 3 == Gamma)
*
*Chris Cahill 2019
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
  using namespace Eigen;

  // Data
  DATA_INTEGER(Nobs);          //Observations in the dataset
  DATA_VECTOR(length_i);       //length fish i
  DATA_VECTOR(age_i);          //age of fish i
  DATA_IVECTOR(lake_i);        //grouping factor for lake
  DATA_MATRIX(X_ij_omega);     //covariate matrix for omega
  DATA_VECTOR(sex_i);          //sex of fish i
  DATA_INTEGER(Nlakes);        //Number of lakes

  DATA_INTEGER( CTL );         //Control for likelihood
  DATA_VECTOR( predTF_i );     //indicator vector for fold

  //Parameters
  PARAMETER(ln_global_omega);
  PARAMETER(ln_sd_omega);

  PARAMETER(ln_global_linf);
  PARAMETER(ln_sd_linf);

  PARAMETER(global_tzero);
  PARAMETER(ln_sd_tzero);
  PARAMETER(ln_b_sex);

  PARAMETER_VECTOR(b_j_omega);

  //Random coefficients
  PARAMETER_VECTOR(eps_omega);
  PARAMETER_VECTOR(eps_linf);
  PARAMETER_VECTOR(eps_t0);

  //Likelihood noise term
  PARAMETER(ln_cv);

  //calculate the fixed effects:
  vector<Type> eta_fixed_i = X_ij_omega * b_j_omega;

  vector<Type> length_pred(Nobs);
  vector<Type> sigma (Nobs);
  vector<Type> eps_i (Nobs);

  // Objective function
  Type jnll = 0;
  Type pred_jnll = 0;

  vector<Type> jnll_i(Nobs);
  jnll_i.setZero();

  // Probability of random coefficients:
  for(int l=0; l<Nlakes; l++){
    jnll -= dnorm(eps_linf(l), Type(0.0), exp(ln_sd_linf), true);
    jnll -= dnorm(eps_t0(l), Type(0.0), exp(ln_sd_tzero), true);
    jnll -= dnorm(eps_omega(l), Type(0.0), exp(ln_sd_omega), true);
  }

  for(int i=0; i<Nobs; i++){
    Type omega = 0;
    Type linf  = 0;
    Type t0    = 0;

    omega = exp(ln_global_omega) +            //intercept
      eta_fixed_i(i) +                        //fixed effects
      eps_omega(lake_i(i));                   //std ran eff

    linf  = exp(ln_global_linf) +             //intercept
      exp(ln_b_sex)*sex_i(i) +                //sex effect
      eta_fixed_i(i) +
      eps_linf(lake_i(i));                    //std ran eff

    t0    = global_tzero +                    //intercpet
      eps_t0(lake_i(i));                      //std ran eff

    length_pred(i) = linf*(1-exp(-(omega/linf) * (age_i (i) - t0 )));

    //CYO likelihood
    if(CTL == 1){
      //Normal
      if( !isNA(length_i(i)) ) jnll_i(i) -= dnorm( length_i(i), length_pred(i), exp(ln_cv)*length_pred(i), true );
    }
    if(CTL == 2){
      //Lognormal
      if( !isNA(length_i(i)) ) jnll_i(i) -= dlnorm(length_i(i), log(length_pred(i)) - pow(exp(ln_cv), 2)/2, exp(ln_cv), true );
    }
    if(CTL == 3){
      //Gamma
      if( !isNA(length_i(i)) ) jnll_i(i) -= dgamma( length_i(i), 1/pow(exp(ln_cv),2), length_pred(i)*pow(exp(ln_cv),2), true );
    }
    // Running counter
    if( predTF_i(i)==0 ) jnll += jnll_i(i); //estimation
    if( predTF_i(i)==1 ) pred_jnll += jnll_i(i); //prediction
  }

  // Reporting
  REPORT(length_pred);
  REPORT(pred_jnll);
  REPORT(predTF_i);

  return jnll;
}
