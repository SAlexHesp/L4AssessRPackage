//  SADA ASsessment Team, February 2022
// Biomass Dynamics Model, state-space framework cpp code. 
// # create cpp code of BDM - switches for Fox and Schaefer and PT models 
// # and additions of env and dep
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <TMB.hpp>


// BDM.cpp


template<class Type>
Type x_is_gt_bound(Type x, Type bound, Type slope)
{
  // Returns 0 if x << bound and 1 if x >> bound.
  // The parameter, slope, determines the range over which the return value changes. The function
  // used is the logistic function, where the return value is 0.5 when x == bound and 0.95
  // when x95 == bound + log(19)/slope  
  Type result;
  
  result = Type(1.0) / (Type(1.0) + exp(slope * (bound - x)));
  
  return result;
}

template<class Type>
Type x_is_lt_bound(Type x, Type bound, Type slope)
{
  // Returns 0 if x >> bound and 1 if x << bound.
  // The parameter, slope, determines the range over which the return value changes. The function
  // used is the logistic function, where the return value is 0.5 when x == bound and 0.05
  // when x05 == bound + log(19)/slope.
  Type result;
  
  result = Type(1.0) - (Type(1.0) / (Type(1.0) + exp(slope * (bound - x))));
  
  return result;
}

template<class Type>
Type x_is_gt_lower_bound_and_lt_upper_bound(Type x, Type lower_bound, Type upper_bound, Type slope)
{
  // Returns 0 if x << lower_bound or >> upper_bound and 1 if x >> lower_bound and x << upper_bound.
  // The parameter, slope, determines the range over The function used is the logistic function,
  // where the return value is 0.5 when x == lower_bound or x == upper_bound, 
  // 0.95 when x95 == lower_bound + log(19)/slope, and 0.05 when x05 == upper_bound + log(19)/slope. 
  Type result;
  
  result = (Type(1.0) / (Type(1.0) + exp(slope * (lower_bound - x)))) * 
    (Type(1.0) - (Type(1.0) / (Type(1.0) + exp(slope * (upper_bound - x)))));
  
  return result;
}


template<class Type>
Type posfun(Type x, Type eps, Type pen)
{
  // The argument, pen, is not used but included for consistency with original posfun function from ADMB 
  Type bound;
  Type result;
  
  // Assume eps = log(19)/slope
  Type slope = log(Type(19.0)) / eps;
  
  bound = eps;  
  
  result = x *  x_is_gt_bound(x, bound, slope);
  
  return result;
}


template<class Type>
Type penfun(Type x, Type eps, Type pen)
{
  // The argument, pen, is not used but included for consistency with the posfun function above.
  Type bound = Type(0.0);
  Type result;
  
  // Assume eps = log(19)/slope
  Type slope = log(Type(19.0)) / eps;
  
  bound = eps;
  
  result = (x-eps) * (x-eps) * x_is_lt_bound(x, bound, slope);
  
  return result;
}


template <class Type> Type square(Type x){return x*x;}


template<class Type> Type objective_function<Type>::operator() () {
  // Biomass Dynamics Model
  // The model may be fitted without process error (scenario 1), in which case all error is
  // assumed to be contributed by observation error of the CPUE indices, and
  // the observed environmental indices is assumed to be an explanatory variable. 
  // The model may also be fitted with process error (scenario 2), in which case the annual
  // biomass is calculated with  the error associated with a random effect, where, in
  // models with environmental effects, the observed environmental variable is assumed 
  // to be an index of the annual random effect, i.e. environment (temp/chl) is assumed to 
  // have a relationship with fish production.
  
  // Declaration of data to be input from R
  // ======================================
  
  DATA_INTEGER(mod_type);       // model type - 1=Schaefer, 2=Fox, 3=Pella-Tomlinson
  DATA_INTEGER(mod_option);     // model option - 1=simple bdm, 2=bdm + env, 3=bdm + depensation, 4=bdm + both
  DATA_INTEGER(mod_scenario);   // model scenario - 1=without process error, 2=with process error
  DATA_INTEGER(nyrs);
  DATA_INTEGER(cpue1_nyrs);
  DATA_INTEGER(cpue2_nyrs);
  DATA_INTEGER(cpue3_nyrs);
  DATA_INTEGER(cpue4_nyrs);
  DATA_INTEGER(first_cpue1_yr);
  DATA_INTEGER(last_cpue1_yr);
  DATA_INTEGER(first_cpue2_yr);
  DATA_INTEGER(last_cpue2_yr);
  DATA_INTEGER(first_cpue3_yr);
  DATA_INTEGER(last_cpue3_yr);
  DATA_INTEGER(first_cpue4_yr);
  DATA_INTEGER(last_cpue4_yr);
  DATA_VECTOR(season);
  DATA_VECTOR(tot_catch);           // Total catch 
  DATA_VECTOR(obs_ln_cpue1);             // Natural log of observed CPUE 1
  DATA_VECTOR(obs_ln_cpue2);           // Natural log of observed CPUE 2
  DATA_VECTOR(obs_ln_cpue3);            // Natural log of observed CPUE 3
  DATA_VECTOR(obs_ln_cpue4);            // Natural log of observed CPUE 4
  DATA_VECTOR(cpue1_se);               // SE of observed values of the natural logs of cpue1, obs_ln_cpue1
  DATA_VECTOR(cpue2_se);             // SE of observed values of the natural logs of cpue2, obs_ln_cpue2
  DATA_VECTOR(cpue3_se);              // SE of observed values of the natural logs of cpue3, obs_ln_cpue3
  DATA_VECTOR(cpue4_se);              // SE of observed values of the natural logs of cpue4, obs_ln_cpue4
  DATA_VECTOR(obs_env);             // Observed environmental effect (chl or temp)
  DATA_VECTOR(env_se);         // SE of observed enviro values
  DATA_SCALAR(mean_obs_env); // read in mean enviro value from dat file
  DATA_SCALAR(SigmaR); // if specifying sigma R as single value and not estimating
  DATA_SCALAR(Sigma_env); // if specifying sigma env as single value, and not estimating
  DATA_SCALAR(max_currBrel); // for setting a penalty for current level of stock depletion
  
  
  // The weighting factors by which the penalties are multiplied in order to constrain the
  // parameter estimates to feasible ranges
  DATA_SCALAR(wt_param_pen);        // parameter penalty
  DATA_SCALAR(wt_depl_pen);         // depletion penalty
  DATA_SCALAR(wt_biom_pen);         // biomass penalty
  DATA_SCALAR(wt_harv_pen);         // harvest penalty
  DATA_SCALAR(wt_cpue1);  
  DATA_SCALAR(wt_cpue2);
  DATA_SCALAR(wt_cpue3);
  DATA_SCALAR(wt_cpue4);
  
  // Bounds for the parameters
  DATA_VECTOR(uppbound);
  DATA_VECTOR(lowbound);
  
  // Declaration of parameters to be input from R
  // ============================================
  
  // Parameters for the Schaefer or Fox versions of the Biomass Dynamics Models, both of which are > 0
  PARAMETER(pInit);  // Proportion unfished, for initial population
  PARAMETER(ln_K);
  PARAMETER(ln_r);
  PARAMETER(ln_pt_parm);   // Pella-Tomlinson shape parameter
  
  // Parameters for the catchabilities (> 0) where estimates of CPUE = q * Biomass 
  // errors in log space are assumed to be normally distributed.  
  PARAMETER(ln_q1);
  PARAMETER(ln_q2);
  PARAMETER(ln_q3);
  PARAMETER(ln_q4);
  
  // Parameter representing the standard deviation (> 0) of the deviations of the natural logs of the CPUES
  // from the natural logs of the values of the predicted CPUEs
  PARAMETER(ln_sd1);
  PARAMETER(ln_sd2);
  PARAMETER(ln_sd3);
  PARAMETER(ln_sd4);
  PARAMETER(env_p); // Parameter determining the effect of the environment 
  PARAMETER(lndep); // Parameter determining the effect of depensation 
  PARAMETER_VECTOR(FF);      // Estimated exploitation parameters.
  PARAMETER_VECTOR(EpsR);   // // Parameters for state-space random effect - process error (size == nyrs)
  
  // Declaration and initializaion of penalty terms
  Type param_pen = Type(0.0);
  Type K_pen = Type(0.0);
  Type r_pen = Type(0.0);
  Type q1_pen = Type(0.0);
  Type q2_pen = Type(0.0);
  Type q3_pen = Type(0.0);
  Type q4_pen = Type(0.0);
  Type sd1 = Type(0.0);
  Type sd2 = Type(0.0);
  Type sd3 = Type(0.0);
  Type sd4 = Type(0.0);
  Type sd1_pen = Type(0.0);
  Type sd2_pen = Type(0.0);
  Type sd3_pen = Type(0.0);
  Type sd4_pen = Type(0.0);
  Type env_pen = Type(0.0);
  Type dep_pen = Type(0.0);
  Type pt_pen = Type(0.0);
  Type depl_pen = Type(0.0);
  Type harv_pen = Type(0.0);
  Type biom_pen = Type(0.0);
  Type pInit_pen = Type(0.0);  
  
  
  // Calculate penalties if the parameters lie outside specified feasible ranges, and reset
  // to feasible ranges for population being modeled
  Type eps = Type(0.001);
  Type temp_parm;
  
  K_pen = penfun(ln_K - lowbound(0), eps, K_pen);
  temp_parm = lowbound(0) + posfun(ln_K - lowbound(0), eps, K_pen);
  K_pen += penfun(uppbound(0) - temp_parm, eps, K_pen);
  temp_parm = uppbound(0) - posfun(uppbound(0) - temp_parm, eps, K_pen);
  Type K = exp(temp_parm);
  param_pen += K_pen;
  
  r_pen = penfun(ln_r - lowbound(1), eps, r_pen);
  temp_parm = lowbound(1) + posfun(ln_r - lowbound(1), eps, r_pen);
  r_pen += penfun(uppbound(1) - temp_parm, eps, r_pen);
  temp_parm = uppbound(1) - posfun(uppbound(1) - temp_parm, eps, r_pen);
  Type r = exp(temp_parm);
  param_pen += r_pen;
  
  q1_pen = penfun(ln_q1 - lowbound(2), eps, q1_pen);
  temp_parm = lowbound(2) + posfun(ln_q1 - lowbound(2), eps, q1_pen);
  q1_pen += penfun(uppbound(2) - temp_parm, eps, q1_pen);
  temp_parm = uppbound(2) - posfun(uppbound(2) - temp_parm, eps, q1_pen);  
  Type q1 = exp(temp_parm);
  param_pen += q1_pen;
  
  q2_pen = penfun(ln_q2 - lowbound(3), eps, q2_pen);
  temp_parm = lowbound(3) + posfun(ln_q2 - lowbound(3), eps, q2_pen);
  q2_pen += penfun(uppbound(3) - temp_parm, eps, q2_pen);
  temp_parm = uppbound(3) - posfun(uppbound(3) - temp_parm, eps, q2_pen);
  Type q2 = exp(temp_parm);
  param_pen += q2_pen;
  
  q3_pen = penfun(ln_q3 - lowbound(4), eps, q3_pen);
  temp_parm = lowbound(4) + posfun(ln_q3 - lowbound(4), eps, q3_pen);
  q3_pen += penfun(uppbound(4) - temp_parm, eps, q3_pen);
  temp_parm = uppbound(4) - posfun(uppbound(4) - temp_parm, eps, q3_pen);
  Type q3 = exp(temp_parm);
  param_pen += q3_pen;
  
  q4_pen = penfun(ln_q3 - lowbound(5), eps, q4_pen);
  temp_parm = lowbound(5) + posfun(ln_q3 - lowbound(5), eps, q4_pen);
  q4_pen += penfun(uppbound(5) - temp_parm, eps, q4_pen);
  temp_parm = uppbound(5) - posfun(uppbound(5) - temp_parm, eps, q4_pen);
  Type q4 = exp(temp_parm);
  param_pen += q4_pen;

  if (cpue1_nyrs>0) {  
    sd1_pen = penfun(ln_sd1 - lowbound(6), eps, sd1_pen);
    temp_parm = lowbound(6) + posfun(ln_sd1 - lowbound(6), eps, sd1_pen);
    sd1_pen += penfun(uppbound(6) - temp_parm, eps, sd1_pen);
    temp_parm = uppbound(6) - posfun(uppbound(6) - temp_parm, eps, sd1_pen);
    sd1 = exp(temp_parm);
    param_pen += sd1_pen;
  }
  
  if (cpue2_nyrs>0) {
    sd2_pen = penfun(ln_sd2 - lowbound(7), eps, sd2_pen);
    temp_parm = lowbound(7) + posfun(ln_sd2 - lowbound(7), eps, sd2_pen);
    sd2_pen += penfun(uppbound(7) - temp_parm, eps, sd2_pen);
    temp_parm = uppbound(7) - posfun(uppbound(7) - temp_parm, eps, sd2_pen);
    sd2 = exp(temp_parm);
    param_pen += sd2_pen;    
  }

  if (cpue3_nyrs>0) {
    sd3_pen = penfun(ln_sd3 - lowbound(8), eps, sd3_pen);
    temp_parm = lowbound(8) + posfun(ln_sd3 - lowbound(8), eps, sd3_pen);
    sd3_pen += penfun(uppbound(8) - temp_parm, eps, sd3_pen);
    temp_parm = uppbound(8) - posfun(uppbound(8) - temp_parm, eps, sd3_pen);
    sd3 = exp(temp_parm);
    param_pen += sd3_pen;
  }
  
  if (cpue4_nyrs>0) {
    sd4_pen = penfun(ln_sd4 - lowbound(9), eps, sd4_pen);
    temp_parm = lowbound(9) + posfun(ln_sd4 - lowbound(9), eps, sd4_pen);
    sd4_pen += penfun(uppbound(9) - temp_parm, eps, sd4_pen);
    temp_parm = uppbound(9) - posfun(uppbound(9) - temp_parm, eps, sd4_pen);
    sd4 = exp(temp_parm);
    param_pen += sd4_pen;
  }
  
  env_pen = penfun(env_p - lowbound(10), eps, env_pen);
  temp_parm = lowbound(10) + posfun(env_p - lowbound(10), eps, env_pen);
  env_pen += penfun(uppbound(10) - temp_parm, eps, env_pen);
  temp_parm= uppbound(10) - posfun(uppbound(10) - temp_parm, eps, env_pen);
  Type env_param = temp_parm;
  param_pen += env_pen;
  
  dep_pen = penfun(lndep - lowbound(11), eps, dep_pen);
  temp_parm = lowbound(11) + posfun(lndep - lowbound(11), eps, dep_pen);
  dep_pen += penfun(uppbound(11) - temp_parm, eps, dep_pen);
  temp_parm = uppbound(11) - posfun(uppbound(11) - temp_parm, eps, dep_pen);
  Type dep = exp(temp_parm);
  param_pen += dep_pen;
  
  pt_pen = penfun(ln_pt_parm - lowbound(12), eps, pt_pen);
  temp_parm = lowbound(12) + posfun(ln_pt_parm - lowbound(12), eps, pt_pen);
  pt_pen += penfun(uppbound(12) - temp_parm, eps, pt_pen);
  temp_parm = uppbound(12) - posfun(uppbound(12) - temp_parm, eps, pt_pen);
  Type pt_parm = exp(temp_parm);
  param_pen += pt_pen;
  
  pInit_pen = penfun(pInit - lowbound(13), eps, pInit_pen);
  temp_parm = lowbound(13) + posfun(pInit - lowbound(13), eps, pInit_pen);
  pInit_pen += penfun(uppbound(13) - temp_parm, eps, pInit_pen);
  temp_parm = uppbound(13) - posfun(uppbound(13) - temp_parm, eps, pInit_pen);
  Type pInit_parm = temp_parm;
  param_pen += pInit_pen;
  
  std::cout << "pInit " << pInit << " lbnd " << lowbound(13) << " ubnd " << uppbound(13) <<
   " pInit_parm " << pInit_parm << std::endl;
  
  // Declare variables in which to store outputs
  vector<Type> biom(nyrs+1); // biomass
  vector<Type> prodn(nyrs+1); // production
  vector<Type> relbiom(nyrs); // relative biomass
  vector<Type> harv(nyrs); // harvest rate
  vector<Type> Fmort(nyrs); //fishing mortality
  vector<Type> est_ln_cpue1(nyrs); //estimated/est cpue 1
  vector<Type> est_ln_cpue2(nyrs); //estimated cpue 2
  vector<Type> est_ln_cpue3(nyrs); //estimated cpue 3
  vector<Type> est_ln_cpue4(nyrs); //estimated cpue 4
  vector<Type> est_env(nyrs);  
  
  // when est depensation
  vector<Type>prod_without_dep(nyrs+1);
  vector<Type>depensation_factor(nyrs+1);
  
  // cpue time series
  vector<Type> sq1(nyrs); // SS  - could probably use single SS for each one and zero before running each
  vector<Type> sq2(nyrs); //
  vector<Type> sq3(nyrs); //
  vector<Type> sq4(nyrs); //
  vector<Type> lnresid1(nyrs); // residuals
  vector<Type> lnresid2(nyrs); //
  vector<Type> lnresid3(nyrs); //
  vector<Type> lnresid4(nyrs); //
  vector<Type> var_cpue1(nyrs); // variance
  vector<Type> var_cpue2(nyrs);
  vector<Type> var_cpue3(nyrs);
  vector<Type> var_cpue4(nyrs);
  
  //  when including envt
  vector<Type> env_var(nyrs); 
  vector<Type> resid_env(nyrs); 
  vector<Type> sq_resid_env(nyrs); 
  vector<Type> Expl(nyrs); // exploitation rate
  vector<Type> Chat(nyrs); // expected catch 
  
  Type pi = 4*atan(1); // pi
  Type tempBiom; // temporary biomass variable
  Type NLL_CPUE1;
  Type NLL_CPUE2;
  Type NLL_CPUE3;
  Type NLL_CPUE4;
  Type NLL_env;
  Type EpsR_pen;
  Type NLL_Chat;
  Type final_depl;
  Type NLL;  
    
    

  // Calculate the predicted values of environment in year i for year from 0 to nyrs - 1
  if (mod_scenario == 1) {
    // For scenario 1, the predicted value of environment is just the mean as we assume
    // no process error
    for(int i=0; i<nyrs; i++) {
      est_env(i) = mean_obs_env;
    }
  } else {
    
    // Scenario 2
    // std::cout << "mod_scenario: " << mod_scenario << "  mod_option: " << mod_option << std::endl;
    if ((mod_option == 1) || (mod_option == 3)) {
      // For scenario 2 and model types 1 and 3 the predicted values of environment are just the mean
      // as no relationship with process errors are assumed. 
      // We should set the parameter to zero and use the map to avoid fitting it.
      for(int i=0; i<nyrs; i++) {
        est_env(i) = mean_obs_env;
      }
    } else {
      for(int i=0; i<nyrs; i++) {
        // For scenario 2 and model types 2 and 4 the predicted values of environment are calculated from the
        // random process errors.
        est_env(i) = mean_obs_env + env_param * EpsR(i);
      }
    }
    // } // mod_scenario==2 loop
  }


  // Calculate the expected biomass at the start of each year
  
  // Set the biomass at the start of the first year, i.e.
  // year 0, to the carrying capacity
  biom(0) = K * exp(env_param * obs_env(0)) * pInit_parm;
  relbiom(0) = pInit_parm;
  
  // Loop over the remaining years, updating the biomass for each year using the
  // biomass dynamics model and the effect of process error (if any)
  for(int i=0; i<nyrs; i++){
    Expl(i) = Type(1.0) / (Type(1.0) + exp(-FF(i))); //exploitation rate 1-1+exp(F) - estimated F
    
    // Calculate the depletion, harvest fraction, and fishing mortality
    if(i>0) {
      relbiom(i) = biom(i) / (K * exp(env_param * obs_env(i)));
      harv(i) = tot_catch(i) / biom(i);
      harv_pen += penfun(Type(1.0) - harv(i), eps, harv_pen);
      harv(i)  = Type(1.0) - posfun(Type(1.0) - harv(i), Type(0.0001), biom_pen);
      if (harv(i) <= Type(0.0001)) {
        std::cout << "harv(i) <= eps " << harv(i) << std::endl;
      }
      Fmort(i) = -log(Type(1.0)-harv(i));
    }    
    
    //////////// PRODUCTION FUNCTION for non-state space model ////////////
    if (mod_scenario == 1) { // No process error
      
      //  SCHAEFER production function
      if(mod_type==1){ // Schaefer production function
        if(mod_option==1){ // standard schaefer model
          prodn(i) = (r * biom(i)) * (Type(1.0) - (biom(i)/K));
        }
        
        if(mod_option==2) { // environment
          prodn(i) = (r*biom(i)) * (Type(1.0)-(biom(i)/K)) * exp(env_param * (obs_env(i) - mean_obs_env)); 
        }
        
        if(mod_option==3){ // depensation
          prod_without_dep(i) = (r * biom(i)) * (Type(1.0) - (biom(i) / K)); 
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
        if(mod_option==4){ // depensation and environment
          prod_without_dep(i) = (r * biom(i)) * (Type(1.0) - (biom(i) / K)) * exp(env_param * (obs_env(i) - mean_obs_env));
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
      } // Schaefer production function
      
      // FOX production function
      if(mod_type==2){ // Fox production function
        if(mod_option==1){ // standard Fox model
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }          
          prodn(i) = log(K) * r * biom(i) * (Type(1.0) -(log(biom(i)) / log(K)));  
        }
        
        if(mod_option==2){  // plus envt
          prodn(i) = log(K) * r * biom(i) * (Type(1.0) -(log(biom(i)) / log(K))) * exp(env_param * (obs_env(i) - mean_obs_env));
        }
        
        if(mod_option==3){ // depensation
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }   
          prod_without_dep(i) = log(K) * r * biom(i) * (Type(1.0) - (log(biom(i)) / log(K))); 
          depensation_factor(i) = Type(1.0) - exp((log(0.5) * biom(i)) / (dep * K)); 
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
        if(mod_option==4){ // depensation plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }             
          prod_without_dep(i) = log(K) * r * biom(i) * (Type(1.0) -(log(biom(i)) / log(K))) * exp(env_param * (obs_env(i) - mean_obs_env));  
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
      } // Fox production function
      
      // PELLA-TOMLINSON production function
      if(mod_type==3){ // Pella-Tomlinson production function
        if(mod_option==1){ // standard PT model
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prodn(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
        }
        
        if(mod_option==2){  // plus envt
          prodn(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm)) * exp(env_param * (obs_env(i) - mean_obs_env));
        }
        
        if(mod_option==3){ // depensation
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prod_without_dep(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
          depensation_factor(i) = Type(1.0) - exp((log(0.5) * biom(i)) / (dep * K));
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
        if(mod_option==4){ // depensation plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prod_without_dep(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm)) * exp(env_param * (obs_env(i) - mean_obs_env));
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
      } // Pella-Tomlinson production function
      
      
    } else {
      // Process error is to be factored into the calculation.
      
      // SCHAEFER production function
      if(mod_type==1){ // Schaefer production function
        
        if(mod_option==1){ // standard schaefer model
          prodn(i) = (r * biom(i)) * (Type(1.0) - (biom(i)/K));
        }
        
        if(mod_option==2) { // envt
          prodn(i) = (r * biom(i)) * (Type(1.0) - (biom(i)/K)); 
        }
        
        if(mod_option==3){ // depensation
          prod_without_dep(i) = (r * biom(i)) * (Type(1.0) - (biom(i) / K)); 
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        
        if(mod_option==4){ // depensation and envt
          prod_without_dep(i) = (r * biom(i)) * (Type(1.0) - (biom(i) / K));
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
      } // Schaefer production function
      
      // FOX production function
      if(mod_type==2){ // Fox production function
        if(mod_option==1){ // standard Fox model
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }   
          prodn(i) = log(K) * r * biom(i) * (Type(1.0) - (log(biom(i)) / log(K)));  
        }
        if(mod_option==2){  // plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }   
          prodn(i) = log(K) * r * biom(i) * (Type(1.0) - (log(biom(i)) / log(K)));
        }
        if(mod_option==3){ // depensation
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }   
          prod_without_dep(i) = log(K) * r * biom(i) * (Type(1.0) -(log(biom(i)) / log(K))); 
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); 
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        if(mod_option==4){ // depensation plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }   
          prod_without_dep(i) = log(K) * r * biom(i) * (Type(1.0) -(log(biom(i)) / log(K)));  
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
      } //  Fox production function      
      
      // PELLA-TOMLINSON production function
      if(mod_type==3){ // Pella-Tomlinson production function
        if(mod_option==1){ // standard PT model
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prodn(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
        }
        if(mod_option==2){  // plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prodn(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
        }
        if(mod_option==3){ // depensation
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prod_without_dep(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K));
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
        if(mod_option==4){ // depensation plus envt
          if (K <= eps) {
            std::cout << "K <= eps: " << K << std::endl;
          }
          if (biom(i) <= eps) {
            std::cout << "biom(i) <= eps: " << biom(i) << std::endl;
          }
          prod_without_dep(i) = (r/pt_parm) * biom(i) * (Type(1.0) - pow((biom(i) / K),pt_parm));
          depensation_factor(i) = Type(1.0) - exp((log(Type(0.5)) * biom(i)) / (dep * K)); // annual depensation
          prodn(i) = depensation_factor(i) * prod_without_dep(i);
        }
      } // Pella-Tomlinson production function
    } // outer loop with annual exploitation rate
    
    
    // Calculate the expected catch - using the exploitation rate and biomass calc above 
    Chat(i) = Expl(i) * biom(i);
    
    // Calculate the biomass at the start of the next year, allowing for production and catch and, if present,
    // process error.
    if (mod_scenario == 1) {
      
      biom(i+1) = biom(i) + prodn(i) - Chat(i); 
    }
    if (mod_scenario == 2) {
      // Norms version
      // biom(i+1) = biom(i) + prodn(i) * exp(EpsR(i)) - Expl(i) * biom(i);
      
      // Andre's version 
      biom(i+1) = (biom(i) + prodn(i) - Chat(i)) * exp(EpsR(i));      
    }    
    
    // Ensure that the biomass does not fall below zero
    tempBiom = biom(i+1);
    biom_pen += penfun(tempBiom - Type(1.0), Type(1.0), biom_pen);
    biom(i+1) = Type(1.0) + posfun(tempBiom - Type(1.0), Type(1.0), biom_pen);
    
  }   // End of loop over years
  
  // depletion penalty and final depletion penalty
  depl_pen = penfun(max_currBrel - relbiom(nyrs-1), Type(0.01), depl_pen);
  final_depl = max_currBrel - posfun(max_currBrel - relbiom(nyrs-1), Type(0.01), depl_pen);
  
  
  // NLL calcs
  NLL_CPUE1 = Type(0.0);
  NLL_CPUE2 = Type(0.0);
  NLL_CPUE3 = Type(0.0);
  NLL_CPUE4 = Type(0.0);
  NLL_env = Type(0.0);
  EpsR_pen = Type(0.0);
  NLL_Chat = Type(0.0);
  NLL=Type(0.0);

  // Calculate the contribution to the variance that is due to the deviation of the expected value
  // of each CPUE from the expected value predicted by the model
  Type var1 = sd1 * sd1;
  Type var2 = sd2 * sd2;
  Type var3 = sd3 * sd3;
  Type var4 = sd4 * sd4;
  Type sum_ratio1 = Type(0.0);
  Type sum_ratio2 = Type(0.0);
  Type sum_ratio3 = Type(0.0);
  Type sum_ratio4 = Type(0.0);
  
  for(int i=0; i<nyrs; i++){
    
    // Calculate expected CPUE for each series
    est_ln_cpue1(i) = log(q1 * biom(i));
    est_ln_cpue2(i) = log(q2 * biom(i));
    est_ln_cpue3(i) = log(q3 * biom(i));
    est_ln_cpue4(i) = log(q4 * biom(i));
    
    // Calculate the NLL for CPUEs,
    if (cpue1_nyrs>0) {
      if((season(i) >= first_cpue1_yr) & (season(i) <= last_cpue1_yr)){ //specifying years for cpue 1
        
        if(exp(obs_ln_cpue1(i)) > Type(0.0)){
          lnresid1(i) = obs_ln_cpue1(i) - est_ln_cpue1(i); // observed - exp
          sq1(i) = lnresid1(i) * lnresid1(i); // Sum squares
          var_cpue1(i) = cpue1_se(i) * cpue1_se(i); // Calculate the contribution to the variance that is due to sampling 
          
          // The distribution of observation errors of the natural logs of the CPUES is assumed to be
          // normal, with mean of zero and variance equal to the combination of sampling error and
          // deviation of the expected value of the log of CPUE from the log of the value predicted
          // by the model.
          NLL_CPUE1 += Type(0.5) * log(var1 + var_cpue1(i)) + Type(0.5) * log(Type(2.0)*pi) +
            (sq1(i)/(Type(2.0) * (var1 + var_cpue1(i))));
          
        } // if obs_cpue >0 loop
      } // cpue1 loop
    }
    
    // Calculate the NLL for CPUE 2
    if (cpue2_nyrs>0) {
      if((season(i)>=first_cpue2_yr) & (season(i)<= last_cpue2_yr)){
        
        lnresid2(i) = obs_ln_cpue2(i)- est_ln_cpue2(i);
        sq2(i) = lnresid2(i) * lnresid2(i);
        var_cpue2(i) = cpue2_se(i) * cpue2_se(i);
        
        NLL_CPUE2 += Type(0.5) * log(var2 + var_cpue2(i)) + Type(0.5) * log(Type(2.0)*pi) +
          sq2(i) /(Type(2.0) * (var2 + var_cpue2(i)));
        
      } // cpue 2 loop
    }
    
    // Calculate the NLL for the CPUE 3
    if (cpue3_nyrs>0) {
      if((season(i)>=first_cpue3_yr) & (season(i)<= last_cpue3_yr)){
        
        lnresid3(i) = obs_ln_cpue3(i) - est_ln_cpue3(i);
        sq3(i) = lnresid3(i) * lnresid3(i);
        var_cpue3(i) = cpue3_se(i) * cpue3_se(i);
        
        NLL_CPUE3 += Type(0.5) * log(var3 + var_cpue3(i))+ Type(0.5) * log(Type(2.0)*pi) +
          sq3(i) /(Type(2.0) *(var3 + var_cpue3(i)));
        
      } //loop for cpue 3
    }
    
    // Calculate the NLL for the CPUE 4
    if (cpue4_nyrs>0) {
      if((season(i)>=first_cpue4_yr) & (season(i)<= last_cpue4_yr)){
        
        lnresid4(i) = obs_ln_cpue4(i) - est_ln_cpue4(i);
        sq4(i) = lnresid4(i) * lnresid4(i);
        var_cpue4(i) = cpue4_se(i) * cpue4_se(i);
        
        NLL_CPUE4 += Type(0.5) * log(var4 + var_cpue4(i))+ Type(0.5) * log(Type(2.0)*pi) +
          sq4(i) /(Type(2.0) *(var4 + var_cpue4(i)));
        
        
      } //loop for cpue 4
    }
    
    // calculating NLL associated with environmental data
    resid_env(i) = obs_env(i) - est_env(i);
    sq_resid_env(i) = resid_env(i) * resid_env(i);
    env_var(i) = env_se(i) * env_se(i);
    NLL_env += Type(0.5) * log((Sigma_env * Sigma_env) + env_var(i)) + Type(0.5) * log(Type(2.0) * pi) +
      sq_resid_env(i) / (Type(2.0) * (Sigma_env * Sigma_env) + env_var(i));
    
  } // yr loop
  
  //  apply penalties
  param_pen = param_pen * wt_param_pen;
  depl_pen = depl_pen * wt_depl_pen;
  biom_pen = biom_pen * wt_biom_pen;
  harv_pen = harv_pen * wt_harv_pen;
  
  NLL_CPUE1 = NLL_CPUE1 * wt_cpue1;
  NLL_CPUE2 = NLL_CPUE2 * wt_cpue2;
  NLL_CPUE3 = NLL_CPUE3 * wt_cpue3;
  NLL_CPUE4 = NLL_CPUE4 * wt_cpue4;
  
  EpsR_pen = - sum(dnorm(EpsR, Type(0.0), SigmaR, true));
  
  NLL_Chat = - sum(dnorm(log(tot_catch), log(Chat), Type(0.025), true)); // catch NLL
  
  //  OBJECTIVE FUNCTION
  NLL = NLL_CPUE1 + NLL_env + NLL_Chat + NLL_CPUE2 + NLL_CPUE3 + NLL_CPUE4 + 
    biom_pen + depl_pen + param_pen + harv_pen + EpsR_pen;
  
  std::cout << "NLL: " << NLL << "NLL_CPUE1: " << NLL_CPUE1 << "  NLL_CPUE2: " << NLL_CPUE2 << "  NLL_CPUE3: " << NLL_CPUE3 <<  "  NLL_CPUE4: " << NLL_CPUE4 << 
    "  NLL_env: " << NLL_env << " NLL_Chat " << NLL_Chat << std::endl;
  std::cout << "biom_pen: " << biom_pen << "  depl_pen: " << depl_pen << "  param_pen: " << param_pen << "  harv_pen: " << 
    harv_pen << "  EpsR_pen: " << EpsR_pen << std::endl;
  std::cout << "" << std::endl;
  

  REPORT(ln_K);
  REPORT(ln_r);
  REPORT(ln_q1);
  if (cpue2_nyrs>0) REPORT(ln_q2);
  if (cpue3_nyrs>0) REPORT(ln_q3);
  if (cpue4_nyrs>0) REPORT(ln_q4);
  REPORT(ln_sd1);
  if (cpue2_nyrs>0) REPORT(ln_sd2);
  if (cpue3_nyrs>0) REPORT(ln_sd3);
  if (cpue4_nyrs>0) REPORT(ln_sd4);
  REPORT(env_param);
  REPORT(lndep);  
  REPORT(dep);
  REPORT(ln_pt_parm);
  REPORT(FF);
  REPORT(EpsR);
  
  // Back-transformed parameters in arithmetic space
  REPORT(pInit);
  REPORT(K);
  REPORT(r);
  REPORT(q1);
  if (cpue2_nyrs>0) REPORT(q2);
  if (cpue3_nyrs>0) REPORT(q3);
  if (cpue4_nyrs>0) REPORT(q4);
  REPORT(sd1);
  if (cpue2_nyrs>0) REPORT(sd2);
  if (cpue3_nyrs>0) REPORT(sd3);
  if (cpue4_nyrs>0) REPORT(sd4);
  REPORT(env_param);
  if (mod_option==3 || mod_option==4) REPORT(dep);
  REPORT(pt_parm);
  REPORT(SigmaR);
  REPORT(Sigma_env);
  REPORT(Expl);
  REPORT(Chat);
  
  // Data used in the model
  REPORT(season);
  REPORT(tot_catch);
  REPORT(obs_ln_cpue1);
  REPORT(obs_ln_cpue2);
  REPORT(obs_ln_cpue3);
  REPORT(obs_ln_cpue4);
  REPORT(obs_env);
  REPORT(cpue1_se);
  REPORT(cpue2_se);
  REPORT(cpue3_se);
  REPORT(cpue4_se);
  REPORT(env_se);
  
  // Derived variables
  REPORT(biom);
  REPORT(relbiom);
  REPORT(harv);
  REPORT(Fmort);
  REPORT(prodn);
  REPORT(prod_without_dep);
  REPORT(depensation_factor);
  REPORT(est_ln_cpue1);
  REPORT(est_ln_cpue2);
  REPORT(est_ln_cpue3);
  REPORT(est_ln_cpue4);
  REPORT(est_env);
  REPORT(final_depl);
  
  // Contributions to NLL
  REPORT(NLL);
  REPORT(NLL_CPUE1);
  REPORT(NLL_CPUE2);
  REPORT(NLL_CPUE3);
  REPORT(NLL_CPUE4);
  REPORT(NLL_env);
  REPORT(NLL_Chat);
  
  // Penalties
  REPORT(K_pen);
  REPORT(r_pen);
  REPORT(q1_pen);
  REPORT(q2_pen);
  REPORT(q3_pen);
  REPORT(q4_pen);
  REPORT(sd1_pen);
  REPORT(sd2_pen);
  REPORT(sd3_pen);
  REPORT(sd4_pen);
  REPORT(env_pen);
  REPORT(dep_pen);
  REPORT(param_pen);
  REPORT(biom_pen);
  REPORT(depl_pen);
  REPORT(harv_pen);
  REPORT(EpsR_pen);
  REPORT(uppbound);
  REPORT(lowbound);
   
  // asymptotic standard errors 
  ADREPORT(K);
  ADREPORT(r);
  if (cpue1_nyrs>0) ADREPORT(q1);
  if (cpue2_nyrs>0) ADREPORT(q2);
  if (cpue3_nyrs>0) ADREPORT(q3);
  if (cpue4_nyrs>0) ADREPORT(q4);
  if (mod_option==3 || mod_option==4) ADREPORT(dep);
  if (mod_option==2 || mod_option==4) ADREPORT(env_param);
  
  // Derived variables
  ADREPORT(Chat);
  ADREPORT(biom);
  ADREPORT(relbiom);
  ADREPORT(Expl);
  ADREPORT(Fmort);
  ADREPORT(prodn);
  if (cpue1_nyrs>0) ADREPORT(est_ln_cpue1);
  if (cpue2_nyrs>0) ADREPORT(est_ln_cpue2);
  if (cpue3_nyrs>0) ADREPORT(est_ln_cpue3);
  if (cpue4_nyrs>0) ADREPORT(est_ln_cpue4);
  ADREPORT(est_env);
  if (mod_option==3 || mod_option==4) ADREPORT(prod_without_dep);
  if (mod_option==3 || mod_option==4) ADREPORT(depensation_factor);
  if (mod_type==3) ADREPORT(pt_parm);
  
  return NLL;
  
}
