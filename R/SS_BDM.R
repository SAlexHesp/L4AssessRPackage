#' @useDynLib L4Assess, .registration = TRUE
NULL
#
# L4Assess - State space biomass dynamics model
# Implementation of model described in Marks et al. (2021)
# DPIRD 2023
#
#' Get maximum and minimum x values and x interval
#'
#' Used for setting plot defaults
#'
#' @keywords internal
#'
#' @param x_data x axis data for plot
#'
#' @return list object with minimum and maximum x axis values and x axis interval
Get_xaxis_scale <- function(x_data) {

  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

  xmax_data = max(x_data)
  xmin_data = min(x_data)

  xpow = log10(xmax_data-xmin_data)
  xint = 10 ^ (xpow-round(xpow,0))

  if (xint>=0 & xint<2.5) {
    xint = 0.2
  }
  if (xint>=2.5 & xint<5) {
    xint = 0.5
  }
  if (xint>=5 & xint<7.5) {
    xint = 1
  }
  if (xint>=7.5) {
    xint = 2
  }

  xint = xint * 10^round(xpow,0) # major ticks
  xmin = xint * round(xmin_data/xint,0)
  xmax = xint * (round(xmax_data / xint,0) + 1)

  results = list(xmin = xmin,
                 xmax = xmax,
                 xint = xint)

}

#' Get maximum and minimum y values and y interval
#'
#' Used for setting plot defaults
#'
#' @keywords internal
#'
#' @param y_data y axis data for plot
#'
#' @return list object with minimum and maximum y axis values and y axis interval
Get_yaxis_scale <- function(y_data) {

  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

  ymax_data = max(1.1 * y_data)
  ymin_data = min(y_data)

  ypow = log10(ymax_data-ymin_data)
  yint = 10 ^ (ypow-round(ypow,0))

  if (yint>=0 & yint<2.5) {
    yint = 0.2
  }
  if (yint>=2.5 & yint<5) {
    yint = 0.5
  }
  if (yint>=5 & yint<7.5) {
    yint = 1
  }
  if (yint>=7.5) {
    yint = 2
  }

  yint = yint * 10^round(ypow,0) # major ticks
  ymin = yint * round(ymin_data/yint,0)
  ymax = yint * (round(ymax_data / yint,0) + 1)

  results = list(ymin = ymin,
                 ymax = ymax,
                 yint = yint)

}

#
#' Set initial values for estimated model parameters
#'
#' This function produces a list of initial values for estimated model parameters.
#' The parameters listed depends on number of CPUE series. Note, when fitting models
#' without any link to an environmental index, the value should be set to zero. Set
#' initial values for all parameters, whether estimated or not.
#'
#' @param ln_K natural log of initial value for population carrying capacity.
#' @param ln_r natural log of initial value for population intrinsic increase.
#' @param ln_q1 natural log of initial value for catchability param for 1st CPUE series.
#' @param ln_q2 natural log of initial value for catchability param for 2nd CPUE series.
#' @param ln_q3 natural log of initial value for catchability param for 3rd CPUE series.
#' @param ln_q4 natural log of initial value for catchability param for 4th CPUE series.
#' @param ln_sdmod1 natural log of initial value for sd for 1st CPUE series.
#' @param ln_sdmod2 natural log of initial value for sd for 2nd CPUE series.
#' @param ln_sdmod3 natural log of initial value for sd for 3rd CPUE series.
#' @param ln_sdmod4 initial value for sd for 4th CPUE series. Set to zero for no direct effect of environmental parameter.
#' @param env_param natural log of initial value for environmental param linked to biomass.
#' @param lndep depensation parameter. In normal space, 0=no depensation, 1=full depensation.
#' @param FF exploitation parameters.
#' @param EpsR random effects parameters.
#' @param ln_SigmaR level of process error.
#'
#' @return object containing initial values of model parameters for input into the model
#'
#' @examples
#' nyrs = 44
#' ln_K = log(8000)
#' ln_r = log(0.5)
#' ln_q1 = -10
#' ln_q2 = -10
#' ln_q3 = -10
#' ln_q4 = -10
#' ln_sdmod1 = NA
#' ln_sdmod2 = NA
#' ln_sdmod3 = NA
#' ln_sdmod4 = NA
#' env_param = 0
#' lndep = -1
#' FF = rep(-2, nyrs)
#' EpsR = rep(0, times=nyrs)
#' ln_SigmaR = log(0.3)
#' bdm_params = Set_InitValsForBDM_Params(ln_K, ln_r, ln_q1, ln_q2, ln_q3, ln_q4, ln_sd1, ln_sd2,
#'                                        ln_sd3,ln_sd4, env_p, lndep, ln_pt_parm, FF, EpsR)
#' @export
Set_InitValsForBDM_Params <- function(ln_K, ln_r, ln_q1, ln_q2, ln_q3, ln_q4, ln_sd1, ln_sd2,
                                      ln_sd3,ln_sd4, env_p, lndep, ln_pt_parm, FF, EpsR) {


  bdm_parameters <- list(ln_K=ln_K, ln_r=ln_r, ln_q1=ln_q1, ln_q2=ln_q2, ln_q3=ln_q3,ln_q4=ln_q4,
                         ln_sd1=ln_sd1, ln_sd2=ln_sd2, ln_sd3=ln_sd3,ln_sd4=ln_sd4,
                         env_p=env_p, lndep=lndep, ln_pt_parm=ln_pt_parm,
                         FF=FF, EpsR=EpsR)


  result = bdm_parameters

  return(result)

}

#' Set lower and upper bounds for estimated model parameters
#'
#' This function produces a lists for lower and upper bounds estimated model parameters.
#' Values are specified for all parameters, whether estimated or not.
#'
#' @param ln_K_bnds natural log of initial value for population carrying capacity.
#' @param ln_r_bnds natural log of initial value for population intrinsic increase.
#' @param ln_q_params_bnds natural log of initial value for catchability param for 1st CPUE series.
#' @param ln_sdmod_params_bnds natural log of initial value for catchability param for 2nd CPUE series. Set to NA if only 1 CPUE series.
#' @param env_param_bnds natural log of initial value for catchability param for 3rd CPUE series. Set to NA if 2 CPUE series.
#' @param ln_dep_param_bnds natural log of initial value for catchability param for 4th CPUE series. Set to NA if 3 CPUE series.
#'
#' @return object containing lower and upper bounds of model parameters for input into the model
#'
#' @examples
#' ln_K_bnds <- c(6,9)
#' ln_r_bnds <- log(c(0.01,1.5))
#' ln_q_bnds <- c(-20,0)
#' ln_sd_bnds <- c(-20,1)
#' env_param_bnds <- c(-20,20)
#' ln_dep_bnds <- c(-20,0)
#' ln_pt_bnds <- c(-20,0)
#' bdm_param_bounds = Set_BoundsForBDM_Params(ln_K_bnds, ln_r_bnds, ln_q_bnds, ln_sd_bnds,
#'                                          env_param_bnds, ln_dep_bnds, ln_pt_bnds)
#' @export
Set_BoundsForBDM_Params <- function(ln_K_bnds, ln_r_bnds, ln_q_bnds, ln_sd_bnds,
                                    env_param_bnds, ln_dep_bnds, ln_pt_bnds) {

  # create lists of upper and lower bounds
  low_bound_list <- list(ln_K_low=ln_K_bnds[1],
                         ln_r_low=ln_r_bnds[1],
                         ln_q1_low=ln_q_bnds[1],
                         ln_q2_low=ln_q_bnds[1],
                         ln_q3_low=ln_q_bnds[1],
                         ln_q4_low=ln_q_bnds[1],
                         ln_sd1_low=ln_sd_bnds[1],
                         ln_sd2_low=ln_sd_bnds[1],
                         ln_sd3_low=ln_sd_bnds[1],
                         ln_sd4_low=ln_sd_bnds[1],
                         env_param_low=env_param_bnds[1],
                         lndep_low=ln_dep_bnds[1],
                         ln_pt_parm_low=ln_pt_bnds[1])

  upp_bound_list <- list(ln_K_high=ln_K_bnds[2],
                         ln_r_high=ln_r_bnds[2],
                         ln_q1_high=ln_q_bnds[2],
                         ln_q2_high=ln_q_bnds[2],
                         ln_q3_high=ln_q_bnds[2],
                         ln_q4_high=ln_q_bnds[2],
                         ln_sd1_high=ln_sd_bnds[2],
                         ln_sd2_high=ln_sd_bnds[2],
                         ln_sd3_high=ln_sd_bnds[2],
                         ln_sd4_high=ln_sd_bnds[2],
                         env_param_high=env_param_bnds[2],
                         lndep_high=ln_dep_bnds[2],
                         ln_pt_parm_high=ln_pt_bnds[2])

  result = list(low_bound_list = low_bound_list,
                upp_bound_list = upp_bound_list)

  return(result)

}

#' Get input data for the biomass dynamics model
#'
#' This function produces a list containing all of the data objects required as input for the model
#'
#' @param DatFromCSVFile Data object read in from a csv file, containing all model time series data
#' @param wt_param_pen lambda for parameter penalty
#' @param wt_depl_pen lambda for final depletion penalty
#' @param wt_biom_pen lambda for biomass penalty
#' @param wt_harv_pen lambda for harvest rate penalty
#' @param wt_cpue1 lambda for obs_ln_cpue1 objective function component
#' @param wt_cpue2 lambda for obs_ln_cpue2 objective function component
#' @param wt_cpue3 lambda for obs_ln_cpue3 objective function component
#' @param wt_cpue4 lambda for obs_ln_cpue4 objective function component
#' @param SigmaR specified level of process error associated with biomass production
#' @param Sigma_env specified level of process error associated with environment
#' @param max_currBrel maximum level of relative biomass for current year, used in penalty for final depletion level
#'
#' @return object containing all data required as model input
#'
#' @export
Get_BDM_Data <- function(DatFromCSVFile, wt_param_pen, wt_depl_pen, wt_biom_pen, wt_harv_pen,
                         wt_cpue1, wt_cpue2, wt_cpue3, wt_cpue4, SigmaR, Sigma_env, max_currBrel) {

  dat <- DatFromCSVFile

  # Get the number of years of data
  nyrs <- nrow(dat)

  # Mean of observed environmental index
  mean_obs_env <- mean(DatFromCSVFile$obs_env[which(!is.na(DatFromCSVFile$obs_env))])

  ## Get number of years for each cpue time series
  cpue1_nyrs = length(which(!is.na(dat$cpue1)))
  cpue2_nyrs = length(which(!is.na(dat$cpue2)))
  cpue3_nyrs = length(which(!is.na(dat$cpue3)))
  cpue4_nyrs = length(which(!is.na(dat$cpue4)))

  ## Get first and last year of each cpue time series
  first_cpue1_yr = ifelse(cpue1_nyrs>0, min(dat$year[!is.na(dat$cpue1)]), 0)
  last_cpue1_yr = ifelse(cpue1_nyrs>0, max(dat$year[!is.na(dat$cpue1)]), 0)
  first_cpue2_yr = ifelse(cpue2_nyrs>0, min(dat$year[!is.na(dat$cpue2)]), 0)
  last_cpue2_yr = ifelse(cpue2_nyrs>0, max(dat$year[!is.na(dat$cpue2)]), 0)
  first_cpue3_yr = ifelse(cpue3_nyrs>0, min(dat$year[!is.na(dat$cpue3)]), 0)
  last_cpue3_yr = ifelse(cpue3_nyrs>0, max(dat$year[!is.na(dat$cpue3)]), 0)
  first_cpue4_yr = ifelse(cpue4_nyrs>0, min(dat$year[!is.na(dat$cpue4)]), 0)
  last_cpue4_yr = ifelse(cpue4_nyrs>0, max(dat$year[!is.na(dat$cpue4)]), 0)

  bdm_param_bounds = Set_BoundsForBDM_Params(ln_K_bnds, ln_r_bnds, ln_q_bnds, ln_sd_bnds,
                                             env_param_bnds, ln_dep_bnds, ln_pt_bnds)
  temp_uppbound=bdm_param_bounds$upp_bound_list
  temp_lowbound=bdm_param_bounds$low_bound_list
  uppbound = unlist(temp_uppbound, use.names = F)
  lowbound = unlist(temp_lowbound, use.names = F)

  bdm_data <- list(mod_type=NA,
                   mod_option=NA,
                   mod_scenario=NA,
                   nyrs=nyrs,
                   cpue1_nyrs=cpue1_nyrs,
                   cpue2_nyrs=cpue2_nyrs,
                   cpue3_nyrs=cpue3_nyrs,
                   cpue4_nyrs=cpue4_nyrs,
                   first_cpue1_yr=first_cpue1_yr,
                   last_cpue1_yr=last_cpue1_yr,
                   first_cpue2_yr=first_cpue2_yr,
                   last_cpue2_yr=last_cpue2_yr,
                   first_cpue3_yr=first_cpue3_yr,
                   last_cpue3_yr=last_cpue3_yr,
                   first_cpue4_yr=first_cpue4_yr,
                   last_cpue4_yr=last_cpue4_yr,
                   season=dat$year,
                   tot_catch=dat$tot_catch,
                   obs_ln_cpue1=dat$obs_ln_cpue1,
                   obs_ln_cpue2=dat$obs_ln_cpue2,
                   obs_ln_cpue3=dat$obs_ln_cpue3,
                   obs_ln_cpue4=dat$obs_ln_cpue4,
                   cpue1_se=dat$cpue1_se,
                   cpue2_se=dat$cpue2_se,
                   cpue3_se=dat$cpue3_se,
                   cpue4_se=dat$cpue4_se,
                   obs_env=dat$obs_env,
                   env_se=dat$env_se,
                   mean_obs_env=mean_obs_env,
                   SigmaR=SigmaR,
                   Sigma_env=Sigma_env,
                   max_currBrel=max_currBrel,
                   wt_param_pen=wt_param_pen,
                   wt_depl_pen=wt_depl_pen,
                   wt_biom_pen=wt_biom_pen,
                   wt_harv_pen=wt_harv_pen,
                   wt_cpue1=wt_cpue1,
                   wt_cpue2=wt_cpue2,
                   wt_cpue3=wt_cpue3,
                   wt_cpue4=wt_cpue4,
                   uppbound=uppbound,
                   lowbound=lowbound)

  return(bdm_data)

}

#' Get model template model builder (TMB) parameter map
#'
#' This function produces an list, listing model parameters that should not be estimated using TMB. Note
#' that for the second scenario, EpsR are 'estimated' as random effects.
#'
#' @param DatFromCSVFile Data object read in from a csv file, containing all model time series data
#' @param mod_scenario 1=non-state space, 2=state space
#' @param mod_option 1=traditional, 2= +chl, 3= +dep, 3 = +chl and dep
#'
#' @return List of model parameters that should not be estimated using TMB
#'
#' @examples
#' mod_scenario = 2 # 1=non-state space, 2=state space
#' model_option = 1 # 1=traditional, 2= +chl, 3= +dep, 3 = +chl and dep
#' nyrs = 44
#' Get_BDM_Map(DatFromCSVFile, mod_option, mod_scenario)
#' @export
Get_BDM_Map <- function(DatFromCSVFile, mod_scenario, mod_option) {


  bdm_map <- NULL
  # if model type is schaefer or Fox
  if (mod_scenario == 1 & mod_type %in% c(1,2)) {
    bdm_map <- list(ln_pt_parm=factor(NA),env_p=factor(NA), lndep=factor(NA), ln_sd1=factor(NA),
                    ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                    EpsR=rep(factor(NA), times=nyrs))
    if(mod_option==2 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA),lndep=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
    if(mod_option==3 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA),env_p=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
    if(mod_option==4 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
  }
  # if model type is PT
  if(mod_scenario == 1 & mod_type==3){
    bdm_map <- list(env_p=factor(NA), lndep=factor(NA),ln_sd1=factor(NA), ln_pt_parm=factor(NA),
                    ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                    EpsR=rep(factor(NA), times=nyrs))
    if(mod_option==2  & mod_type==3) {
      bdm_map <- list(lndep=factor(NA), ln_sd1=factor(NA), ln_pt_parm=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
    if(mod_option==3  & mod_type==3) {
      bdm_map <- list(env_p=factor(NA), ln_sd1=factor(NA), ln_pt_parm=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
    if(mod_option==4  & mod_type==3) {
      bdm_map <- list(ln_sd1=factor(NA),
                      ln_pt_parm=factor(NA),
                      ln_q2=factor(NA),ln_q3=factor(NA),ln_q4=factor(NA),
                      EpsR=rep(factor(NA), times=nyrs))
    }
  }
  # if model type is schaefer or Fox
  if(mod_scenario==2 & mod_type %in% c(1,2)) {

    bdm_map <- list(ln_pt_parm=factor(NA), ln_sd1=factor(NA),
                    ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                    env_p=factor(NA), lndep=factor(NA))

    if(mod_option==2 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      lndep=factor(NA))
    }
    if(mod_option==3 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
                      env_p=factor(NA))
    }
    if(mod_option==4 & mod_type %in% c(1,2)) {
      bdm_map <- list(ln_pt_parm=factor(NA), ln_sd1=factor(NA),
                      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA)) # all params included
    }
  }
  # if model type is PT
  if(mod_scenario==2 & mod_type==3) {

    # For mod_option 1 and 3, set env_param to zero
    bdm_map <- list(
      ln_sd1=factor(NA),ln_pt_parm=factor(NA),
      ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
      env_p=factor(NA), lndep=factor(NA))
    if(mod_option==2 & mod_type==3) {
      bdm_map <- list(
        ln_sd1=factor(NA), ln_pt_parm=factor(NA),
        ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
        lndep=factor(NA))
    }
    if(mod_option==3 & mod_type==3) {
      bdm_map <- list(
        ln_sd1=factor(NA), ln_pt_parm=factor(NA),
        ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA),
        env_p=factor(NA))
    }
    if(mod_option==4 & mod_type==3) {
      bdm_map <- list(
        ln_sd1=factor(NA), ln_pt_parm=factor(NA),
        ln_sd2=factor(NA),ln_sd3=factor(NA), ln_sd4=factor(NA)) # all params included
    }
  }

  dat = DatFromCSVFile
  nCPUE1Obs = length(which(!is.na(dat$obs_ln_cpue1)))
  if (nCPUE1Obs==0) bdm_map$ln_q1=factor(NA)
  nCPUE2Obs = length(which(!is.na(dat$obs_ln_cpue2)))
  if (nCPUE2Obs==0) bdm_map$ln_q2=factor(NA)
  nCPUE3Obs = length(which(!is.na(dat$obs_ln_cpue3)))
  if (nCPUE3Obs==0) bdm_map$ln_q3=factor(NA)
  nCPUE4Obs = length(which(!is.na(dat$obs_ln_cpue4)))
  if (nCPUE4Obs==0) bdm_map$ln_q4=factor(NA)

  return(bdm_map)

}

#' Get values calculated by ADREPORT statements in the TMB C++ code
#'
#' This function sets up the TMB function to avoid NLL from fn() differing from
#' report and failing to match calculated value with penalties. The values calculated by ADREPORT
#' statements in the TMB C++ code are returned. These may be erroneous if the model has not converged.
#' Note, obj = obj_TMB.
#'
#' @param obj object containing TMB output
#' @return object containing TMB output
#' @export
GetSDrep <- function(obj) {

  # Set up TMB function to avoid NLL from fn() differing from report and failing to
  # match calculated value with penalties. Collect the values calculated by ADREPORT
  # statements in the TMB C++ code. These may be erroneous if the model has not converged
  # obj = obj_TMB

  result <- tryCatch(
    {
      # 'tryCatch()' will return the last evaluated expression in case the "try" part was completed successfully
      sdreport(obj)

      # The return value of `readLines()` is the actual value that will be returned in case there is no condition
      # (e.g. warning or error). You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped inside a function (unlike that for the condition handlers
      # for warnings and error below)
    },
    error=function(cond)
    {
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond)
    {
      # Choose a return value in case of warning
      return(NULL)
    },
    finally=
      {
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>'
        # message(paste("Processed URL:", url))
        # message("Some other message at the end")
      })
  return(result)
}

#' Update the full set of parameters for the model with the values estimated in TMB
#'
#' This function updates the full set of parameters for the model with the values estimated in TMB,
#' from the TMB object obj$env$last.par.best
#'
#' @keywords internal
#' @param obj object containing TMB output
#' @return object containing TMB output
update_full_parameter_set <- function(full_parameter_set, best_param_estimates) {

  new_full_parameter_set <- full_parameter_set
  names_of_new_full_parameter_set <- names(new_full_parameter_set)
  names_of_best_param_estimates <- names(best_param_estimates)
  for (j in 1:length(names_of_new_full_parameter_set)) {
    ptr <- which(names_of_best_param_estimates == names_of_new_full_parameter_set[j])
    if (length(ptr) == 0) {
      # This parameter from the full set of parameters is not being estimated.
    } else {
      if (length(ptr) == 1) {
        # The current estimate of parameter j in the full parameter set is a scalar and is entry ptr
        # in the vector, best_param_estimates
        new_full_parameter_set[[j]] <- best_param_estimates[ptr]
      } else {
        new_full_parameter_set[[j]] <- as.numeric(best_param_estimates[ptr[1]:ptr[length(ptr)]])
      }
    }
  }
  return(new_full_parameter_set)
}


#' This function returns a range of outputs for the fitted model
#'
#' This function returns a range of outputs for the fitted model (required for plotting)
#'
#' @param result results of fitted TMB model
#' @param nyrs number of years of data
#'
#' @return model_outputs, object containing TMB outputs for plotting
#'
#' @export
Get_model_outputs <- function(result, nyrs) {

  result_TMB <- result$rep_TMB
  obj_TMB <- result$obj_TMB
  sdrep <- GetSDrep(obj_TMB)
  sdrep_outputs <- as.data.frame(summary(sdrep)) # shows all values specified in cpp file
  sdrep_outputs$name <- row.names(sdrep_outputs)
  sdrep_outputs$upp <- sdrep_outputs$Estimate + (1.96*sdrep_outputs$`Std. Error`)
  sdrep_outputs$low <- sdrep_outputs$Estimate - (1.96*sdrep_outputs$`Std. Error`)

  model_outputs <- data.frame(season=result_TMB$season,
                        Chat=result_TMB$Chat[1:nyrs],
                        Chat_upp=sdrep_outputs$upp[grep("Chat",sdrep_outputs$name)][1:nyrs],
                        Chat_low=sdrep_outputs$low[grep("Chat",sdrep_outputs$name)][1:nyrs],
                        biom=result_TMB$biom[1:nyrs],
                        biom_upp=sdrep_outputs$upp[grep("biom",sdrep_outputs$name)][1:nyrs],
                        biom_low=sdrep_outputs$low[grep("biom",sdrep_outputs$name)][1:nyrs],
                        relbiom=result_TMB$relbiom[1:nyrs],
                        relbiom_upp=sdrep_outputs$upp[grep("relbiom",sdrep_outputs$name)][1:nyrs],
                        relbiom_low=sdrep_outputs$low[grep("relbiom",sdrep_outputs$name)][1:nyrs],
                        prodn=result_TMB$prodn[1:nyrs],
                        prodn_upp=sdrep_outputs$upp[grep("prodn",sdrep_outputs$name)][1:nyrs],
                        prodn_low=sdrep_outputs$low[grep("prodn",sdrep_outputs$name)][1:nyrs],
                        obs_ln_cpue1=result_TMB$obs_ln_cpue1[1:nyrs],
                        obs_ln_cpue2=result_TMB$obs_ln_cpue2[1:nyrs],
                        obs_ln_cpue3=result_TMB$obs_ln_cpue3[1:nyrs],
                        obs_ln_cpue4=result_TMB$obs_ln_cpue3[1:nyrs],
                        est_cpue1=result_TMB$est_ln_cpue1[1:nyrs],
                        est_cpue1_upp=sdrep_outputs$upp[grep("est_ln_cpue1",sdrep_outputs$name)][1:nyrs],
                        est_cpue1_low=sdrep_outputs$low[grep("est_ln_cpue1",sdrep_outputs$name)][1:nyrs],
                        est_cpue2=result_TMB$est_ln_cpue2[1:nyrs],
                        est_cpue2_upp=sdrep_outputs$upp[grep("est_ln_cpue2",sdrep_outputs$name)][1:nyrs],
                        est_cpue2_low=sdrep_outputs$low[grep("est_ln_cpue2",sdrep_outputs$name)][1:nyrs],
                        est_cpue3=result_TMB$est_ln_cpue3[1:nyrs],
                        est_cpue3_upp=sdrep_outputs$upp[grep("est_ln_cpue3",sdrep_outputs$name)][1:nyrs],
                        est_cpue3_low=sdrep_outputs$low[grep("est_ln_cpue3",sdrep_outputs$name)][1:nyrs],
                        est_cpue4=result_TMB$est_ln_cpue4[1:nyrs],
                        est_cpue4_upp=sdrep_outputs$upp[grep("est_ln_cpue4",sdrep_outputs$name)][1:nyrs],
                        est_cpue4_low=sdrep_outputs$low[grep("est_ln_cpue4",sdrep_outputs$name)][1:nyrs],
                        Fmort=sdrep_outputs$Estimate[grep("Fmort",sdrep_outputs$name)][1:nyrs],
                        Fmort_upp=sdrep_outputs$upp[grep("Fmort",sdrep_outputs$name)][1:nyrs],
                        Fmort_low=sdrep_outputs$low[grep("Fmort",sdrep_outputs$name)][1:nyrs],
                        Expl=sdrep_outputs$Estimate[grep("Expl",sdrep_outputs$name)][1:nyrs],
                        Expl_upp=sdrep_outputs$upp[grep("Expl",sdrep_outputs$name)][1:nyrs],
                        Expl_low=sdrep_outputs$low[grep("Expl",sdrep_outputs$name)][1:nyrs],
                        obs_chl=result_TMB$obs_env[1:nyrs],
                        est_env=sdrep_outputs$Estimate[grep("est_env",sdrep_outputs$name)][1:nyrs],
                        est_env_upp=sdrep_outputs$upp[grep("est_env",sdrep_outputs$name)][1:nyrs],
                        est_env_low=sdrep_outputs$low[grep("est_env",sdrep_outputs$name)][1:nyrs],
                        EpsR=sdrep_outputs$Estimate[grep("EpsR",sdrep_outputs$name)][1:nyrs],
                        EpsR_upp=sdrep_outputs$upp[grep("EpsR",sdrep_outputs$name)][1:nyrs],
                        EpsR_low=sdrep_outputs$low[grep("EpsR",sdrep_outputs$name)][1:nyrs],
                        FF=sdrep_outputs$Estimate[grep("FF",sdrep_outputs$name)][1:nyrs],
                        FF_upp=sdrep_outputs$upp[grep("FF",sdrep_outputs$name)][1:nyrs],
                        FF_low=sdrep_outputs$low[grep("FF",sdrep_outputs$name)][1:nyrs])

  return(model_outputs)

}


#' This function returns a range of model parameter outputs for the fitted model
#'
#' This function returns a range of model parameter outputs for the fitted model
#'
#' @param result results of fitted TMB model
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#'
#' @return outputs, object containing TMB outputs for plotting
#'
#' @export
#'
Get_model_parameter_outputs <- function(result, mod_type, mod_option) {

  fit_TMB <- result$fit_TMB
  obj_TMB <- result$obj_TMB
  sdrep <- GetSDrep(obj_TMB)
  sdrep_outputs <- as.data.frame(summary(sdrep)) # shows all values specified in cpp file
  sdrep_outputs$name <- row.names(sdrep_outputs)
  sdrep_outputs$upp <- sdrep_outputs$Estimate + (1.96*sdrep_outputs$`Std. Error`)
  sdrep_outputs$low <- sdrep_outputs$Estimate - (1.96*sdrep_outputs$`Std. Error`)

  param_results <- data.frame(NLL=fit_TMB$objective,
                              NLL_CPUE1=obj_TMB$report()$NLL_CPUE1,
                              NLL_CPUE2=obj_TMB$report()$NLL_CPUE2,
                              NLL_CPUE3=obj_TMB$report()$NLL_CPUE3,
                              NLL_CPUE4=obj_TMB$report()$NLL_CPUE4,
                              NLL_Chat=obj_TMB$report()$NLL_Chat,
                              NLL_EpsR=obj_TMB$report()$EpsR_pen,
                              ln_K=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_K"],
                              K=sdrep_outputs$Estimate[sdrep_outputs$name=="K"],
                              K_upp=sdrep_outputs$upp[sdrep_outputs$name=="K"],
                              K_low=sdrep_outputs$low[sdrep_outputs$name=="K"],
                              ln_r=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_r"],
                              r=sdrep_outputs$Estimate[sdrep_outputs$name=="r"],
                              r_upp=sdrep_outputs$upp[sdrep_outputs$name=="r"],
                              r_low=sdrep_outputs$low[sdrep_outputs$name=="r"],
                              ln_q1=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q1"],
                              ln_q1_upp=sdrep_outputs$upp[sdrep_outputs$name=="ln_q1"],
                              ln_q1_low=sdrep_outputs$low[sdrep_outputs$name=="ln_q1"],
                              q1=sdrep_outputs$Estimate[sdrep_outputs$name=="q1"],
                              q1_upp=sdrep_outputs$upp[sdrep_outputs$name=="q1"],
                              q1_low=sdrep_outputs$low[sdrep_outputs$name=="q1"])

  if (length(sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q2"])==1) {
    param_results$ln_q2=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q2"]
    param_results$ln_q2_upp=sdrep_outputs$upp[sdrep_outputs$name=="ln_q2"]
    param_results$ln_q2_low=sdrep_outputs$low[sdrep_outputs$name=="ln_q2"]
    param_results$q2=sdrep_outputs$Estimate[sdrep_outputs$name=="q2"]
    param_results$q2_upp=sdrep_outputs$upp[sdrep_outputs$name=="q2"]
    param_results$q2_low=sdrep_outputs$low[sdrep_outputs$name=="q2"]
  }
  if (length(sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q3"])==1) {
    param_results$ln_q3=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q3"]
    param_results$ln_q3_upp=sdrep_outputs$upp[sdrep_outputs$name=="ln_q3"]
    param_results$ln_q3_low=sdrep_outputs$low[sdrep_outputs$name=="ln_q3"]
    param_results$q3=sdrep_outputs$Estimate[sdrep_outputs$name=="q3"]
    param_results$q3_upp=sdrep_outputs$upp[sdrep_outputs$name=="q3"]
    param_results$q3_low=sdrep_outputs$low[sdrep_outputs$name=="q3"]
  }
  if (length(sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q4"])==1) {
    param_results$ln_q4=sdrep_outputs$Estimate[sdrep_outputs$name=="ln_q4"]
    param_results$ln_q4_upp=sdrep_outputs$upp[sdrep_outputs$name=="ln_q4"]
    param_results$ln_q4_low=sdrep_outputs$low[sdrep_outputs$name=="ln_q4"]
    param_results$q4=sdrep_outputs$Estimate[sdrep_outputs$name=="q4"]
    param_results$q4_upp=sdrep_outputs$upp[sdrep_outputs$name=="q4"]
    param_results$q4_low=sdrep_outputs$low[sdrep_outputs$name=="q4"]
  }


  if (mod_option == 2 | mod_option == 4) {
    param_results$env_p=sdrep_outputs$Estimate[sdrep_outputs$name=="env_param"]
    param_results$env_p_upp=sdrep_outputs$upp[sdrep_outputs$name=="env_param"]
    param_results$env_p_low=sdrep_outputs$low[sdrep_outputs$name=="env_param"]
  }

  if (mod_option == 3 | mod_option == 4) {
    param_results$dep=sdrep_outputs$Estimate[sdrep_outputs$name=="dep"]
    param_results$dep_upp=sdrep_outputs$upp[sdrep_outputs$name=="dep"]
    param_results$dep_low=sdrep_outputs$low[sdrep_outputs$name=="dep"]
  }

  return(param_results)

}

#' Fit the biomass dynamics model using TMB
#'
#' The is a wrapper function used to fit the biomass dynamics model using TMB.
#' The model is always first fitted as a non-state space model, to get better
#' starting values of key parameters before fitting the state space model.
#'
#' @param DatFromCSVFile Data object read in from a csv file, containing all model time series data
#' @param mod_scenario 1=non-state space, 2=state space
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#' @param bdm_data object containing input data for TMB model
#' @param bdm_params object containing parameters for TMB model
#' @param bdm_map object containing parameter map for TMB model
#' @param initial_params optional input to input starting values when fitting state space model
#'
#' @return multiple objects with output from TMB  (rep_TMB, fit_TMB, obj_TMB, best_par_TMB, best_fitted_par_TMB, best_val_TMB, new_bdm_map)
#' @export
fit_the_model <- function(DatFromCSVFile, mod_scenario, mod_type, mod_option,
                          bdm_data, bdm_params, bdm_map=NULL,
                          initial_params=NULL, nyrs) {

  new_bdm_data <- bdm_data

  # Set up the model scenario, type and option
  new_bdm_data$mod_type <- mod_type
  new_bdm_data$mod_option <- mod_option
  new_bdm_data$mod_scenario <- mod_scenario

  # Get the map if it has not been provided when the function was called. Note that,
  # if provided, it is possible that the list is empty, i.e. bdm_map==list()
  new_bdm_map <- bdm_map
  if (is.null(bdm_map)) {
    # A map list was not provided. Determine the list from the model option and scenario
    new_bdm_map <- Get_BDM_Map(DatFromCSVFile, mod_scenario, mod_option)
  }

  cat("Map used in fitting the model: \n")
  print(new_bdm_map)
  cat("\n")

  # Set up the TMB model, assuming that all parameters are fixed effects
  # add new bdm data and new bdm map here
  BDM_model <- MakeADFun(data=new_bdm_data,
                         parameters=bdm_params,
                         DLL="SSBDM", silent=T, map=new_bdm_map,
                         control=list(eval.max=10000,iter.max=2000, rel.tol=1e-15))

  # Set up the initial values of the parameters, noting that, if the actual model to be fitted
  # includes the random effects, EpsR, then it is assumed that these will have been
  # provided in the argument, initial_params, if this is supplied when the function,
  # fit_the_model, is called.
  newpar <- BDM_model$par
  if (!is.null(initial_params)) {
    newpar <- initial_params
  }

  cat("Initial parameter estimates: \n")
  print(newpar)
  cat("\n")

  # Get improved estimates for the parameters. This will assist when
  # fitting a random effects model
  fit_TMB <- suppressWarnings(nlminb(newpar, BDM_model$fn, BDM_model$gr,
                                     control=list(trace=100,eval.max=10000,iter.max=2000)))

  cat("Results of first nlminb fit using TMB functions: \n")
  print(fit_TMB)
  cat("\n")

  # Set newpar to the improved parameter estimates
  newpar <- fit_TMB$par

  # Set the indices that point to the parameters for the current fixed effects, i.e.
  # all current parameters
  Index <- 1:length(newpar)

  # Update the full set of parameters for the model with the values from obj$env$last.par.best
  best_param_estimates <- BDM_model$env$last.par.best
  new_bdm_parameters <- update_full_parameter_set(bdm_params,best_param_estimates)

  # If the model to be fitted is actually a random effects model, recreate the TMB
  # model using the improved parameter estimates

  # new map and new data here
  if (mod_scenario == 1) {
    BDM_model <- MakeADFun(data=new_bdm_data,
                           parameters=new_bdm_parameters,
                           DLL="SSBDM",
                           silent=T, map=new_bdm_map,
                           control=list(eval.max=10000,iter.max=2000, rel.tol=1e-15))
  } else {
    # Re-create the TMB model
    BDM_model <- MakeADFun(data=new_bdm_data,
                           parameters=new_bdm_parameters,
                           DLL="SSBDM",
                           silent=T, map=new_bdm_map, random="EpsR",
                           control=list(eval.max=10000,iter.max=2000, rel.tol=1e-15))
  }

  # Re-create newpar
  newpar <- BDM_model$par

  # Evaluate the function to set up a new value of obj$env$last.par.best
  BDM_model$fn(newpar)

  if (mod_scenario == 2) {
    # Reset the indices that point to the parameters for the fixed effects in
    # the random effects model
    Index <- which(names(BDM_model$env$last.par.best)!="EpsR")
  }

  # Fit the model using the TMB functions -- excl EpsR
  for (i in 1:1) {
    fit_TMB <- suppressWarnings(nlminb(BDM_model$env$last.par.best[Index],
                                       BDM_model$fn, BDM_model$gr,
                                       control=list(trace=100,eval.max=10000, iter.max=2000)))
  }

  cat("\n")
  cat("Result of fitting the model a further 10 times using the TMB functions: \n")
  print(fit_TMB)
  cat("\n")

  # Get the best parameters
  best_par_TMB <- BDM_model$env$last.par.best
  best_fitted_par_TMB <- BDM_model$env$last.par.best[Index]
  best_val_TMB <- BDM_model$env$value.best   # NLL
  obj_TMB <- BDM_model

  best_param_estimates <- BDM_model$env$last.par.best
  new_bdm_parameters <- update_full_parameter_set(new_bdm_parameters,best_param_estimates)

  # report from TMB with biomass estimates etc. NLL NOT CORRECT
  rep_TMB <- BDM_model$report(BDM_model$env$last.par.best)

  cat("The parameters that provided the best fit with TMB functions were : \n")
  print(best_par_TMB)
  cat("Value of objective function: ", best_val_TMB, "\n")
  cat("Value of objective function returned by fn(): ",
      BDM_model$fn(best_fitted_par_TMB), "\n")
  cat("\n")

  # return list of interest to output
  retlist <- list(rep_TMB=rep_TMB,
                  fit_TMB=fit_TMB,
                  obj_TMB=obj_TMB,
                  best_par_TMB=best_par_TMB,
                  best_fitted_par_TMB=best_fitted_par_TMB,
                  best_val_TMB=best_val_TMB,
                  new_bdm_map=new_bdm_map)

  return(retlist)
}


#' Calculate MSY and BMSY
#'
#' This function calculates maximum sustainable yield (MSY) and biomass at MSY (BMSY),
#' for each of the model types and model options
#'
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#' @param param_results object containing results outputted by for model by TMB
#' @param env standardised environmental index value. 0=average conditions.
#'
#' @return Bmsy1, MSY1, F_lim, F_thresh, F_targ, B_lim, B_thresh, B_targ
#'
#' @export
Calc_MSYAndBiolRefPoints <- function(mod_type, mod_option, param_results, env){

  # specify model number
  model_num <- (mod_type-1) * 4 + mod_option

  r <- param_results$r
  K <- param_results$K
  pt_parm <- param_results$pt_parm
  g <- param_results$env_p
  d <- param_results$dep
  PellaMSY_a <- (pt_parm + 1.0)
  PellaMSY_b <- -(pt_parm + 1.0)/pt_parm

  ## Define function for production
  ## Note x here means B i.e. biomass
  prod_func = switch(model_num,
                     function(x) {r*x*(1-x/K)},
                     function(x) {r*x*(1-x/K)*exp(g*env)},
                     function(x) {r*x*(1-x/K)*(1-exp(log(0.5)*x/(d*K)))},
                     function(x) {r*x*(1-x/K)*(1-exp(log(0.5)*x/(d*K)))*exp(g*env)},
                     function(x) {r*log(K)*x*(1-log(x)/log(K))},
                     function(x) {r*log(K)*x*(1-log(x)/log(K))*exp(g*env)},
                     function(x) {r*log(K)*x*(1-log(x)/log(K))*(1-exp(log(0.5)*x/(d*K)))},
                     function(x) {r*log(K)*x*(1-log(x)/log(K))*(1-exp(log(0.5)*x/(d*K)))*exp(g*env)},
                     function(x) {(r/pt_parm)*x*(1-(x/K)^pt_parm)},
                     function(x) {r*x*(1-x/K)*exp(g*env)},
                     function(x) {(r/pt_parm)*x*(1-(x/K)^pt_parm)*(1-exp(log(0.5)*x/(d*K)))},
                     function(x) {(r/pt_parm)*x*(1-(x/K)^pt_parm)*(1-exp(log(0.5)*x/(d*K)))*exp(g*env)})


  ## Define function for derivative of production function (i.e. function that defines slope)
  ## Note x here means B i.e. biomass
  dprod_func = switch(model_num,
                      function(x) {r*(1-2*x/K)},
                      function(x) {r*(1-2*x/K)*exp(g*env)},
                      function(x) {(r*(1-2*x/K) + exp(log(0.5)*x/(d*K)) * (r*(2*x/K - 1) + log(0.5)*r*x*(x-K) / (d*K^2)))},
                      function(x) {(r*(1-2*x/K) + exp(log(0.5)*x/(d*K)) * (r*(2*x/K - 1) + log(0.5)*r*x*(x-K) / (d*K^2))) *exp(g*env)},
                      function(x) {r*(log(K/x)-1)},
                      function(x) {r*(log(K/x)-1)*exp(g*env)},
                      function(x) {(r*(log(K/x)-1)*(1-exp(log(0.5)*x/(d*K))) -
                                      log(K)*log(0.5)*r*x/(d*K)*(1-log(x)/log(K))*exp(log(0.5)*x/(d*K)))},
                      function(x) {(r*(log(K/x)-1)*(1-exp(log(0.5)*x/(d*K))) -
                                      log(K)*log(0.5)*r*x/(d*K)*(1-log(x)/log(K))*exp(log(0.5)*x/(d*K)))*exp(g*env)})

  Bmsy1 = switch(model_num,
                 K/2,
                 K/2,
                 uniroot(dprod_func, lower = K/10, upper = K, tol = 1e-9)$root,
                 uniroot(dprod_func, lower = K/10, upper = K, tol = 1e-9)$root,
                 K/exp(1),
                 K/exp(1),
                 uniroot(dprod_func, lower = K/10, upper = K, tol = 1e-9)$root,
                 uniroot(dprod_func, lower = K/10, upper = K, tol = 1e-9)$root)


  MSY1 = switch(model_num,
                r*K/4,
                r*K/4*exp(g*env),
                prod_func(Bmsy1),
                prod_func(Bmsy1),
                r*K/exp(1),
                r*K/exp(1)*exp(g*env),
                prod_func(Bmsy1),
                prod_func(Bmsy1),
                (r*K)*PellaMSY_a^PellaMSY_b,
                ((r*K)*PellaMSY_a^PellaMSY_b)*exp(g*env))

  if(model_num %in% c(1:4,9:12)){ # check - 9:12 for PT model
    F_thresh=param_results$r*0.5
    F_targ <- param_results$r*0.4
    F_lim <- param_results$r*0.75
    B_thresh <- 0.5
    B_lim <- 0.25
    B_targ <- 1.2*B_thresh
  }

  if(model_num %in% c(5:8)){
    F_thresh=param_results$r
    F_targ <- param_results$r*0.82
    F_lim <- param_results$r*1.69
    B_thresh <- 1/exp(1)
    B_lim <- 0.5*B_thresh
    B_targ <- 1.2*B_thresh
  }

  results = list(Bmsy = Bmsy1,
                 MSY = MSY1,
                 F_lim = F_lim,
                 F_thresh = F_thresh,
                 F_targ = F_targ,
                 B_lim = B_lim,
                 B_thresh = B_thresh,
                 B_targ = B_targ)

  return(results)

}

#' Output  MSY and BMSY estimates from parametric resampling
#'
#' This function outputs MSY and BMSY estimates from parametric resampling, using model parameter values
#' derived by passing  random values of model parameters, from the estimated variance-covariance matrix, through
#' MSY and BMSY calculutions, assuming the estimated parameters conform to a multivariate normal distribution
#'
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#' @param ResampRes object containing results outputted by for model by TMB, after resampling
#' @param env standardised environmental index value. 0=average conditions.
#'
#' @return MSY_sim_res, estimates of MSY and BMSY, including mean, median, lower and upper 60 and 95 percentiles
#'
#' @export
Get_MSY_and_BMSY_with_err <- function(mod_type, mod_option, ResampRes, env) {

  MSY_sim <- rep(NA,nreps)
  Bmsy_sim <- rep(NA,nreps)
  for (i in 1:nreps) {
    r <- ResampRes$params$r[i]
    K <- ResampRes$params$K[i]
    pt_parm <- ResampRes$params$pt_parm[i]
    env_p <- ResampRes$params$env_p[i]
    dep <- ResampRes$params$dep[i]
    param_results = list(r=r,K=K,pt_parm=pt_parm,env_p=env_p,dep=dep)
    MSYRes = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env)
    MSY_sim[i] <- MSYRes$MSY
    Bmsy_sim[i] <- MSYRes$Bmsy
  }
  MSY_resamp_mean = mean(MSY_sim)
  MSY_resamp_med = median(MSY_sim)
  MSY_resamp_lw95 = quantile(MSY_sim,0.025)
  MSY_resamp_up95 = quantile(MSY_sim,0.975)
  MSY_resamp_lw60 = quantile(MSY_sim,0.2)
  MSY_resamp_up60 = quantile(MSY_sim,0.8)
  Bmsy_resamp_mean = mean(Bmsy_sim)
  Bmsy_resamp_med = median(Bmsy_sim)
  Bmsy_resamp_lw95 = quantile(Bmsy_sim,0.025)
  Bmsy_resamp_up95 = quantile(Bmsy_sim,0.975)
  Bmsy_resamp_lw60 = quantile(Bmsy_sim,0.2)
  Bmsy_resamp_up60 = quantile(Bmsy_sim,0.8)

  MSY_sim_res = data.frame(MSY_mean=MSY_resamp_mean,
                           MSY_med=MSY_resamp_med,
                           MSY_lw95=MSY_resamp_lw95,
                           MSY_up95=MSY_resamp_up95,
                           MSY_lw60=MSY_resamp_lw60,
                           MSY_up60=MSY_resamp_up60,
                           Bmsy_mean=Bmsy_resamp_mean,
                           Bmsy_med=Bmsy_resamp_med,
                           Bmsy_lw95=Bmsy_resamp_lw95,
                           Bmsy_up95=Bmsy_resamp_up95,
                           Bmsy_lw60=Bmsy_resamp_lw60,
                           Bmsy_up60=Bmsy_resamp_up60,
                           row.names = '')

  MSY_sim_res <- round(MSY_sim_res,2)

  return(MSY_sim_res)

}

#' Plot biomass, catches and exploitation (asymptotic error)
#'
#' This function plots relative and absolute biomass, catch and exploitation
#' and associated 95 percent confidence limits (derived from asymptotic errors)
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#' @param model_outputs object containing outputs from fitted TMB model
#' @param param_results object containing outputs from fitted TMB model
#' @param xaxis_lab x axis label, to overide default
#'
#' @return plots of biomass, catch and exploitation
#'
#' @export
Plot_Biomass_Catch_And_Exploitation <- function(DatFromCSVFile, mod_type, mod_option, model_outputs, param_results, xaxis_lab) {

  dat=DatFromCSVFile
  year <- model_outputs$season

  # specify model number
  model_num <- (mod_type-1) * 4 + mod_option

  red = as.numeric(col2rgb(2)) / 255
  orange = as.numeric(col2rgb("orange")) / 255
  yellow = as.numeric(col2rgb("lightgoldenrod2")) / 255
  green = as.numeric(col2rgb(3)) / 255

  par(mfrow=c(2,2), mar=c(3,4,1,2), mgp=c(2, 0.5, 0), oma=c(0,0,2,0))

  # (1) Biomass----
  if (is.na(xaxis_lab)) xaxis_lab = "Year"
  ylims = Get_yaxis_scale(model_outputs$biom_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(model_outputs$season, model_outputs$biom, type="l", lwd=1, ylim=c(0,ymax),
       xlab=xaxis_lab, ylab="", las=1)
  x <- c(year,rev(year))
  y <- c(model_outputs$biom_low,rev(model_outputs$biom_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  mtext("Biomass, t",side=2,line=3,cex=0.7)
  legend("topleft",paste("a)"),bty="n")

  # (2) Relative Biomass----
  B_lim = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_lim
  B_thresh = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_thresh
  B_targ = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_targ

  ylims = Get_yaxis_scale(model_outputs$relbiom_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(year,model_outputs$relbiom, type="l", ylim=c(0,ymax), lwd=1,
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Relative Biomass", side=2,line=2.5,cex=0.7)
  x <- c(year,rev(year))
  y <- c(model_outputs$relbiom_low,rev(model_outputs$relbiom_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  abline(h=c(B_lim, B_thresh, B_targ), lty=1, col=c(2, "lightgoldenrod2", 3))
  legend("topleft",paste("b)"),bty="n")

  ## (3) catch and MSY ----
  ylims = Get_yaxis_scale(model_outputs$Chat_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(dat$year,dat$tot_catch, type="p", ylim=c(0,ymax), pch=21,
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Catch, t",side=2,line=2.5,cex=0.7)
  lines(model_outputs$season,model_outputs$Chat, type="l", lwd=1)
  x <- c(year,rev(year))
  y <- c(model_outputs$Chat_low,rev(model_outputs$Chat_upp))
  polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
  abline(h=c(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$MSY), lty=c(1,3,3), col=c(2,2,2))
  legend("topright",paste("MSY ", round(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$MSY),"t"), bty="n")
  legend("topleft",paste("c)"),bty="n")

  ## (4) exploitation ----
  ylims = Get_yaxis_scale(model_outputs$Expl_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(model_outputs$season,model_outputs$Expl, type="l", lwd=1, ylim=c(0,ymax),
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Exploitation",side=2,line=2.5,cex=0.7)
  x <- c(year,rev(year))
  y <- c(model_outputs$Expl_low,rev(model_outputs$Expl_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  # not plotting, as high uncertainty for r estimate may mean these values are misleading
  # abline(h=c(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_lim,
  #            Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_thresh,
  #            Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_targ), lty=1, col=c(2, "lightgoldenrod2", 3))
  legend("topleft",paste("d)"),bty="n")

}

#' Plot biomass, catches and exploitation (parametric resampling)
#'
#' This function plots relative and absolute biomass, catch and exploitation
#' and associated 95 percent confidence limits (derived from parametric resampling)
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param mod_type 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
#' @param mod_option 1=traditional, 2=chl, 3=dep, 4=chl and dep
#' @param ResampRes object containing outputs from parametric resampling
#' @param param_results object containing outputs from fitted TMB model
#' @param xaxis_lab x axis label, to overide default
#'
#' @return plots of biomass, catch and exploitation
#'
#' @export
Plot_Biomass_Catch_And_Exploitation_Resamp <- function(DatFromCSVFile, mod_type, mod_option, ResampRes, param_results, xaxis_lab) {

  dat=DatFromCSVFile
  year <- ResampRes$model_outputs2$season

  # specify model number
  model_num <- (mod_type-1) * 4 + mod_option

  red = as.numeric(col2rgb(2)) / 255
  orange = as.numeric(col2rgb("orange")) / 255
  yellow = as.numeric(col2rgb("lightgoldenrod2")) / 255
  green = as.numeric(col2rgb(3)) / 255

  par(mfrow=c(2,2), mar=c(3,4,1,2), mgp=c(2, 0.5, 0), oma=c(0,0,2,0))

  # (1) Biomass----
  if (is.na(xaxis_lab)) xaxis_lab = "Year"
  ylims = Get_yaxis_scale(model_outputs$biom_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(year, ResampRes$model_outputs2$biom, type="l", lwd=1, ylim=c(0,ymax),
       xlab=xaxis_lab, ylab="", las=1)
  x <- c(year,rev(year))
  y <- c(ResampRes$model_outputs2$biom_low,rev(ResampRes$model_outputs2$biom_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  mtext("Biomass, t",side=2,line=3,cex=0.7)
  legend("topleft",paste("a)"),bty="n")

  # (2) Relative Biomass----
  B_lim = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_lim
  B_thresh = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_thresh
  B_targ = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$B_targ

  ylims = Get_yaxis_scale(ResampRes$model_outputs2$relbiom_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(year,ResampRes$model_outputs2$relbiom, type="l", ylim=c(0,ymax), lwd=1,
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Relative Biomass", side=2,line=2.5,cex=0.7)
  x <- c(year,rev(year))
  y <- c(ResampRes$model_outputs2$relbiom_low,rev(ResampRes$model_outputs2$relbiom_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  abline(h=c(B_lim, B_thresh, B_targ), lty=1, col=c(2, "lightgoldenrod2", 3))
  legend("topleft",paste("b)"),bty="n")

  ## (3) catch and MSY ----
  ylims = Get_yaxis_scale(ResampRes$model_outputs2$Chat_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(dat$year,dat$tot_catch, type="p", ylim=c(0,ymax), pch=21,
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Catch, t",side=2,line=2.5,cex=0.7)
  lines(ResampRes$model_outputs2$season,ResampRes$model_outputs2$Chat, type="l", lwd=1)
  x <- c(year,rev(year))
  y <- c(ResampRes$model_outputs2$Chat_low,rev(ResampRes$model_outputs2$Chat_upp))
  polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
  abline(h=c(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$MSY), lty=c(1,3,3), col=c(2,2,2))
  legend("topright",paste("MSY ", round(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$MSY1),"t"), bty="n")
  legend("topleft",paste("c)"),bty="n")

  ## (4) exploitation ----
  ylims = Get_yaxis_scale(ResampRes$model_outputs2$Expl_upp)
  ymax = ylims$ymax
  yint = ylims$yint
  plot(ResampRes$model_outputs2$season,ResampRes$model_outputs2$Expl, type="l", lwd=1, ylim=c(0,ymax),
       xlab=xaxis_lab, ylab="", las=1)
  mtext("Exploitation",side=2,line=2.5,cex=0.7)
  x <- c(year,rev(year))
  y <- c(ResampRes$model_outputs2$Expl_low,rev(ResampRes$model_outputs2$Expl_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  # not plotting, as high uncertainty for r estimate may mean these values are misleading
  # abline(h=c(Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_lim,
  #            Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_thresh,
  #            Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$F_targ), lty=1, col=c(2, "lightgoldenrod2", 3))
  legend("topleft",paste("d)"),bty="n")

}


#' Plot observed vs expected CPUE (asymptotic errors)
#'
#' This function plots observed vs expected CPUE and associated 95 percent confidence limits
#' based on asymptotic errors
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param model_outputs object containing outputs from fitted TMB model
#' @param xaxis_lab x axis label, to overide default
#'
#' @return plots of observed vs expected cpue
#'
#' @export
Plot_Obs_vs_Exp_CPUE <- function(DatFromCSVFile, model_outputs, xaxis_lab="Year") {

  dat=DatFromCSVFile
  year <- model_outputs$season
  x <- c(year,rev(year))

  # (5) catch rates----
  par(mfrow=c(2,2), mar=c(3,4,1,2), mgp=c(2, 0.5, 0), oma=c(0,0,2,0))

  cpue1_nyrs = length(which(!is.na(dat$cpue1)))
  cpue2_nyrs = length(which(!is.na(dat$cpue2)))
  cpue3_nyrs = length(which(!is.na(dat$cpue3)))
  cpue4_nyrs = length(which(!is.na(dat$cpue4)))
  first_cpue1_yr = ifelse(cpue1_nyrs>0, min(dat$year[!is.na(dat$cpue1)]), 0)
  last_cpue1_yr = ifelse(cpue1_nyrs>0, max(dat$year[!is.na(dat$cpue1)]), 0)
  first_cpue2_yr = ifelse(cpue2_nyrs>0, min(dat$year[!is.na(dat$cpue2)]), 0)
  last_cpue2_yr = ifelse(cpue2_nyrs>0, max(dat$year[!is.na(dat$cpue2)]), 0)
  first_cpue3_yr = ifelse(cpue3_nyrs>0, min(dat$year[!is.na(dat$cpue3)]), 0)
  last_cpue3_yr = ifelse(cpue3_nyrs>0, max(dat$year[!is.na(dat$cpue3)]), 0)
  first_cpue4_yr = ifelse(cpue4_nyrs>0, min(dat$year[!is.na(dat$cpue4)]), 0)
  last_cpue4_yr = ifelse(cpue4_nyrs>0, max(dat$year[!is.na(dat$cpue4)]), 0)

  if (is.na(xaxis_lab)) xaxis_lab = "Year"
  if(cpue1_nyrs>0){
    CLs = c(dat$obs_ln_cpue1-1.96*dat$cpue1_se, dat$obs_ln_cpue1+1.96*dat$cpue1_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue1, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue1+1.96*dat$cpue1_se,
           dat$year,dat$obs_ln_cpue1-1.96*dat$cpue1_se,angle=90,code=3,length=0)
    x <-c(first_cpue1_yr:last_cpue1_yr,last_cpue1_yr:first_cpue1_yr)
    y <- c(model_outputs$est_cpue1_low[model_outputs$season %in% c(first_cpue1_yr:last_cpue1_yr)],
           rev(model_outputs$est_cpue1_upp[model_outputs$season %in% c(first_cpue1_yr:last_cpue1_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("a"),bty="n")
    legend("topright","CPUE_1",bty="n")
  }

  if(cpue2_nyrs>0){
    CLs = c(dat$obs_ln_cpue2-1.96*dat$cpue2_se, dat$obs_ln_cpue2+1.96*dat$cpue2_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue2, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue2+1.96*dat$cpue2_se,
           dat$year,dat$obs_ln_cpue2-1.96*dat$cpue2_se,angle=90,code=3,length=0)
    x <-c(first_cpue2_yr:last_cpue2_yr,last_cpue2_yr:first_cpue2_yr)
    y <- c(model_outputs$est_cpue2_low[model_outputs$season %in% c(first_cpue2_yr:last_cpue2_yr)],
           rev(model_outputs$est_cpue2_upp[model_outputs$season %in% c(first_cpue2_yr:last_cpue2_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("b"),bty="n")
    legend("topright","CPUE_2",bty="n")
  }

  if(cpue3_nyrs>0){
    CLs = c(dat$obs_ln_cpue3-1.96*dat$cpue3_se, dat$obs_ln_cpue3+1.96*dat$cpue3_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue3, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue3+1.96*dat$cpue3_se,
           dat$year,dat$obs_ln_cpue3-1.96*dat$cpue3_se,angle=90,code=3,length=0)
    x <-c(first_cpue3_yr:last_cpue3_yr,last_cpue3_yr:first_cpue3_yr)
    y <- c(model_outputs$est_cpue3_low[model_outputs$season %in% c(first_cpue3_yr:last_cpue3_yr)],
           rev(model_outputs$est_cpue3_upp[model_outputs$season %in% c(first_cpue3_yr:last_cpue3_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("c"),bty="n")
    legend("topright","CPUE_3",bty="n")
  }

  if(cpue4_nyrs>0){
    CLs = c(dat$obs_ln_cpue4-1.96*dat$cpue4_se, dat$obs_ln_cpue4+1.96*dat$cpue4_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue4, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue4+1.96*dat$cpue4_se,
           dat$year,dat$obs_ln_cpue4-1.96*dat$cpue4_se,angle=90,code=3,length=0)
    x <-c(first_cpue4_yr:last_cpue4_yr,last_cpue4_yr:first_cpue4_yr)
    y <- c(model_outputs$est_cpue4_low[model_outputs$season %in% c(first_cpue4_yr:last_cpue4_yr)],
           rev(model_outputs$est_cpue4_upp[model_outputs$season %in% c(first_cpue4_yr:last_cpue4_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("d"),bty="n")
    legend("topright","CPUE_4",bty="n")
  }
}

#' Plot observed vs expected CPUE (parametric resampling)
#'
#' This function plots observed vs expected CPUE and associated 95 percent confidence limits
#' based on parametric resampling
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param ResampRes object containing outputs from parametric resampling
#' @param xaxis_lab x axis label, to overide default
#'
#' @return plots of observed vs expected cpue
#'
#' @export
Plot_Obs_vs_Exp_CPUE_resamp <- function(DatFromCSVFile, ResampRes, xaxis_lab="Year") {

  dat=DatFromCSVFile
  year <- ResampRes$model_outputs2$season
  x <- c(year,rev(year))

  # (5) catch rates----
  par(mfrow=c(2,2), mar=c(3,4,1,2), mgp=c(2, 0.5, 0), oma=c(0,0,2,0))

  cpue1_nyrs = length(which(!is.na(dat$cpue1)))
  cpue2_nyrs = length(which(!is.na(dat$cpue2)))
  cpue3_nyrs = length(which(!is.na(dat$cpue3)))
  cpue4_nyrs = length(which(!is.na(dat$cpue4)))
  first_cpue1_yr = ifelse(cpue1_nyrs>0, min(dat$year[!is.na(dat$cpue1)]), 0)
  last_cpue1_yr = ifelse(cpue1_nyrs>0, max(dat$year[!is.na(dat$cpue1)]), 0)
  first_cpue2_yr = ifelse(cpue2_nyrs>0, min(dat$year[!is.na(dat$cpue2)]), 0)
  last_cpue2_yr = ifelse(cpue2_nyrs>0, max(dat$year[!is.na(dat$cpue2)]), 0)
  first_cpue3_yr = ifelse(cpue3_nyrs>0, min(dat$year[!is.na(dat$cpue3)]), 0)
  last_cpue3_yr = ifelse(cpue3_nyrs>0, max(dat$year[!is.na(dat$cpue3)]), 0)
  first_cpue4_yr = ifelse(cpue4_nyrs>0, min(dat$year[!is.na(dat$cpue4)]), 0)
  last_cpue4_yr = ifelse(cpue4_nyrs>0, max(dat$year[!is.na(dat$cpue4)]), 0)

  if (is.na(xaxis_lab)) xaxis_lab = "Year"
  if(cpue1_nyrs>0){
    CLs = c(dat$obs_ln_cpue1-1.96*dat$cpue1_se, dat$obs_ln_cpue1+1.96*dat$cpue1_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue1, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue1+1.96*dat$cpue1_se,
           dat$year,dat$obs_ln_cpue1-1.96*dat$cpue1_se,angle=90,code=3,length=0)
    x <-c(first_cpue1_yr:last_cpue1_yr,last_cpue1_yr:first_cpue1_yr)
    y <- c(ResampRes$model_outputs2$est_cpue1_low[ResampRes$model_outputs2$season %in% c(first_cpue1_yr:last_cpue1_yr)],
           rev(ResampRes$model_outputs2$est_cpue1_upp[ResampRes$model_outputs2$season %in% c(first_cpue1_yr:last_cpue1_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("a"),bty="n")
    legend("topright","CPUE_1",bty="n")
  }

  if(cpue2_nyrs>0){
    CLs = c(dat$obs_ln_cpue2-1.96*dat$cpue2_se, dat$obs_ln_cpue2+1.96*dat$cpue2_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue2, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue2+1.96*dat$cpue2_se,
           dat$year,dat$obs_ln_cpue2-1.96*dat$cpue2_se,angle=90,code=3,length=0)
    x <-c(first_cpue2_yr:last_cpue2_yr,last_cpue2_yr:first_cpue2_yr)
    y <- c(ResampRes$model_outputs2$est_cpue2_low[ResampRes$model_outputs2$season %in% c(first_cpue2_yr:last_cpue2_yr)],
           rev(ResampRes$model_outputs2$est_cpue2_upp[ResampRes$model_outputs2$season %in% c(first_cpue2_yr:last_cpue2_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("b"),bty="n")
    legend("topright","CPUE_2",bty="n")
  }

  if(cpue3_nyrs>0){
    CLs = c(dat$obs_ln_cpue3-1.96*dat$cpue3_se, dat$obs_ln_cpue3+1.96*dat$cpue3_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue3, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue3+1.96*dat$cpue3_se,
           dat$year,dat$obs_ln_cpue3-1.96*dat$cpue3_se,angle=90,code=3,length=0)
    x <-c(first_cpue3_yr:last_cpue3_yr,last_cpue3_yr:first_cpue3_yr)
    y <- c(ResampRes$model_outputs2$est_cpue3_low[ResampRes$model_outputs2$season %in% c(first_cpue3_yr:last_cpue3_yr)],
           rev(ResampRes$model_outputs2$est_cpue3_upp[ResampRes$model_outputs2$season %in% c(first_cpue3_yr:last_cpue3_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("c"),bty="n")
    legend("topright","CPUE_3",bty="n")
  }

  if(cpue4_nyrs>0){
    CLs = c(dat$obs_ln_cpue4-1.96*dat$cpue4_se, dat$obs_ln_cpue4+1.96*dat$cpue4_se)
    ylims = Get_yaxis_scale(CLs[!is.na(CLs)])
    y_max = ylims$ymax
    y_min = ylims$ymin
    plot(dat$year,dat$obs_ln_cpue4, type="p", ylim=c(y_min,y_max),
         xlab=xaxis_lab, ylab="", las=1)
    mtext("lnCPUE",side=2,line=2.5,cex=0.7)
    arrows(dat$year,dat$obs_ln_cpue4+1.96*dat$cpue4_se,
           dat$year,dat$obs_ln_cpue4-1.96*dat$cpue4_se,angle=90,code=3,length=0)
    x <-c(first_cpue4_yr:last_cpue4_yr,last_cpue4_yr:first_cpue4_yr)
    y <- c(ResampRes$model_outputs2$est_cpue4_low[ResampRes$model_outputs2$season %in% c(first_cpue4_yr:last_cpue4_yr)],
           rev(ResampRes$model_outputs2$est_cpue4_upp[ResampRes$model_outputs2$season %in% c(first_cpue4_yr:last_cpue4_yr)]))
    polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
    legend("topleft",paste("d"),bty="n")
    legend("topright","CPUE_4",bty="n")
  }
}


#' Plot observed vs expected environmental index
#'
#' This function plots observed vs expected environmental index and associated 95 percent confidence limits
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param model_outputs object containing outputs from fitted TMB model
#' @param xaxis_lab x axis label, to overide default
#' @param y_max y axis maximum, to overide default
#' @param y_min y axis minimum, to overide default
#'
#' @return plots of observed vs expected environmental index
#'
#' @export
Plot_Obs_vs_Exp_Env_Index <- function(DatFromCSVFile, model_outputs, xaxis_lab=NA, y_max=NA, y_min=NA) {

  dat=DatFromCSVFile
  year <- model_outputs$season
  x <- c(year,rev(year))
  if (is.na(xaxis_lab)) xaxis_lab = "Year"

  if (is.na(y_max)) {
    MaxCL = max(c(dat$obs_env+1.96*dat$env_se,model_outputs$est_env_upp))
    ylims = Get_yaxis_scale(MaxCL)
    y_max = ylims$ymax
    y_min = -y_max
  }

  plot(dat$year,dat$obs_env, type="p", ylim=c(y_min,y_max),
       xlab=xaxis_lab, ylab="", las=1)
  lines(dat$year, dat$obs_env, lty="solid")
  lines(dat$year, model_outputs$est_env, lty="dotted")
  mtext("env. index",side=2,line=2.5,cex=0.7)
  arrows(dat$year,dat$obs_env+1.96*dat$env_se,
         dat$year,dat$obs_env-1.96*dat$env_se,angle=90,code=3,length=0)
  x <-c(dat$year,rev(dat$year))
  y <- c(model_outputs$est_env_low, rev(model_outputs$est_env_upp))
  polygon(x,y, col=rgb(0.2,0.2,0.2,0.3), border=NA)
  legend("topright",c("obs. env","exp. env"),bty="n", lty=c("solid","dotted"), pch=c(1,-1))

}

#' Plot exploitation vs relative biomass (i.e. Kobe plot)
#'
#' This function plots exploitation vs relative biomass (i.e. Kobe plot)
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param model_outputs object containing outputs from fitted TMB model
#'
#' @return Kobe plot
#'
#' @export
Plot_Kobe_Plot <- function(DatFromCSVFile, model_outputs) {

  plot(model_outputs$relbiom,model_outputs$Expl, type="l",xlim=c(0,max(model_outputs$relbiom)), ylim=c(0, max(model_outputs$Expl)),
       xlab="Relative Biomass", ylab="", las=1)
  mtext("Exploitation",side=2,line=2.5,cex=1)
  zeroval = -0.05
  uppval <- 1.5

  res=Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)
  B_lim = res$B_lim; B_thresh=res$B_thresh; B_targ=res$B_targ
  F_lim=res$F_lim; F_thresh=res$F_thresh; F_targ=res$F_targ

  red = as.numeric(col2rgb(2)) / 255
  orange = as.numeric(col2rgb("orange")) / 255
  yellow = as.numeric(col2rgb("lightgoldenrod2")) / 255
  green = as.numeric(col2rgb(3)) / 255

  polygon(c(zeroval, B_lim, B_lim, zeroval, zeroval), c(zeroval, zeroval, uppval, uppval, zeroval), col=rgb(red[1], red[2], red[3], 0.3), border=NA,fillOddEven = F)
  polygon(c(B_lim, uppval, uppval, B_lim, B_lim), c(F_lim, F_lim, uppval, uppval, F_lim), col=rgb(red[1], red[2], red[3], 0.3), border=NA)
  polygon(c(B_lim, B_thresh, B_thresh, B_lim, B_lim), c(zeroval, zeroval, F_lim, F_lim, zeroval), col=rgb(orange[1], orange[2], orange[3], 0.3), border=NA)
  polygon(c(B_thresh, uppval, uppval, B_thresh, B_thresh), c(F_lim, F_lim, F_thresh, F_thresh, F_lim), col=rgb(orange[1], orange[2], orange[3], 0.3), border=NA)
  polygon(c(B_thresh, B_targ, B_targ, B_thresh, B_thresh), c(zeroval, zeroval, F_thresh, F_thresh, zeroval), col=rgb(yellow[1], yellow[2], yellow[3], 0.3), border=NA)
  polygon(c(B_targ, uppval, uppval, B_targ, B_targ), c(F_targ, F_targ, F_thresh, F_thresh, F_targ), col=rgb(yellow[1], yellow[2], yellow[3], 0.3), border=NA)
  polygon(c(B_targ, uppval, uppval, B_targ, B_targ), c(zeroval, zeroval, F_targ, F_targ, zeroval), col=rgb(green[1], green[2], green[3], 0.3), border=NA)
  points(model_outputs$relbiom[1],model_outputs$Expl[1], pch=16, cex=1.5, col="green")
  points(model_outputs$relbiom[nyrs],model_outputs$Expl[nyrs], pch=16, cex=1.5, col=2)
  text(model_outputs$relbiom[1],model_outputs$Expl[1],cex=0.8,paste(model_outputs$season[1]),pos=2)
  text(model_outputs$relbiom[nyrs],model_outputs$Expl[nyrs],cex=0.8,paste(model_outputs$season[nyrs]),pos=1)

}

#' Plot estimated random effects
#'
#' This function plots estimated random effects from the state space production model
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param model_outputs object containing outputs from fitted TMB model
#'
#' @return plot of estimated random effects
#'
#' @export
Plot_Estimated_Random_Effects <- function(DatFromCSVFile, model_outputs) {

  dat=DatFromCSVFile
  year <- model_outputs$season
  y_min <- min(model_outputs$EpsR_low) # set y-axis limit)
  y_max <- max(model_outputs$EpsR_upp) # set y-axis limit
  plot(dat$year,model_outputs$EpsR, xlab="Year",ylab="EpsR", ylim=c(y_min,y_max), type="l")
  x <- c(dat$year, rev(dat$year))
  y <- c(model_outputs$EpsR_low,rev(model_outputs$EpsR_upp))
  polygon(x,y,col=rgb(0.2,0.2,0.2,0.2), border=NA)
  abline(h=0, lty="dotted")

}

#' Get variance-covariance matrix, for use in resampling
#'
#' This function outputs variance-covariance matrix, dervied from Hessian matrix outputted by TMB,
#' for use in parametric resampling of model parameters to generate estimates of uncertainty, i.e. without
#' use of asymptotic error approximation. The Params and vcov.Params objects can be passed to the
#' function Get_model_outputs_from_resampling()
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param bdm_params object containing parameters for TMB model
#' @param obj_TMB outputs outputted directly from fitted TMB biomass dynamics model
#' @param model_outputs key model outputs and associated asymptotic variances
#' @param nreps specified number of random values for each parameter produced by parametric resampling
#' @param use60CL_for_Biom flag to output 95 or 60 percent confidence limits. False = 95 percent, True = 60 percent
#'
#' @return Params, vcov.Params
#'
#' @export
Calc_variance_covariance_matrix <- function(obj_TMB) {

  ## select all parameters not including random effect i.e. Eps
  Index = which(names(obj_TMB$env$last.par.best)!="EpsR")
  obj_TMB$env$last.par.best[Index]

  # run optim, to get parameter estimates, along with the Hessian matrix,
  # required to calculate 95% confidence limits - takes a while to run!
  opt <- optim(par=obj_TMB$env$last.par.best[Index],fn=obj_TMB$fn,
               method="BFGS",hessian=TRUE,control=list(trace=1))
  Params <- opt$par

  # get standard errors of the parameter estimates from the hessian matrix, as calculated by optim.
  vcov.Params = NULL
  vcov.Params <- solve(opt$hessian)

  result = list(Params = Params,
                vcov.Params = vcov.Params)

  return(result)

}



#' Output model results from parametric resampling
#'
#' This function outputs model results derived by parametric resampling, i.e. by passing
#' random values of model parameters, from the estimated variance-covariance matrix, through
#' the model and retrieving results, assuming the estimated parameters conform to a multivariate normal distribution
#'
#' @param DatFromCSVFile model data read in from csv file
#' @param bdm_params object containing parameters for TMB model
#' @param obj_TMB outputs outputted directly from fitted TMB biomass dynamics model
#' @param VarCov_res estimated params and variance covariance matrix, from inverted Hessian matrix
#' @param model_outputs key model outputs and associated asymptotic variances
#' @param nreps specified number of random values for each parameter produced by parametric resampling
#' @param use60CL_for_Biom flag to output 95 or 60 percent confidence limits. False = 95 percent, True = 60 percent
#'
#' @return plot of estimated random effects
#'
#' @export
Get_model_outputs_from_resampling <- function(DatFromCSVFile, bdm_params, obj_TMB, VarCov_res, model_outputs, nreps, use60CL_for_Biom) {

  dat=DatFromCSVFile

  # get variance covariance matrix and associated parameters
  if (is.list(VarCov_res)) {
    res =  VarCov_res
  } else {
    res = Calc_variance_covariance_matrix(obj_TMB)
  }
  Params = res$Params
  vcov.Params = res$vcov.Params

  # generate resampled parameter values from multivariate normal distribution using MASS package
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nreps, Params, vcov.Params))
  names(sims) <- names(Params)

  # set up storage
  samp_NLL <- rep(NA, times=nreps)
  season <- c(dat$year,max(dat$year)+1) # years in dataframe + 1
  BiomResults <- data.frame(matrix(ncol=length(season),nrow=nreps))
  colnames(BiomResults) <- season
  prodnResults <- BiomResults
  cpue1Results <- data.frame(matrix(ncol=length(dat$year),nrow=nreps))
  colnames(cpue1Results) <- dat$year
  relResults <- cpue1Results
  cpue2Results <- cpue1Results
  cpue3Results <- cpue1Results
  cpue4Results <- cpue1Results
  FF_Results <- cpue1Results
  Fmort_Results <- cpue1Results
  Expl_Results <- cpue1Results
  Chat_Results <- cpue1Results
  EpsR_Results <- cpue1Results
  vparams <- names(bdm_params[grep("ln_q", names(bdm_params))])
  params <- data.frame(matrix(ncol=length(vparams) + 9,nrow=nreps))
  colnames(params) <- c("NLL", "K","r",vparams,"env_param","lndep","dep",
                        "ln_pt_parm","pt_parm","MSY")


  # do the resampling
  jj <- 0 # row counter
  for (j in 1:nreps) {
    # new parameters based on sims outputs
    if(!is.nan(obj_TMB$fn(sims[j,]))){
      jj <- jj + 1
      obj_TMB$report(obj_TMB$env$last.par)

      sims_results <- obj_TMB$report(obj_TMB$env$last.par)
      BiomResults[jj,] <- sims_results$biom
      relResults[jj,] <- sims_results$relbiom
      prodnResults[jj,] <- sims_results$prodn

      if("ln_q1" %in% names(fit_TMB$par)){
        cpue1Results[jj,] <- sims_results$est_ln_cpue1
        params[jj,'ln_q1'] <- sims_results$ln_q1
      }
      if("ln_q2" %in% names(fit_TMB$par)){
        cpue2Results[jj,] <- sims_results$est_ln_cpue2
        params[jj,'ln_q2'] <- sims_results$ln_q2
      }
      if("ln_q3" %in% names(fit_TMB$par)){
        cpue3Results[jj,] <- sims_results$est_ln_cpue3
        params[jj,'ln_q3'] <- sims_results$ln_q3
      }
      if("ln_q4" %in% names(fit_TMB$par)){
        cpue3Results[jj,] <- sims_results$est_ln_cpue4
        params[jj,'ln_q4'] <- sims_results$ln_q4
      }

      Fmort_Results[jj,] <- sims_results$Fmort
      Expl_Results[jj,] <- sims_results$Expl
      FF_Results[jj,] <- sims_results$FF
      Chat_Results[jj,] <- sims_results$Chat
      EpsR_Results[jj,] <- sims_results$EpsR

      # parameter outputs
      params[jj,'NLL'] <- sims_results$NLL
      params[jj,'K'] <- sims_results$K
      params[jj,"r"] <- sims_results$r

      if("env_p" %in% names(fit_TMB$par)){
        params[jj,"env_param"] <- sims_results$env_param
      }
      if("lndep" %in% names(fit_TMB$par)){
        params[jj,"lndep"] <- sims_results$lndep
        params[jj,"dep"] <- sims_results$dep
      }
      if("ln_pt_parm"%in%names(fit_TMB$par)){ # if using PT model
        params[jj,"ln_pt_parm"] <- sims_results$ln_pt_parm
        params[jj,"pt_parm"] <- sims_results$pt_parm
      }

      # MSY
      param_results = list(r = params[jj,"r"],
                           K = params[jj,'K'],
                           pt_parm = params[jj,"pt_parm"],
                           env_param = params[jj,"env_param"],
                           dep = params[jj,"dep"])
      params[jj,"MSY"] <- Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env=0)$MSY

      cat(j, "\n")
    }
  }

  # ****************************************
  # get uncertainty for estimated parameters
  # ****************************************

  param_est <- data.frame(matrix(ncol=24, nrow=1))
  names(param_est) <- c("K_mn","K_med",'K_low','K_upp',
                        'r_mn','r_med','r_low','r_upp',
                        "env_mn",'env_med',"env_low","env_upp",
                        "dep_mn",'dep_med',"dep_low","dep_upp",
                        "pt_p_mn",'pt_p_med',"pt_p_low","pt_p_upp",
                        "MSY_mn", "MSY_med","MSY_low","MSY_upp")

  param_est$K_mn <- mean(na.omit(params$K))
  param_est$K_med <- quantile(na.omit(params$K),0.5)
  param_est$K_low <- quantile(na.omit(params$K),0.025)
  param_est$K_upp<- quantile(na.omit(params$K),0.975)
  param_est$r_mn <- mean(na.omit(params$r))
  param_est$r_med <- quantile(na.omit(params$r),0.5)
  param_est$r_low <- quantile(na.omit(params$r),0.025)
  param_est$r_upp <- quantile(na.omit(params$r),0.975)
  param_est$env_mn <- NA
  param_est$env_med <- NA
  param_est$env_low <- NA
  param_est$env_upp <- NA
  param_est$dep_mn <- NA
  param_est$dep_med <- NA
  param_est$dep_low <- NA
  param_est$dep_upp <- NA
  param_est$pt_p_mn <- NA
  param_est$pt_p_med <- NA
  param_est$pt_p_low <- NA
  param_est$pt_p_upp <- NA
  param_est$MSY_mn <-   mean(na.omit(params$MSY))
  param_est$MSY_med <- quantile(na.omit(params$MSY),0.5)
  param_est$MSY_low <- quantile(na.omit(params$MSY),0.025)
  param_est$MSY_upp <- quantile(na.omit(params$MSY),0.975)

  if("env_p" %in% names(fit_TMB$par)){
    param_est$ln_env_mn <- mean(na.omit(params$env_p))
    param_est$ln_env_med <- quantile(na.omit(params$env_p),0.5)
    param_est$ln_env_low <- quantile(na.omit(params$env_p),0.025)
    param_est$ln_env_upp <- quantile(na.omit(params$env_p),0.975)
  }
  if("lndep" %in% names(fit_TMB$par)){
    param_est$dep_mn <- mean(na.omit(params$lndep))
    param_est$dep_med <- quantile(na.omit(params$lndep),0.5)
    param_est$dep_low <- quantile(na.omit(params$lndep),0.025)
    param_est$dep_upp <- quantile(na.omit(params$lndep),0.975)
  }
  if("ln_pt_parm" %in% names(fit_TMB$par)){
    param_est$ln_pt_p_mn <- mean(na.omit(params$ln_pt_parm))
    param_est$ln_pt_p_med <- quantile(na.omit(params$ln_pt_parm),0.5)
    param_est$ln_pt_p_low <- quantile(na.omit(params$ln_pt_parm),0.025)
    param_est$ln_pt_p_upp <- quantile(na.omit(params$ln_pt_parm),0.975)
  }

  # ***************************************************
  # get estimates of uncertainty for dervied quantities
  # ***************************************************
  # can get 95 or 60 confidence intervals
  if (use60CL_for_Biom){ low_pc = 0.2; upp_pc = 0.8 } else { low_pc = 0.025; upp_pc = 0.975 }

  # create a copy dataframe of the model outputs
  model_outputs2 <- model_outputs

  # biomass
  model_outputs2$biom = apply(BiomResults[, 1:nyrs], MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$biom_low = apply(BiomResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, low_pc, na.rm=TRUE))
  model_outputs2$biom_upp = apply(BiomResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, upp_pc, na.rm=TRUE))

  # relative biomass
  model_outputs2$relbiom = apply(relResults[, 1:nyrs], MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$relbiom_low = apply(relResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, low_pc, na.rm=TRUE))
  model_outputs2$relbiom_upp = apply(relResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, upp_pc, na.rm=TRUE))

  # production
  model_outputs2$prodn = apply(prodnResults[, 1:nyrs], MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$prodn_low = apply(prodnResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
  model_outputs2$prodn_upp = apply(prodnResults[, 1:nyrs], MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))

  # cpue - these are estiamtes from resampling
  if("ln_q1"%in%names(fit_TMB$par)){
    model_outputs2$cpue1 = apply(cpue1Results, MARGIN=2, FUN=mean, na.rm=TRUE)
    model_outputs2$cpue1_low = apply(cpue1Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
    model_outputs2$cpue1_upp = apply(cpue1Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))
  }
  if("ln_q2"%in%names(fit_TMB$par)){
    model_outputs2$cpue2 = apply(cpue2Results, MARGIN=2, FUN=mean, na.rm=TRUE)
    model_outputs2$cpue2_low = apply(cpue2Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
    model_outputs2$cpue2_upp = apply(cpue2Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))
  }

  if("ln_q3"%in%names(fit_TMB$par)){
    model_outputs2$cpue3 = apply(cpue3Results, MARGIN=2, FUN=mean, na.rm=TRUE)
    model_outputs2$cpue3_low = apply(cpue3Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
    model_outputs2$cpue3_upp = apply(cpue3Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))
  }

  if("ln_q4"%in%names(fit_TMB$par)){
    model_outputs2$cpue4 = apply(cpue4Results, MARGIN=2, FUN=mean, na.rm=TRUE)
    model_outputs2$cpue4_low = apply(cpue4Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
    model_outputs2$cpue4_upp = apply(cpue4Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))
  }
  # exploitation rate
  model_outputs2$Expl = apply(Expl_Results, MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$Expl_low = apply(Expl_Results, MARGIN=2, FUN=function(x) quantile(x, low_pc, na.rm=TRUE))
  model_outputs2$Expl_upp = apply(Expl_Results, MARGIN=2, FUN=function(x) quantile(x, upp_pc, na.rm=TRUE))

  # EpsR
  model_outputs2$EpsR = apply(EpsR_Results, MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$EpsR_low = apply(EpsR_Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
  model_outputs2$EpsR_upp = apply(EpsR_Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))

  # expected catch
  model_outputs2$Chat = apply(Chat_Results, MARGIN=2, FUN=mean, na.rm=TRUE)
  model_outputs2$Chat_low = apply(Chat_Results, MARGIN=2, FUN=function(x) quantile(x, 0.025, na.rm=TRUE))
  model_outputs2$Chat_upp = apply(Chat_Results, MARGIN=2, FUN=function(x) quantile(x, 0.975, na.rm=TRUE))


  result = list(params = params,
                param_est = param_est,
                model_outputs2 = model_outputs2)

  return(result)

}

