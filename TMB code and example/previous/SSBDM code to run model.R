
# L4Assess: State space biomass dynamics model incorporating environmental effects.
# Stock Assessment team, DPIRD (2023)
# Rachel Marks, Alex Hesp, Norm Hall & Ainslie Denham
# Lasted updated 29/11/2023 (allow estimation of pInit - initial relative biomass)

# Note the this model is broadly based on the model described by Marks et al. (2021), but
# is NOT exactly the same, with several modifications. This model is still
# in development, with this implementation closely reflecting the state space biomass dynamics
# models currently in use by the assessment team. 

# clear working directory
rm(list=ls())
# install.packages("C:/~/L4Assess_0.1.0.tar.gz", source = TRUE, repos=NULL)

library(TMB)
library(L4Assess)
citation("L4Assess")

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# compile TMB model----
compile("SSBDM.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("SSBDM"))

# read in dat csv file to get rest of data
DatFromCSVFile = read.csv("SSBDM CS crab data.csv", header=T)
names(DatFromCSVFile)

# Specify the model type and structure to be fitted
mod_type = 1       # 1=Schaefer, 2=Fox, 3=Pella_Tomlinson
mod_option = 2     # 1=traditional, 2=chl, 3=dep, 4=chl and dep
mod_scenario = 2   # 1=non-state space, 2=state space
nyrs = length(DatFromCSVFile$season)
wt_param_pen = 1.0 # parameter penalty
wt_depl_pen = 1000.0 # depletion penalty  # set to zero for most or all stocks except CS crabs
wt_biom_pen = 0.1 # biomass penalty
wt_harv_pen = 100.0 # harvest rate penalty
wt_cpue1 = 1 
wt_cpue2 = 1
wt_cpue3 = 1
wt_cpue4 = 1
SigmaR = 0.2 # assumed level of process error
Sigma_env = 0.1 # assumed strength of environmental link with biomass (low value = strong relationship) 
max_currBrel = 1.0 # upper limit to relative biomass in final year, i.e. final depletion penalty. Turn off by setting to high value.
# specify lower and upper bounds for parameters
pInit_bnds = c(0,1)
ln_K_bnds = c(6,9)
ln_r_bnds = log(c(0.01,1.5))
ln_q_bnds = c(-20,0)
ln_sd_bnds = c(-20,1)
env_param_bnds = c(-20,20)
ln_dep_bnds = c(-20,0)
ln_pt_bnds = c(-20,0)
bdm_param_bounds = Set_BoundsForBDM_Params(ln_K_bnds, ln_r_bnds, ln_q_bnds, ln_sd_bnds,
                                           env_param_bnds, ln_dep_bnds, ln_pt_bnds, pInit_bnds)
length(bdm_param_bounds$low_bound_list[14])
length(bdm_param_bounds$upp_bound_list[14])

# get data inputs for TMB model
bdm_data = Get_BDM_Data(DatFromCSVFile, wt_param_pen, wt_depl_pen, wt_biom_pen, wt_harv_pen, 
                        wt_cpue1, wt_cpue2, wt_cpue3, wt_cpue4, SigmaR, Sigma_env, max_currBrel)

# specify initial parameters for TMB model
pInit = 0.5
ln_K = log(2000)
ln_r = log(0.5) 
ln_q1 = -10  
ln_q2 = -10 
ln_q3 = -10
ln_q4 = -10
ln_sd1 = log(0.02) # additional 'model' error associated with cpue
ln_sd2 = log(0.02)
ln_sd3 = log(0.02)
ln_sd4 = log(0.02)
env_p = 0 # 0 = normal environmental conditions
lndep = log(0.5)  # -10 # 1 = full depensation, close to zero = no depensation (can't be zero)
ln_pt_parm = 1 # 1 is equivalent to Schaefer model
nyrs = length(DatFromCSVFile$year)
FF = rep(-2, nyrs)
EpsR = rep(0, times=nyrs)

# get estimated parameter list
bdm_params = Set_InitValsForBDM_Params(pInit, ln_K, ln_r, ln_q1, ln_q2, ln_q3, ln_q4, ln_sd1, ln_sd2, 
                                       ln_sd3, ln_sd4, env_p, lndep, ln_pt_parm, FF, EpsR)
bdm_params


# set parameter map for TMB model (i.e. list of which model parameters not to estimate)
# setting fix_pInit=1 results in initial relative biomass being fixed at user-specified input value,
# setting to 1, it is estimated.
bdm_map = Get_BDM_Map(DatFromCSVFile, mod_scenario, mod_option, fix_pInit=0)

# fit the model
result = fit_the_model(DatFromCSVFile, mod_scenario, mod_type, mod_option, bdm_data, bdm_params, fix_pInit=0, bdm_map = NULL, initial_params = NULL)
result$fit_TMB$par
# get results
fit_TMB = result$fit_TMB
fit_TMB$par
fit_TMB$objective
fit_TMB$convergence

# Collect results using the report function in TMB
result_TMB = result$rep_TMB
result_TMB$r
result_TMB$K
result_TMB$env_param
result_TMB$pInit

# Collect the BDM_model used when fitting using wrapper function
obj_TMB = result$obj_TMB

# Store key model outputs and associated asymptotic variances
model_outputs = Get_model_outputs(result, nyrs)

# Store model parameter estimates and associated asymptotic variances
param_results = Get_model_parameter_outputs(result, mod_type, mod_option)
param_results

# plot observed vs expected cpue
Plot_Obs_vs_Exp_CPUE(DatFromCSVFile, model_outputs, xaxis_lab=NA)

# plot biomass, catch and exploitation
Plot_Biomass_Catch_And_Exploitation(DatFromCSVFile, mod_type, mod_option, model_outputs, param_results, xaxis_lab=NA)

# plot observed vs expected environmental index
Plot_Obs_vs_Exp_Env_Index(DatFromCSVFile, model_outputs, xaxis_lab=NA, y_max=NA, y_min=NA)

# produce Kobe plot
Plot_Kobe_Plot(DatFromCSVFile, model_outputs)

# plot random effects (linked to environment index)
if (mod_scenario == 2) {
  Plot_Estimated_Random_Effects(DatFromCSVFile, model_outputs)  
}


# get deterministic estimates for MSY and Bmsy (not assuming regime shift)
env = 0 # 'average' env cond. # env = -0.5 # below average env cond.
MSYRes = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env)
MSYRes$MSY
MSYRes$Bmsy

env = -0.5 # 'average' env cond. # env = -0.5 # below average env cond. (not assuming regime shift)
MSYRes = Calc_MSYAndBiolRefPoints(mod_type, mod_option, param_results, env)
MSYRes$MSY
MSYRes$Bmsy

# mod_option = 2 (Schaefer with environmental effect)
# results for model 2
result_TMB$r
result_TMB$K


# ********************************************
# For Simon - dynamic reference point question
# ********************************************

# Current environment at env = -0.5
Orig_K = result_TMB$K
r = result_TMB$r
env_param = result_TMB$env_param
env = -0.5
MSY = MSYRes$MSY

# calculate new carrying capacity
New_K = Orig_K * exp(env_param * env)
New_K

# or...
# New_K check
# MSY = rK/4, thus K=4MSY/r
New_K_check = 4*MSY/r
New_K_check

# calculate new BMSY (i.e. after assumed regime shift)
New_BMSY = New_K / 2
New_BMSY

# so if fishing to threshold, MSY is still 29.55443
# if fishing to a target of 1.2BMSY, long term average equilibrium catch = long term average production at 1.2BMSY
# following Schaefer model...
New_B_targ = 1.2 * New_BMSY
New_targ_catch = r * New_B_targ * (1 - New_B_targ / New_K_check)
New_targ_catch

# Note, the above would not apply if depensation were assumed NOT to occur.
# In the paper by Marks et al, the "best" fitting model was that allowing for both environmental effects 
# and depensatory effects. 

# End of productivity change question
# ----------------------------


# *******************************************************************
# Get estimates of model outputs using parametric resampling approach
# *******************************************************************

# get variance-covariance matrix. Note, this takes a minute or two!
VarCov_res = Calc_variance_covariance_matrix(obj_TMB)

# get model estimates and associated uncertainty from parametric resampling  
nreps = 200
ResampRes = Get_model_outputs_from_resampling(DatFromCSVFile, bdm_params, obj_TMB, VarCov_res, model_outputs, nreps, use60CL_for_Biom=F)
ResampRes$param_est

# plot observed vs expected cpue
Plot_Obs_vs_Exp_CPUE_resamp(DatFromCSVFile, ResampRes, xaxis_lab=NA)

# get estimates of MSY and BMSY, for specified environmental condition, with uncertainty, from resampling
Get_MSY_and_BMSY_with_err(mod_type, mod_option, ResampRes, env)

# plot biomass, catch and exploitation
Plot_Biomass_Catch_And_Exploitation_Resamp(DatFromCSVFile, mod_type, mod_option, ResampRes, param_results, xaxis_lab=NA)


