########################################################################
########################################################################
# This script calculates the real P_Corr.Tr and P_Div.Tr
# based on the variances used to simulate data in 'Trending_Simulation.R
########################################################################
########################################################################

library(mvtnorm)

#Sets the wd to source file location on RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())

source("1.PhenMetric_functions.R")
####

intercept <- 70

## 1. Create scenarios (i.e. sample sizes) and parameters -----

## Sample sizes for species number, year number and number of individuals
medSpeciesN <- 10
highSpeciesN <- 30

medYearN <- 10
highYearN <- 30

medIndN <- 5
highIndN <- 30

speciesN_scen <- c(highSpeciesN, highSpeciesN, medSpeciesN, highSpeciesN)
yearN_scen <- c(highYearN, medYearN, highYearN, highYearN)
indN_scen <- c(medIndN, medIndN, medIndN, highIndN)
scenarios <- cbind(speciesN_scen, 
                   yearN_scen, 
                   indN_scen)

## Parameter inputs
fixedYear_slope <- 0

fixedYear_slope_var <- 0.1
medYear_slope_var <- 0.2
highYear_slope_var <- 0.3

fixedYear_int_var <- 10
fixed_slopeint_cov <- 0
fixedResidual_var <- 10
fixedYearspecies_var <- 10

fixedResidual_slope <- 1
medResidual_slope <- 1.25 
highResidual_slope <- 1.5 

fixedYear_var <- 10
medYear_var <- 20
highYear_var <- 40

year_slope_param <- rep(fixedYear_slope, 9)
year_slope_var_param <- c(rep(fixedYear_slope_var, 3), 0, medYear_slope_var, highYear_slope_var, rep(fixedYear_slope_var, 3))
year_int_var_param <- rep(fixedYear_int_var, 9)
slopeint_cov_param <- rep(fixed_slopeint_cov, 9)
residual_var_param <- rep(fixedResidual_var, 9)
resid_slope_param <- c(rep(fixedResidual_slope, 6), 0, medResidual_slope, highResidual_slope)
yearspecies_var_param <- rep(fixedYearspecies_var, 9)
year_var_param <- c(0, medYear_var, highYear_var, rep(fixedYear_var, 6))


parameter_sets <- cbind(year_slope_param,
                        year_slope_var_param,
                        year_int_var_param,
                        slopeint_cov_param,
                        residual_var_param,
                        resid_slope_param, 
                        yearspecies_var_param,
                        year_var_param)


simulation_inputs <- as.data.frame(cbind(scenario = rep(1:nrow(scenarios), each = nrow(parameter_sets)),
                                  parameter_space = rep(1:nrow(parameter_sets), times = nrow(scenarios)),
                                  scenarios[rep(1:nrow(scenarios), each = nrow(parameter_sets)), ],
                                  parameter_sets[rep(1:nrow(parameter_sets), times = nrow(scenarios)), ]))


## 2. Calculate the P_Corr.tr for each parameter set ------
#Set speciesN to 1000000 --> allows P_Corr.tr to asymptote to its true value
simulation_inputs$P_Corr.tr_real <- NA

for(i in 1:nrow(simulation_inputs)) {
  
  speciesN <- 1000000 
  yearN <- simulation_inputs$yearN_scen[i]
  years <- c(1:yearN)-1
  year_slope <- simulation_inputs$year_slope_param[i] 
  speciesIntercept_var <- simulation_inputs$year_int_var_param[i]
  covar <- simulation_inputs$slopeint_cov_param[i]
  speciesSlope_var <- simulation_inputs$year_slope_var_param[i]
  year_var <- simulation_inputs$year_var_param[i]
  yearspecies_var <- simulation_inputs$yearspecies_var_param[i]
  
  #Variance due to the year fixed effect
  year_effect_var <- year_slope^2*var(years) 
  
  #Get the intercept and slope for each species
  species_IntsSlopes <- rmvnorm(speciesN, 
                                c(intercept, year_slope), #mean for the distributions from which to draw from (intercept and slope)
                                sigma = matrix(nrow = 2, c(speciesIntercept_var, covar, covar, speciesSlope_var)))
  
  
  #Variance due to the species-specific year slopes 
  spyear_effect_var <- rep(NA, length(speciesN))
  for(l in 1:speciesN) {
    spyear_effect_var[l] <- species_IntsSlopes[l, 2]^2*var(years)
  }
  
  #Calculate the P_Corr.tr metric for each species
  P_Corr.tr_realSp <- c()
  for(l in 1:speciesN) {
    P_Corr.tr_realSp[l] <- (year_var + year_effect_var)/(year_var + year_effect_var + yearspecies_var + spyear_effect_var[l])
  }
 
  # Get the mean value across species
  simulation_inputs$P_Corr.tr_real[i] <- mean(P_Corr.tr_realSp)
  
  }


# 3. Calculate the P_Div.tr for each parameter set ----------
colnames_Year <- c(0:(max(yearN_scen)-1))
P_Div.tr_colnames <- paste0("P_Div.tr", "Year", colnames_Year, "_", "input")
simulation_inputs[, P_Div.tr_colnames] <- NA

P_Div.tr_real <- list()

for(i in 1:nrow(simulation_inputs)) {
  speciesIntercept_var <- simulation_inputs$year_int_var_param[i]
  covar <- simulation_inputs$slopeint_cov_param[i]
  speciesSlope_var <- simulation_inputs$year_slope_var_param[i]
  residual_var <- simulation_inputs$residual_var_param[i]
  residvar_slope <- simulation_inputs$resid_slope_param[i]
  yearN <- simulation_inputs$yearN_scen[i]
  years <- c(1:yearN)-1

  sigma_real <- matrix(nrow = 2, ncol = 2,
                       c(speciesIntercept_var, covar, covar, speciesSlope_var))
  
  residvar_add <- rep(NA, yearN)
  specvar_real <- rep(NA, yearN)
  P_Div.tr_real[[i]] <- rep(NA, yearN)
  
  for(l in 1:yearN) {
    residvar_add[l] <- years[l]*residvar_slope
    specvar_real[l] <- getspecvar(sigma_real, as.numeric(as.character(years[l])))
    P_Div.tr_real[[i]][l] <- specvar_real[l]/(specvar_real[l] + residual_var + residvar_add[l])
    simulation_inputs[i, grep(P_Div.tr_colnames[1], colnames(simulation_inputs)):grep(P_Div.tr_colnames[yearN], colnames(simulation_inputs))] <-
      P_Div.tr_real[[i]]
    
  }

}

write.csv(simulation_inputs, "TrendingSim_parameter_sets.csv", row.names = F)


## simulation_inputs is the input file for Trending_Simulation.R