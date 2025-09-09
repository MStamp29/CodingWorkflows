####################################################
####################################################
## This script simulates data based on values from the  
## file created in CreateDF_for_Trending_Simulation.R, 
## and then applies the trending model to test 
## whether the model recovers the input parameters.
####################################################
####################################################

library(mvtnorm)
library(MCMCglmm)

#Sets the wd to source file location on RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("1.PhenMetric_functions.R")

#Clear environment apart from functions 
rm(list = setdiff(ls(), lsf.str()))

simulation_inputs <- read.csv("TrendingSim_parameter_sets.csv")

intercept <- 70

#Model settings - iterations, prior constant 
nittN <- 300000
a <- 1000 # For prior

#Run the simulations
Trending_sim(param_df = simulation_inputs, 
             rowN = 1, #The row of simulation_inputs to use
             reps = 1, #The number of times to simulate the data and run the model
             nittN = nittN)


