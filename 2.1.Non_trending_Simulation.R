####################################################
####################################################
## This script simulates data based on values in 
## Non_trendingSim_parameter_sets.csv, 
## and then applies the non-trending model to test 
## whether the model recovers the input parameters.
####################################################
####################################################

#Sets the wd to source file location on RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MCMCglmm)
#library(tidyverse)

source('1.PhenMetric_functions.R')

#Clear environment apart from functions 
rm(list = setdiff(ls(), lsf.str()))

## simulation_inputs is a df where each row are the sample sizes and variances used to simulate data
simulation_inputs <- read.csv("Non_trendingSim_parameter_sets.csv")

#Model settings - iterations, prior constant 
nittN <- 250000 #160000 for scenario 4
a <- 1000 # For prior

#Run the simulations
Non_trending_sim(param_df = simulation_inputs, 
                 rowN = 1, #The row of simulation_inputs to use
                 reps = 1, #The number of times to simulate the data and run the model
                 nittN = nittN)




