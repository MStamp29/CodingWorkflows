ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
GNU General Public Licence
2. Data derived from other sources
Temperature data used is publicly available at https://harvardforest1.fas.harvard.edu/
exist/apps/datasets/showData.html?id=hf001 (Fisher Meteorological Station data) and  https://harvardforest1.fas.harvard.
edu/exist/apps/datasets/showData.html?id=hf000 (Shaler Meteorological Station data)
3. Recommended citation for this data/code archive
Stamp, M. A. 2025. R Code to accompany "A model-based approach for examining multi-species phenological correlation and divergence". 

DATA & CODE FILE OVERVIEW
This data repository consist of 3 data files, 7 code scripts, and this README document, with the following data and code filenames and variables

Data files and variables
    1. [Harvard_springphen.csv] [year = the year of observation; tree.id = individual-level ID comprised of the first two letters of the genus plus the first two letters of the species followed by a two digit number; species = species code; bb_doy = day of year (ordinal day) on which 50% of the buds were open with visible leaves]
    2. [Harvard_temperature_1990to2023.csv] [date = date of temperature reading; airt = average air temperature (unit: Celsius); ord_day = ordinal day of year; year = year]
    3. [Non_trendingSim_parameter_sets.csv] [scenario = identifies which set of sample sizes (speciesN, yearn, indN) are used; parameter_space = identifies which set of input variances are used; speciesN_scen
 = number of species to be simulated; yearN_scen = number of species to be simulated; indN_scen = number of individuals per species to be simulated; species_var_param = species variance used to simulate data; residual_var_param = residual variance for each species used to simulate data; yearspecies_var_param = species-specific year variance for each species used to simulate data; year_vara_param = year variance used to simulate data]
    


Code scripts and workflow
    1. 1.PhenMetric_functions.R Contains all functions needed to run the rest of the scripts in this repository. This includes functions to calculate the correlation and divergence metrics, and functions to simulate data. This script is sourced in all other scripts in this repository apart from 3.1.Harvard_Temperature_Regression.R. 

Model performance scripts - 2.2.CreateDF_for_Trending_Simulation.R needs to be run before 2.3.Trending.Simulation.R:
    2. 2.1.Non_trending_Simulation.R This script simulates data based on values in Non_trendingSim_parameter_sets.csv, and then applies the non-trending model to test whether the model recovers the input parameters. 
    3. 2.2.CreateDF_for_Trending_Simulation.R This script creates a dataframe similar to Non_trendingSim_parameter_sets.csv but for inputs in to the trending simulation. The key element of this script is the calculation of the 'real' value of P_Corr.Tr given a set of input variances - we calculate this for a speciesN of 1000000 which allows P_Corr.Tr to asymptote to its true value.
    4. 2.3.Trending.Simulation.R This script simulates data based on values from the file created in CreateDF_for_Trending_Simulation.R, and then applies the trending model to test whether the model recovers the input parameters. 

Case study scripts - 3.2.Harvard.R needs to be run before 3.3.Harvard_OverlapChange.R
    5. 3.1.Harvard_Temperature_Regression This models the change in temperature around the Harvard Forest area over time for the spring months. 
    6. 3.2.Harvard.R This runs both the non-trending and trending models on the Harvard Forest dataset (Harvard_springphen.csv) and calculates the correlation and divergence metrics. 
    7. 3.3.Harvard_OverlapChange.R Calculates a) the overlap in phenological distributions (using estimates of the mean and standard deviation for each species from the trending model in Harvard.R) in the first year of the Harvard Forest dataset, and b) the change in overlap between the first and last year. 


SOFTWARE VERSIONS
R v4.4.1

Loaded packages:
MCMCglmm v2.36
mvtnorm v1.2-5
lme4 v1.1-35.4
lmerTest v3.1-3



