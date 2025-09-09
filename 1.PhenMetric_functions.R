##############################################################################
##############################################################################
##Contains all functions needed to run the rest of the scripts in this project. 
##This includes functions to calculate the correlation and divergence metrics, 
##and functions to simulate data.
##############################################################################
##############################################################################


## Metric calculations -----
##For non-trending models

##Calculate P_Corr
#model = the model object (from the non-trending model)
#species_names = a vector of the species names as they appear in the dataset used to run the model
#year_pred = the name of the year random effect variable as a string
#species_pred = the name of the species random effect variable as a string
P_Corr <- function(model, species_names, year_pred, species_pred) { 
  
  P_Corr_post_sp <- list()
  
  for(l in 1:nrow(model$VCV)) {
    P_Corr_post_sp[[l]] <- c(rep(NA, times = length(species_names)))
    for(m in 1:length(species_names)) {
      P_Corr_post_sp[[l]][m] <- model$VCV[l, year_pred]/(model$VCV[l, year_pred] + model$VCV[l, paste0(species_pred, species_names[m], ".", year_pred)])
    }
  }
  
  P_Corr_post <- unlist(lapply(P_Corr_post_sp, mean))
  
  return(P_Corr_post)
}

##Calculate P_Div
#model = the model object (from the non-trending model)
#species_names = a vector of the species names as they appear in the dataset used to run the model
#species_pred = the name of the species random effect variable as a string
P_Div <- function(model, species_names, species_pred) {
  
  P_Div_post_sp <- list()
  
  for(l in 1:nrow(model$VCV)) {
    P_Div_post_sp[[l]] <- rep(NA, times = length(species_names))
    for(m in 1:length(species_names)) {
      P_Div_post_sp[[l]][m] <- model$VCV[l, species_pred]/(model$VCV[l, species_pred] +  model$VCV[l, paste0(species_pred, species_names[m], ".units")])
    }
  }
  
  P_Div_post <- unlist(lapply(P_Div_post_sp, mean))
  
  return(P_Div_post)
  
}


##For trending models 
#Get the species variance at a point Xval from the sigma matrix
#sigma = a variance covariance matrix consisting of the variance in intercepts, 
#the covariance between intercept and slopes, and the variance in slopes
#xval = any value from the continuous predictor
getspecvar <- function(sigma, xval) {
  specvar <- sigma[1, 1] + sigma[2, 2]*xval^2 + 2*xval*sigma[2, 1]
}

#Calculate P_Div.tr
#model = the model object (from the trending model)
#driver_vector = a vector of unique values of the continuous predictor from the model
#species_names = a vector of the species names as they appear in the dataset used to run the model
#species_pred = the name of the species random effect variable as a string
#year_pred = the name of the year random effect variable as a string
getP_Div.tr <- function(model, driver_vector, species_names, species_pred, year_pred) {
  
  species_number <- length(species_names)
  
  sigma_post <- list()
  specvar_post <- list()
  residvar_speciesest_post <- list()
  P_Div.tr_species_post <- list()
  P_Div.tr_post <- list()
  residvar_est_post <- list()
  
  sigma_colnames <- c(paste0("(Intercept):(Intercept).", species_pred),
                      paste0(year_pred, ":(Intercept).", species_pred),
                      paste0("(Intercept):", year_pred, ".", species_pred),
                      paste0(year_pred, ":", year_pred, ".", species_pred))
  
  for(i in 1:nrow(model$VCV)) { #
    sigma_post[[i]] <- matrix(nrow = 2, ncol = 2,
                              c(model$VCV[i, sigma_colnames]))
    specvar_post[[i]] <- rep(NA, length(driver_vector))
    residvar_speciesest_post[[i]] <- list()
    P_Div.tr_species_post[[i]] <- list()
    P_Div.tr_post[[i]] <- rep(NA, length(driver_vector))
    residvar_est_post[[i]] <- rep(NA, length(driver_vector))
    
    for(j in 1:length(driver_vector)) {
      specvar_post[[i]][j] <- getspecvar(sigma_post[[i]], driver_vector[j])
      residvar_speciesest_post[[i]][[j]] <- rep(NA, species_number)
      P_Div.tr_species_post[[i]][[j]] <- rep(NA, species_number)
      
      for(k in 1:species_number) {
        #Get the estimated variance for the species in question at each year
        #This consists of the residual variance term + the additional variance estimated by the variance slope term
        residvar_speciesest_post[[i]][[j]][k] <- model$VCV[i, paste(species_pred, species_names[k], ".units", sep = "")] + #speciespredX.units
          model$VCV[i, grep('sqrt', colnames(model$VCV))]*driver_vector[j] #sqrt(yearpred).units * sqrt(year) 
        
        #P_Div.tr for each species in each year - speciesVariance/(speciesVariance + residualVariance)
        P_Div.tr_species_post[[i]][[j]][k] <- specvar_post[[i]][j]/(specvar_post[[i]][j] + residvar_speciesest_post[[i]][[j]][k])###
        
      }
      #P_Div.tr for each year - the mean of the ICC spread of each species in that year
      P_Div.tr_post[[i]][j] <- mean(P_Div.tr_species_post[[i]][[j]])
      residvar_est_post[[i]][[j]] <- mean(residvar_speciesest_post[[i]][[j]])
    }
  }
  
  metrics <- list('Species var' = specvar_post, 'P_Div.tr' = P_Div.tr_post, 'Resid var' = residvar_est_post)
  return(metrics)
  
}

#Calculate P_Corr.tr
#model = the model object (from the trending model)
#driver_name = the name of the continuous predictor variable as a string (likely to be the same as 'year_pred' but not necessarily)
#driver_vector = a vector of unique values of the continuous predictor from the model
#species_names = a vector of the species names as they appear in the dataset used to run the model
#species_pred = the name of the species random effect variable as a string
#year_pred = the name of the year random effect variable as a string
getP_Corr.tr <- function(model, driver_name, driver_vector, year_pred, species_names, species_pred) {
  var_pred_values_post <- c()
  sp_var_pred_values_post <- list() #Get variance in species year effects for each species for the posterior
  P_Corr.tr_species_post <- list()
  P_Corr.tr_post <- c()
  
  speciesN <- length(species_names)
  
  for(i in 1:nrow(model$Sol)) {#nrow(model$Sol)
    #Get the variance in the values predicted by the estimated year slope across years
    var_pred_values_post[i] <- model$Sol[i, driver_name]^2*var(driver_vector) #var(model$Sol[i, "(Intercept)"] + model$Sol[i, driver_name]*driver_vector) 
    sp_var_pred_values_post[[i]] <- c(rep(NA, times = speciesN))
    P_Corr.tr_species_post[[i]] <- rep(NA, times = speciesN)
    
    for(j in 1:speciesN) {
      sp_var_pred_values_post[[i]][j] <- (model$Sol[i, driver_name] + model$Sol[i, paste(driver_name, species_pred, species_names[j], sep = '.')])^2*var(driver_vector)
      P_Corr.tr_species_post[[i]][j] <- (model$VCV[i, year_pred] + var_pred_values_post[i])/(model$VCV[i, year_pred] + var_pred_values_post[i] + model$VCV[i, paste0(species_pred, species_names[j], '.', year_pred)] + sp_var_pred_values_post[[i]][j])
      
    }
    #P_Corr.tr - the mean P_Corr for all species in a given post row 
    P_Corr.tr_post[i] <- mean(P_Corr.tr_species_post[[i]])
  }
  return(P_Corr.tr_post)
}


#Caluculate the mean, median, mode, LCI and UCI of the metric estimates across the posterior
#variable_type = name of metric as a string
#posterior list = posterior estimates of the focal metric - the output of P_Corr/P_Div/getP_Corr.tr/getP_Div.tr
#driver_vector = a vector of unique values of the continuous predictor from the model
metric_descriptivestats <- function(variable_type, posterior_list, driver_vector) {
  
  if(typeof(posterior_list) == 'list') {
    
    mean <- apply(do.call(rbind, posterior_list), 2, mean) 
    median <- apply(do.call(rbind, posterior_list), 2, median)
    
    mode <- c()
    LCI <- c()
    UCI <- c()
    for(i in 1:length(driver_vector)){
      mode[i] <- posterior.mode(as.mcmc(do.call(rbind, posterior_list))[ ,i])
      LCI[i] <- HPDinterval(as.mcmc(do.call(rbind, posterior_list))[ ,i])[, "lower"]
      UCI[i] <- HPDinterval(as.mcmc(do.call(rbind, posterior_list))[ ,i])[, "upper"]
    }
    
    
  } else {
    
    mean <- mean(posterior_list)
    median <- median(posterior_list)
    mode <- posterior.mode(as.mcmc(posterior_list))
    LCI <- HPDinterval(as.mcmc(posterior_list))[, 'lower']
    UCI <- HPDinterval(as.mcmc(posterior_list))[, 'upper']
    
  }
  
  outputList <- list()
  outputList[[paste(variable_type, 'mean', sep = '_')]] <- mean
  outputList[[paste(variable_type, 'median', sep = '_')]] <- median
  outputList[[paste(variable_type, 'mode', sep = '_')]] <- mode
  outputList[[paste(variable_type, 'LCI', sep = '_')]] <- LCI
  outputList[[paste(variable_type, 'UCI', sep = '_')]] <- UCI
  
  return(outputList)
  
}

## Power analysis functions --------
#Converts to a mcmc object first and then takes the mode or HPDinterval (removes errors)
safe_posterior_mode <- function(x, na.rm = FALSE) {
  if(na.rm == TRUE){
    x <- na.omit(x)
  }
  posterior.mode(as.mcmc(x))
}

safe_HPDinterval <- function(x, na.rm = FALSE) {
  if(na.rm == TRUE){
    x <- na.omit(x)
  }
  HPDinterval(as.mcmc(x))
}

## Non-trending simulation ----
#param_df = the dataframe of sample sizes and variances with which to run the simulations
#rowN = the row of param_df which contains the sample sizes and variance with which to run the simulations
#reps = how many simulations to run
#nittN = how many iterations to run the model for
Non_trending_sim <- function(param_df, rowN, reps, nittN) {
  speciesN <- param_df[rowN, "speciesN_scen"]
  yearN <- param_df[rowN, "yearN_scen"]
  indN <- param_df[rowN, "indN_scen"]
  
  species_var <- param_df[rowN, "species_var_param"]
  residual_var <- param_df[rowN, "residual_var_param"]
  yearspecies_var <- param_df[rowN, "yearspecies_var_param"]
  year_var <- param_df[rowN, "year_var_param"]
  
  #Calculate the real value of the metrics based on the input variances 
  P_Corr_real <- year_var/(year_var + yearspecies_var)
  P_Div_real <- species_var/(species_var + residual_var)
  
  
  for(k in 1:reps) {
    
    #Simulate each effect
    species_effects <- rnorm(speciesN, 
                             0, 
                             sqrt(species_var))
    
    year_effects <- rnorm(yearN, 
                          0, 
                          sqrt(year_var))
    
    #Different residual for each species
    residual_effects_list <- list()
    
    for(l in 1:yearN) {
      residual_effects_list[[l]] <- list()
      for(p in 1:speciesN) {
        residual_effects_list[[l]][[p]] <- rnorm(indN,
                                                 0,
                                                 sqrt(residual_var))
      }
    }
    residual_effects <- unlist(residual_effects_list)
    
    
    #Species-specific year variance
    yearspecies_effects <- list()
    for(l in 1:speciesN){
      yearspecies_effects[[l]] <- rnorm(yearN, 
                                        0, 
                                        sqrt(yearspecies_var))
    }
    
    
    #Year effects, species effects, and year-species effects for each individual in the df
    year_effects_ind <- c()
    species_effects_ind <- c()
    yearspecies_effects_ind <- c()
    
    #Year and species as predictors
    yearpred <- c()
    speciespred <- c()
    
    
    for(m in 1:yearN) {
      newvalue <- year_effects[rep(m, 
                                   times = indN*speciesN)]
      year_effects_ind <- c(year_effects_ind, newvalue)
      
      newvalueyear <- rep(m, 
                          times = indN*speciesN)
      
      yearpred <- as.factor(c(yearpred, newvalueyear))
      
      for(l in 1:speciesN) {
        newvalue2 <- species_effects[rep(l, 
                                         times = indN)]
        species_effects_ind <- c(species_effects_ind, newvalue2)
        
        newvaluespec <- rep(l, 
                            times = indN)
        
        speciespred <- as.factor(c(speciespred, newvaluespec))
        
        ###
        new_value3 <- rep(yearspecies_effects[[l]][m], 
                          times = indN)
        yearspecies_effects_ind <- c(yearspecies_effects_ind, new_value3)
        
      }
    }
    
    # Create a data point for each individual in each year
    phenology <-
      species_effects_ind + 
      year_effects_ind + 
      yearspecies_effects_ind +
      residual_effects
    
    df <- data.frame(speciespred, yearpred, phenology)
    
    #parameter expanded priors
    pa_prior <- list(R = list(V = diag(speciesN), nu = 0.02), 
                     G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a), 
                              G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                              G1 = list(V = diag(speciesN), nu = 1, alpha.mu = rep(0, speciesN), alpha.V = diag(speciesN)*a)))
    
    
    model <- MCMCglmm(phenology ~ 1, 
                      random = ~ speciespred + 
                        yearpred + 
                        idh(speciespred):yearpred, #Species-specific year variance
                      rcov = ~ idh(speciespred):units, #Species-specific residual 
                      data = df, 
                      pr = TRUE, 
                      nitt = nittN,
                      burnin = nittN*0.25,
                      prior = pa_prior)
    
    #If the ESSs of the focal terms are below 1000, run it for twice as many iterations
    if(summary(model)$Gcovariance["speciespred", "eff.samp"] < 1000 |
       summary(model)$Gcovariance["yearpred", "eff.samp"] < 1000) {
      
      rm(model)
      
      model <- MCMCglmm(phenology ~ 1, 
                        random = ~ speciespred + 
                          yearpred + 
                          idh(speciespred):yearpred, 
                        rcov = ~ idh(speciespred):units, 
                        data = df, 
                        pr = TRUE, 
                        nitt = nittN*2,
                        burnin = (nittN*2)*0.25,
                        prior = pa_prior)
      
    }
    
    #Calculate P_corr across the posterior    
    P_Corr_post <- P_Corr(model = model, 
                          species_names = unique(speciespred),
                          year_pred = "yearpred", 
                          species_pred = "speciespred")
    
    
    #Calculate P_Spread across the posterior
    P_Div_post <- P_Div(model = model,
                        species_names = unique(speciespred),
                        species_pred = "speciespred")
    
    #Identifier key to connect outputs in different output files
    key <- paste(param_df[rowN, "scenario"], 
                 param_df[rowN, "parameter_space"], 
                 k,
                 sep = "_")
    
    #Save the outputs
    df_save <- data.frame(key = key,
                          input_speciesVar = species_var,
                          speciesVar_mean = summary.MCMCglmm(model)$Gcovariances["speciespred", "post.mean"],
                          speciesVar_median = median(model$VCV[, "speciespred"]),
                          speciesVar_mode = posterior.mode(model$VCV[, "speciespred"]),
                          speciesVar_LCI = summary(model)$Gcovariances["speciespred", "l-95% CI"],
                          speciesVar_UCI = summary(model)$Gcovariances["speciespred", "u-95% CI"],
                          speciesVar_ESS = summary(model)$Gcovariances["speciespred", "eff.samp"],
                          input_yearVar = year_var,
                          yearVar_mean = summary.MCMCglmm(model)$Gcovariances["yearpred", "post.mean"],
                          yearVar_median = median(model$VCV[, "yearpred"]),
                          yearVar_mode = posterior.mode(model$VCV[, "yearpred"]),
                          yearVar_LCI = summary(model)$Gcovariances["yearpred", "l-95% CI"],
                          yearVar_UCI = summary(model)$Gcovariances["yearpred", "u-95% CI"],
                          yearVar_ESS = summary(model)$Gcovariances["yearpred", "eff.samp"],
                          input_P_Corr = P_Corr_real,
                          P_Corr_mean = mean(P_Corr_post, na.rm = T),
                          P_Corr_median = median(P_Corr_post, na.rm = T),
                          P_Corr_mode = posterior.mode(as.mcmc(P_Corr_post)),
                          P_Corr_LCI = HPDinterval(as.mcmc(P_Corr_post))[, "lower"],
                          P_Corr_UCI = HPDinterval(as.mcmc(P_Corr_post))[, "upper"],
                          input_P_Div = P_Div_real,
                          P_Div_mean = mean(P_Div_post, na.rm = T),
                          P_Div_median = median(P_Div_post, na.rm = T),
                          P_Div_mode = posterior.mode(as.mcmc(P_Div_post)),
                          P_Div_LCI = HPDinterval(as.mcmc(P_Div_post))[, "lower"],
                          P_Div_UCI = HPDinterval(as.mcmc(P_Div_post))[, "upper"])
    
    
    
    write.table(df_save, 
                file = "Non_trendingSimulation_Output.csv", 
                append = T, 
                col.names=!file.exists("Non_trendingSimulation_Output.csv"), 
                row.names = F,
                sep = ",")
    
    
    
    ##Save output for other parameters
    
    #Create column names for the df
    desc_stats <- c("mean", "median", "mode", "LCI", "UCI", "ESS")
    
    columns_spyear <- list()
    columns_resid <- list()
    for (p in 1:length(desc_stats)) {
      columns_spyear[[p]] <- paste("speciespred", c(1:max(param_df$yearN_scen)), ".yearpred_", desc_stats[p], sep = "")
      columns_resid[[p]] <- paste("speciespred", c(1:max(param_df$yearN_scen)), ".units_", desc_stats[p], sep = "")
    }
    names(columns_spyear) <- desc_stats
    names(columns_resid) <- desc_stats
    
    columns <- c("key", 
                 "speciespred.yearpred_input",
                 "speciespred.units_input",
                 unlist(columns_spyear),
                 unlist(columns_resid))  
    
    df2_save <- data.frame(matrix(nrow = 1, ncol = length(columns))) 
    colnames(df2_save) <- columns
    
    
    df2_save[, "key"] <- key
    df2_save[, "speciespred.yearpred_input"] <- param_df[rowN, "yearspecies_var_param"]
    df2_save[, "speciespred.units_input"] <- param_df[rowN, "residual_var_param"]
    
    ##The names of the model terms to save
    speciesyear_model_names <- paste0("speciespred", 1:speciesN, ".yearpred")
    residual_model_names <- paste0("speciespred", 1:speciesN, ".units")
    
    #Species-specific year estimates 
    df2_save[, columns_spyear[["mean"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "post.mean"]
    df2_save[, columns_spyear[["median"]]] <- apply(model$VCV[, speciesyear_model_names], 2, median)
    df2_save[, columns_spyear[["mode"]]] <- apply(model$VCV[, speciesyear_model_names], 2, safe_posterior_mode)
    df2_save[, columns_spyear[["LCI"]]] <-  summary(model)$Gcovariance[speciesyear_model_names, "l-95% CI"]
    df2_save[, columns_spyear[["UCI"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "u-95% CI"]
    df2_save[, columns_spyear[["ESS"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "eff.samp"]
    
    ##Species-specific residual estimates
    df2_save[, columns_resid[["mean"]]] <- summary(model)$Rcovariance[residual_model_names, "post.mean"]
    df2_save[, columns_resid[["median"]]] <- apply(model$VCV[ ,residual_model_names], 2, median)
    df2_save[, columns_resid[["mode"]]] <- apply(model$VCV[, residual_model_names], 2, safe_posterior_mode)
    df2_save[, columns_resid[["LCI"]]] <- summary(model)$Rcovariance[residual_model_names, "l-95% CI"]
    df2_save[, columns_resid[["UCI"]]] <-  summary(model)$Rcovariance[residual_model_names, "u-95% CI"]
    df2_save[, columns_resid[["ESS"]]] <- summary(model)$Rcovariance[residual_model_names, "eff.samp"]
    
    
    write.table(df2_save, 
                file = "Non_trendingSimulation_Output_additional.csv", 
                append = T, 
                col.names=!file.exists("Non_trendingSimulation_Output_additional.csv"), 
                row.names = F,
                sep = ",")
  }
}
## Trending simulation -----
#param_df = the dataframe of sample sizes and variances with which to run the simulations
#rowN = the row of param_df which contains the sample sizes and variance with which to run the simulations
#reps = how many simulations to run
#nittN = how many iterations to run the model for
Trending_sim <- 
  function(param_df, rowN, reps, nittN) { 
    
    ##Vector of max number of years 
    max_years <- 0:(max(param_df$yearN_scen - 1))
    
    #Set the number of species, number of years and number of individuals per species
    speciesN <- param_df[rowN, "speciesN_scen"]
    yearN <- param_df[rowN, "yearN_scen"]
    indN <- param_df[rowN, "indN_scen"]
    
    #The years go from 0 to yearN-1
    years <- c(1:yearN)-1
    
    #Set the parameters 
    year_slope <- param_df[rowN, "year_slope_param"] #Slope in phenology across years 
    speciesSlope_var <- param_df[rowN, "year_slope_var_param"] #Variance in the slope across years for each species 
    speciesIntercept_var <- param_df[rowN, "year_int_var_param"] #Variance in the intercept for each species
    covar <- param_df[rowN, "slopeint_cov_param"] #Covariance between these slopes and intercepts 
    residual_var <- param_df[rowN, "residual_var_param"] #The residual variance
    
    
    residvar_slope <- param_df[rowN, "resid_slope_param"] #A slope that allows the residual variance to change across years 
    residvar_add <- c()
    for(l in 1:yearN) { #Creates this additional variance component for each year
      residvar_add[l] <- years[l]*residvar_slope
    }
    
    year_var <- param_df[rowN, "year_var_param"] #Variance in phenology across years (shared between all species)
    yearspecies_var <- param_df[rowN, "yearspecies_var_param"] #Species specific variance in phenology across years
    
    for (loop_number in 1:reps) {
      
      #Get the effects defined by each of the above parameters 
      
      #Get one year effect for each year
      year_effects <- rnorm(yearN, 
                            0, 
                            sqrt(year_var))
      
      #For each species get one year effect for each year 
      yearspecies_effects <- list()
      for(l in 1:speciesN){
        yearspecies_effects[[l]] <- rnorm(yearN, 
                                          0, 
                                          sqrt(yearspecies_var))
      }
      
      #Get species coefficient for each species - each species has an intercept and slope
      #Draw from a multivariate normal distribution - slopes and intercepts are correlated 
      species_IntsSlopes <- rmvnorm(speciesN, 
                                    c(intercept, year_slope), #mean for the distributions to draw from (intercept and slope)
                                    sigma = matrix(nrow = 2, 
                                                   c(speciesIntercept_var, 
                                                     covar, 
                                                     covar, 
                                                     speciesSlope_var)))
      
      
      #Creates a coefficient for each species for each year
      species_coeff <- list()
      for(l in 1:speciesN) {
        species_coeff[[l]] <- c(rep(NA, times = yearN))
        for(m in 1:yearN) {
          species_coeff[[l]][m] <- years[m]*species_IntsSlopes[l, 2] + species_IntsSlopes[l, 1]
        }
      }
      
      
      #Get residual effect for each individual 
      #The variance differs every year according to the slope set above
      residual_effects_list <- list()
      for(m in 1:yearN) {
        residual_effects_list[[m]] <- list()
        for(l in 1:speciesN) {
          residual_effects_list[[m]][[l]] <- rnorm(indN, 
                                                   0, 
                                                   sqrt(residual_var + residvar_add[m])) 
        } 
      }
      residual_effects <- unlist(residual_effects_list)
      
      year_effects_ind <- c()
      species_coeff_ind <- c()
      yearspecies_effects_ind <- c()
      yearpred <- c()
      speciespred <- c()
      
      
      #Repeat each effect/coefficient term the correct number of times
      #Each individual then has a species coefficient, a year effect, a species specific year effect, and a residual effect
      #Also assign what species each individual is (speciespred) and what year it is for each individual (yearpred)
      for(m in 1:yearN) {
        newvalue <- year_effects[rep(m, 
                                     times = indN*speciesN)]
        year_effects_ind <- c(year_effects_ind, newvalue)
        
        newvalueyear <- rep(years[m], 
                            times = indN*speciesN)
        
        yearpred <- c(yearpred, newvalueyear)
        
        for(l in 1:speciesN) {
          newvalue2 <- rep(species_coeff[[l]][m], 
                           times = indN)
          species_coeff_ind <- c(species_coeff_ind, newvalue2)
          
          new_value3 <- rep(yearspecies_effects[[l]][m], 
                            times = indN)
          yearspecies_effects_ind <- c(yearspecies_effects_ind, new_value3)
          
          newvaluespec <- rep(l, 
                              times = indN)
          speciespred <- as.factor(c(speciespred, newvaluespec))
          
        }
      }
      
      #Combine all the above effects to create the phenology of each individual in each year
      phenology <-
        species_coeff_ind + 
        year_effects_ind + 
        yearspecies_effects_ind +
        residual_effects 
      
      df <- data.frame(speciespred, yearpred, phenology)
      
      
      ## Modelling
      
      #Parameter expanded prior
      pa_prior <- list(R = list(V = diag(speciesN), nu = 0.02),#~ idh(species):units 
                       G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2)*a), #us(1 + yearpred):speciespred
                                G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a), #yearpred
                                G1 = list(V = diag(speciesN), nu = speciesN, alpha.mu = rep(0, speciesN), alpha.V = diag(speciesN)*a),#idh(speciespred):yearpred 
                                G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))#idh(sqrt(yearpred)):units 
      
      model <- MCMCglmm(phenology ~ yearpred, #Set year as a continuous predictor
                        random = ~ us(1 + yearpred):speciespred #Unstructured matrix. Finds the variance in species intercepts ((Intercept):(Intercept).species), the variance in species slopes (scale(year):scale(year).species) and the covariance between the two (scale(year):(Intercept).species ) 
                        + yearpred #How much year variance there is (year as a factor)
                        + idh(speciespred):yearpred #allows year means to vary for different species
                        + idh(sqrt(yearpred)):units, #fits a global slope in residual variance for all species
                        rcov = ~ idh(speciespred):units, #sets a different residual for each species
                        prior = pa_prior,
                        nitt = nittN,
                        burnin = nittN*0.25,
                        data = df, 
                        pr = TRUE)
      
      
      if(summary(model)$Gcovariance["(Intercept):(Intercept).speciespred", "eff.samp"] < 1000 |
         summary(model)$Gcovariance["yearpred:yearpred.speciespred", "eff.samp"] < 1000 |
         summary(model)$Gcovariance["yearpred", "eff.samp"] < 1000) {
        
        rm(model)
        
        model <- MCMCglmm(phenology ~ yearpred, #Set year as a continuous predictor
                          random = ~ us(1 + yearpred):speciespred #Unstructured matrix. Finds the variance in species intercepts ((Intercept):(Intercept).species), the variance in species slopes (scale(year):scale(year).species) and the covariance between the two (scale(year):(Intercept).species ) 
                          + yearpred #How much year variance there is (year as a factor)
                          + idh(speciespred):yearpred #allows year means to vary for different species
                          + idh(sqrt(yearpred)):units, #fits a global slope in residual variance for all species
                          rcov = ~ idh(speciespred):units, #sets a different residual for each species
                          prior = pa_prior,
                          nitt = nittN*2,
                          burnin = (nittN*2)*0.25,
                          data = df, 
                          pr = TRUE)
      }
      
      ##Metrics 
      
      #Calculate the the P_Div.tr and P_Corr.tr metrics 
      
      #P_Div.tr
      #Calculate for each row of the posterior 
      
      ##Calculates species variance, residual variance, and P_Div.tr for each year
      P_Div.tr_list <- getP_Div.tr(model = model,
                                   driver_vector = unique(yearpred),
                                   species_names = unique(speciespred),
                                   species_pred = "speciespred",
                                   year_pred = "yearpred")
      
      P_Div.tr_post <- P_Div.tr_list[["P_Div.tr"]]
      
      #P_Corr.tr 
      P_Corr.tr_post <- getP_Corr.tr(model = model, 
                                     driver_name = "yearpred", 
                                     driver_vector = unique(yearpred), 
                                     year_pred = "yearpred", 
                                     species_names = unique(speciespred),
                                     species_pred = "speciespred")
      
      
      #Identifier key to connect outputs in different output files
      key <- paste(param_df[rowN, "scenario"], 
                   param_df[rowN, "parameter_space"], 
                   loop_number,
                   sep = "_")
      
      #Save the necessary outputs 
      df_save <- data.frame(its = nittN,
                            key = key,
                            scenario = param_df[rowN, "scenario"],
                            parameter_set = param_df[rowN, "parameter_space"],
                            input_yearSlope = as.numeric(year_slope),
                            yearSlope_mean = summary(model)$solutions["yearpred", "post.mean"],
                            yearSlope_median = median(model$Sol[ , "yearpred"]),
                            yearSlope_mode = posterior.mode(model$Sol[ , "yearpred"]),
                            yearSlope_LCI = summary(model)$solutions["yearpred", "l-95% CI"],
                            yearSlope_UCI = summary(model)$solutions["yearpred", "u-95% CI"],
                            yearSlope_ESS = summary(model)$solutions["yearpred", "eff.samp"],
                            input_slopeVar = as.numeric(speciesSlope_var),
                            slopeVar_mean = summary(model)$Gcovariances["yearpred:yearpred.speciespred", "post.mean"],
                            slopeVar_median = median(model$VCV[ , "yearpred:yearpred.speciespred"]),
                            slopeVar_mode = posterior.mode(model$VCV[ , "yearpred:yearpred.speciespred"]),
                            slopeVar_LCI = summary(model)$Gcovariances["yearpred:yearpred.speciespred", "l-95% CI"],
                            slopeVar_UCI = summary(model)$Gcovariances["yearpred:yearpred.speciespred", "u-95% CI"],
                            slopeVar_ESS = summary(model)$Gcovariances["yearpred:yearpred.speciespred", "eff.samp"],
                            input_interceptVar = as.numeric(speciesIntercept_var),
                            interceptVar_mean = summary(model)$Gcovariances["(Intercept):(Intercept).speciespred", "post.mean"],
                            interceptVar_median = median(model$VCV[ , "(Intercept):(Intercept).speciespred"]),
                            interceptVar_mode = posterior.mode(model$VCV[ , "(Intercept):(Intercept).speciespred"]),
                            interceptVar_LCI = summary(model)$Gcovariances["(Intercept):(Intercept).speciespred", "l-95% CI"],
                            interceptVar_UCI = summary(model)$Gcovariances["(Intercept):(Intercept).speciespred", "u-95% CI"],
                            interceptVar_ESS = summary(model)$Gcovariances["(Intercept):(Intercept).speciespred", "eff.samp"],
                            input_cov = as.numeric(covar),
                            cov_mean = summary(model)$Gcovariances["yearpred:(Intercept).speciespred", "post.mean"],
                            cov_median = median(model$VCV[ , "yearpred:(Intercept).speciespred"]),
                            cov_mode = posterior.mode(model$VCV[ , "yearpred:(Intercept).speciespred"]),
                            cov_LCI = summary(model)$Gcovariances["yearpred:(Intercept).speciespred", "l-95% CI"],
                            cov_UCI = summary(model)$Gcovariances["yearpred:(Intercept).speciespred", "u-95% CI"],
                            cov_ESS = summary(model)$Gcovariances["yearpred:(Intercept).speciespred", "eff.samp"],
                            input_residSlope = as.numeric(residvar_slope),
                            residSlope_mean = summary.MCMCglmm(model)$Gcovariances["sqrt(yearpred).units", "post.mean"],
                            residSlope_median = median(model$VCV[, "sqrt(yearpred).units"]),
                            residSlope_mode = posterior.mode(model$VCV[, "sqrt(yearpred).units"]),
                            residSlope_LCI = summary.MCMCglmm(model)$Gcovariances["sqrt(yearpred).units", "l-95% CI"],
                            residSlope_UCI = summary.MCMCglmm(model)$Gcovariances["sqrt(yearpred).units", "u-95% CI"],
                            residSlope_ESS = summary.MCMCglmm(model)$Gcovariances["sqrt(yearpred).units", "eff.samp"],
                            input_yearVar = as.numeric(year_var),
                            yearVar_mean = summary.MCMCglmm(model)$Gcovariances["yearpred", "post.mean"],
                            yearVar_median = median(model$VCV[, "yearpred"]),
                            yearVar_mode = posterior.mode(model$VCV[, "yearpred"]),
                            yearVar_LCI = summary.MCMCglmm(model)$Gcovariances["yearpred", "l-95% CI"],
                            yearVar_UCI = summary.MCMCglmm(model)$Gcovariances["yearpred", "u-95% CI"],
                            yearVar_ESS = summary.MCMCglmm(model)$Gcovariances["yearpred", "eff.samp"],
                            input_P_Corr.tr = param_df[rowN, "P_Corr.tr_input"], 
                            P_Corr.tr_mean = mean(P_Corr.tr_post),
                            P_Corr.tr_median = median(P_Corr.tr_post),
                            P_Corr.tr_mode = posterior.mode(as.mcmc(P_Corr.tr_post)),
                            P_Corr.tr_LCI = HPDinterval(as.mcmc(P_Corr.tr_post))[, "lower"],
                            P_Corr.tr_UCI = HPDinterval(as.mcmc(P_Corr.tr_post))[, "upper"]
      )
      
      write.table(df_save, 
                  file = "TrendingSimulation_Output.csv", 
                  append = T, 
                  col.names =!file.exists("TrendingSimulation_Output.csv"), 
                  row.names = F,
                  sep = ",")
      
      ##Spread dataframe
      ##Years go from 0-(max(yearN-1))
      desc_stats <- c("input", "mean", "median", "mode", "LCI", "UCI")
      columns_P_Div.tr <- setNames(
        lapply(desc_stats, function(stat) paste("P_Div.tr", max_years, stat, sep = "_")),
        desc_stats
      )
      
      P_Div.tr_colnames <- c("key", unlist(columns_P_Div.tr))
      df_P_Div.tr <- setNames(
        as.data.frame(matrix(NA, ncol = length(P_Div.tr_colnames))),
        P_Div.tr_colnames
      )
      
      
      
      ##Fill in the data frame
      df_P_Div.tr$key <- key
      
      ##Convert posterior list to a matrix 
      post_matrix <- do.call(rbind, P_Div.tr_post)
      
      P_Div.tr_LCI <- c()
      P_Div.tr_UCI <- c()
      for(l in 1:length(P_Div.tr_post[[1]])){
        P_Div.tr_LCI[l] <- HPDinterval(as.mcmc(post_matrix)[ ,l])[, "lower"]
        P_Div.tr_UCI[l] <- HPDinterval(as.mcmc(post_matrix)[ ,l])[, "upper"]
      }
      
      df_P_Div.tr[ , columns_P_Div.tr[["input"]]] <- param_df[rowN, grep("P_Div.tr", colnames(param_df))]
      df_P_Div.tr[ , columns_P_Div.tr[["mean"]]] <- apply(post_matrix, 2, mean)
      df_P_Div.tr[ , columns_P_Div.tr[["median"]]] <- apply(post_matrix, 2, median)
      df_P_Div.tr[ , columns_P_Div.tr[["mode"]]] <- apply(post_matrix, 2, safe_posterior_mode)
      df_P_Div.tr[ , columns_P_Div.tr[["LCI"]]] <- P_Div.tr_LCI
      df_P_Div.tr[ , columns_P_Div.tr[["UCI"]]] <- P_Div.tr_UCI
      
      write.table(df_P_Div.tr, 
                  file = "TrendingSimulation_Output_P_Div.tr.csv", 
                  append = T, 
                  col.names =!file.exists("TrendingSimulation_Output_P_Div.tr.csv"), 
                  row.names = F,
                  sep = ",")
      
      ## Save the species year vars and resid vars in extra file 
      
      ##Save output for other parameters
      desc_stats2 <- c("mean", "median", "mode", "LCI", "UCI", "ESS")
      
      columns_spyear <- list()
      columns_resid <- list()
      for (p in 1:length(desc_stats2)) {
        columns_spyear[[p]] <- paste("speciespred", c(1:max(param_df$yearN_scen)), ".yearpred_", desc_stats2[p], sep = "")
        columns_resid[[p]] <- paste("speciespred", c(1:max(param_df$yearN_scen)), ".units_", desc_stats2[p], sep = "")
      }
      names(columns_spyear) <- desc_stats2
      names(columns_resid) <- desc_stats2
      
      columns <- c("key", 
                   "speciespred.yearpred_input",
                   "speciespred.units_input",
                   unlist(columns_spyear),
                   unlist(columns_resid))  
      
      
      #Create empty df with column names above
      df2_save <- data.frame(matrix(nrow = 1, ncol = length(columns))) 
      colnames(df2_save) <- columns
      
      #Input identifier key
      df2_save[, "key"] <- key
      
      #Save the inputs
      df2_save[ , "speciespred.yearpred_input"] <- param_df[rowN, "yearspecies_var_param"]
      df2_save[ , "speciespred.units_input"] <- param_df[rowN, "residual_var_param"]
      
      ##The names of the model terms to save
      speciesyear_model_names <- paste0("speciespred", 1:speciesN, ".yearpred")
      residual_model_names <- paste0("speciespred", 1:speciesN, ".units")
      
      #Species Year Variance 
      #Save the mean, median, CIs and ESS for species specific year variances
      df2_save[ , columns_spyear[["mean"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "post.mean"]
      df2_save[ , columns_spyear[["median"]]] <- apply(model$VCV[ , speciesyear_model_names], 2, median)
      df2_save[, columns_spyear[["mode"]]] <- apply(model$VCV[, speciesyear_model_names], 2, safe_posterior_mode)
      df2_save[ , columns_spyear[["LCI"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "l-95% CI"]
      df2_save[ , columns_spyear[["UCI"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "u-95% CI"]
      df2_save[ , columns_spyear[["ESS"]]] <- summary(model)$Gcovariance[speciesyear_model_names, "eff.samp"]
      
      #Save the mean, median, CIs and ESS for resid variances
      df2_save[ , columns_resid[["mean"]]] <- summary(model)$Rcovariance[residual_model_names, "post.mean"]
      df2_save[ , columns_resid[["median"]]] <- apply(model$VCV[ , residual_model_names], 2, median)
      df2_save[, columns_resid[["mode"]]] <- apply(model$VCV[, residual_model_names], 2, safe_posterior_mode)
      df2_save[ , columns_resid[["LCI"]]] <- summary(model)$Rcovariance[residual_model_names, "l-95% CI"]
      df2_save[ , columns_resid[["UCI"]]] <- summary(model)$Rcovariance[residual_model_names, "u-95% CI"]
      df2_save[ , columns_resid[["ESS"]]] <- summary(model)$Rcovariance[residual_model_names, "eff.samp"]
      
      write.table(df2_save, 
                  file = "TrendingSimulation_Output_additional.csv", 
                  append = T, 
                  col.names=!file.exists("TrendingSimulation_Output_additional.csv"), 
                  row.names = F,
                  sep = ",")
      
      print(paste("I am on repeat number", loop_number, sep = " "))
    }
  }