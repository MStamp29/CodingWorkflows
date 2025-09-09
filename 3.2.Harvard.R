##########################################################
##########################################################
# This runs both the non-trending and trending models 
# on the Harvard Forest dataset (Harvard_springphen.csv) 
# and calculates the correlation and divergence metrics.
##########################################################
##########################################################

#Sets the wd to source file location on RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MCMCglmm)
library(mvtnorm)

rm(list = ls())

#functions
source('1.PhenMetric_functions.R')

#1. Data prep ------
harvard <- read.csv("Harvard_springphen.csv")

#Subset to those species that appear basically every year (17 species)
#Those in 2021 are those that appear almost every year
harvard_sub <- harvard[harvard$species %in% 
                         unique(harvard[which(harvard$year == 2021), ]$species), ]

speciesName <- levels(as.factor(harvard_sub$species))
speciesName_English <- c('Striped maple', 'Red maple', 'Sugar maple',
                         'Shadbush', 'Yellow birch', 'Black birch', 
                         'Paper birch', 'Dogwood',
                         'Hawthorne', 'Beech', 'White ash', 'Witch hazel',
                         'Trembling aspen', 'Black cherry', 'White oak', 
                         'Red oak', 'Black oak')
speciesN <- nlevels(as.factor(harvard_sub$species))



#2. Run the stationary model -----
#parameter expanded priors
a <- 1000
pa_prior <- list(R = list(V = diag(speciesN), nu = 0.02), 
                 G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a), 
                          G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                          G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                          G1 = list(V = diag(speciesN), nu = 1, alpha.mu = rep(0, speciesN), alpha.V = diag(speciesN)*a)))

nittN <- 300000 
burninN <- nittN*0.25

model <- MCMCglmm(bb.doy ~ 1, 
                  random = ~ species 
                  + year 
                  + tree.id  ### acknowledge that flawed
                  + idh(species):year, 
                  rcov = ~ idh(species):units, 
                  data = harvard_sub, 
                  pr = TRUE, 
                  nitt = nittN,
                  burnin = burninN,  
                  prior = pa_prior)

saveRDS(model, file = 'harvard_NonTrendingModel.rds')
#model <- readRDS('harvard_NonTrendingModel.rds')

#2.1. P_Corr ------
P_Corr_post <- P_Corr(model = model, 
                      species_names = speciesName,
                      year_pred = "year", 
                      species_pred = "species")

#Get the mean, median, mode, and CIs across the posterior
P_Corr_descriptivestats <- metric_descriptivestats(variable_type = 'P_Corr', 
                                                   posterior_list = P_Corr_post, 
                                                   driver_vector = NULL)


#2.2. P_Div -------
#Across the posterior distribution 
P_Div_post <- P_Div(model = model,
                    species_names = speciesName,
                    species_pred = "species")

#Get the mean, median, mode, and CIs across the posterior
P_Div_descriptivestats <- metric_descriptivestats(variable_type = 'P_Div',
                                                  posterior_list = P_Div_post,
                                                  driver_vector = NULL)


#4. Run the dynamic model --------
#Flip order of years if want the residual global slope to be decreasing
years <- unique(harvard_sub$year)
yearN <- length(years)
inv_year <- c(0:(yearN-1))[order(0:(yearN-1), decreasing = T)]
asc_year <- c(0:(yearN-1))
harvard_sub$inv_year <- NA

for(i in 1:nrow(harvard_sub)) {
  for(j in 1:length(years)) {
    if(harvard_sub$year[i] == years[j]) {
      harvard_sub$inv_year[i] <- inv_year[j]
    }
  }
}

#Parameter expanded priors
a2 <- 1000
pa_prior2 <- list(R = list(V = diag(speciesN), nu = 0.02),#~ idh(species):units 
                  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2)*a2), #~ us(1 + scaledinv_year):species
                           G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a2), #year,
                           G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a2), #tree.id
                           G1 = list(V = diag(speciesN), nu = speciesN, alpha.mu = rep(0, speciesN), alpha.V = diag(speciesN)*a2),#idh(species):year 
                           G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a2)))#idh(sqrt(inv_year)):units 

nittN_tr <- 1000000
burninN_tr <- nittN_tr*0.25

model2 <- MCMCglmm(bb.doy ~ inv_year, #Set year as a continuous predictor - scale so that it is meaningful (otherwise it sets the intercept as year 0)
                   random = ~ us(1 + inv_year):species #us - unstructured matrix. Finds the variance in species intercepts ((Intercept):(Intercept).species), the variance in species slopes (scale(year):scale(year).species) and the covariance between the two (scale(year):(Intercept).species ) 
                   + inv_year #How much year variance there is (year as a factor)
                   + tree.id
                   + idh(species):inv_year #allows year means to vary for different species
                   + idh(sqrt(inv_year)):units, #fits a global slope in residual variance for all species
                   rcov = ~ idh(species):units, #sets a different residual for each species
                   prior = pa_prior2, 
                   nitt = nittN_tr,
                   burnin = burninN_tr,
                   data = harvard_sub, 
                   pr = TRUE)

saveRDS(model2, file = 'harvard_TrendingModel.rds')
#model2 <- readRDS('harvard_TrendingModel.rds')


## 4.3. Calculate metrics for dynamic models  ------

#4.3.1. P_Div.tr
specvar_P_Div.tr_post <- getP_Div.tr(model = model2,
                                     driver_vector = inv_year,
                                     species_names = speciesName,
                                     species_pred = "species",
                                     year_pred = "inv_year")

#saveRDS(specvar_P_Div.tr_post, 'Speciesvariance_P_DivTr_Resid_posterior.rda')

specvar_post <- specvar_P_Div.tr_post[['Species var']]
P_Div.tr_post <- specvar_P_Div.tr_post[['P_Div.tr']]
residvar_post <- specvar_P_Div.tr_post[['Resid var']]


P_Div.tr_descriptivestats <- metric_descriptivestats(variable_type = 'P_Div.tr',
                                                     posterior_list = P_Div.tr_post,
                                                     driver_vector = inv_year)

specvar_descriptivestats <- metric_descriptivestats(variable_type = 'specvar',
                                                    posterior_list = specvar_post, 
                                                    driver_vector = inv_year)

residvar_descriptivestats <- metric_descriptivestats(variable_type = 'residvar',
                                                     posterior_list = residvar_post, 
                                                     driver_vector = inv_year)


df_P_Div.tr <- cbind(years,
                     asc_year, #This is the 'non-inverted' year
                     as.data.frame(do.call(cbind, P_Div.tr_descriptivestats)),
                     as.data.frame(do.call(cbind, specvar_descriptivestats)),
                     as.data.frame(do.call(cbind, residvar_descriptivestats)))


#write.csv(df_P_Div.tr, 'P_DivTr_year.csv', row.names = F)

##4.3.2. ICC cor ------
P_Corr.tr_post <- getP_Corr.tr(model = model2, 
                               driver_name = "inv_year", 
                               driver_vector = inv_year, 
                               year_pred = "inv_year", 
                               species_names = speciesName,
                               species_pred = "species")


P_Corr.tr_post_descriptivestats <- metric_descriptivestats(variable_type = 'P_Corr.tr',
                                                           posterior_list = P_Corr.tr_post,
                                                           driver_vector = inv_year)

