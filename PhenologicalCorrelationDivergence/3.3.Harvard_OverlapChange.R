############################################################
############################################################
# Calculates a) the overlap in phenological distributions 
# (using estimates of the mean and standard deviation for 
# each species from the trending model in Harvard.R) in 
# the first year of the Harvard Forest dataset, and 
# b) the change in overlap between the first and last year. 
############################################################
############################################################

library(MCMCglmm)

#Sets the wd to source file location on RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())

## 1. Data import and prep ------
#Load the trending Harvard model output
model <- readRDS("harvard_TrendingModel.rds") #Output from 3.2.Harvard.R
fixed_effects_posterior <- model$Sol
random_effects_posterior <- model$VCV

#Data on which the trending model was run
harvard <- read.csv("Harvard_springphen.csv")
harvard_sub <- harvard[harvard$species %in% 
                         unique(harvard[which(harvard$year == 2021), ]$species), ]


species_names <- unique(harvard_sub$species)
names(species_names) <-  c('Striped maple', 'Red maple', 'Sugar maple','Shadbush', 'Yellow birch', 'Black birch', 'Paper birch', 'Dogwood', 'Hawthorne', 'Beech', 'White ash', 'Witch hazel','Aspen', 'Black cherry', 'White oak', 'Red oak', 'Black oak')

## 2. Get species means and residual variance in first and last years for whole posterior ------

#species mean = (Intercept + SpDev in Intercept) + (Slope + SpDev in Slope)*Year
#NB: Years are inverted so Year 0 is the last year and Year 31 is the first year

# residual variance = baseline residual variance + (Residual slope*sqrt(Year))
# I.e. (species).units + sqrt(inv.year)*sqrt(Year)

Intercept <- fixed_effects_posterior[,"(Intercept)"]
Slope <- fixed_effects_posterior[,"inv_year"]
Resid_slope <- random_effects_posterior[, "sqrt(inv_year).units"]

Intercept_deviations_colnames <- paste0("\\(Intercept\\)\\.species\\.", species_names) #\\ escapes the parentheses and periods (when used with grep below)
Slope_deviations_colnames <- paste0("inv_year\\.species\\.", species_names)
Species_resids_colnames <- paste0("species", species_names, "\\.units")

Intercept_deviations <- fixed_effects_posterior[ ,grep(paste(Intercept_deviations_colnames, collapse = "|"), colnames(fixed_effects_posterior))]
Slope_deviations <- fixed_effects_posterior[ ,grep(paste(Slope_deviations_colnames, collapse = "|"), colnames(fixed_effects_posterior))]

#Species_intercepts are species estimates in year 0
#Species_resids are residual variances in year 0
Species_intercepts <- Intercept_deviations + Intercept
Species_slopes <- Slope_deviations + Slope
Species_resids <- random_effects_posterior[ ,grep(paste(Species_resids_colnames, collapse = "|"), colnames(random_effects_posterior))]

colnames(Species_intercepts) <- species_names
colnames(Species_slopes) <- species_names
colnames(Species_resids) <- species_names

#Year 0
Species_means_0 <- Species_intercepts
Species_resids_0 <- Species_resids
Species_stdev_0 <- sqrt(Species_resids)

#Get the estimate and residual variance for each species in year 31
Species_means_31 <- Species_intercepts + (Species_slopes*31)
Species_resids_31 <- Species_resids + Resid_slope*sqrt(31)
Species_stdev_31 <- sqrt(Species_resids_31)

#Get all pairs of species names 
spname_combinations <- as.data.frame(t(combn(species_names, 2)))
combo_names <- apply(spname_combinations, 1, paste, collapse = "_")

#Create lists to store the initial overlap and the change in overlap (Differences)
Initial_overlap <- vector("list", length(combo_names))
names(Initial_overlap) <- combo_names

Differences <- vector("list", length(combo_names))
names(Differences) <- combo_names


#Get the corresponding column index for each species name in the means/residual model output objects (in case they are not in the same order)
species_col_indices <- lapply(species_names, function(species) {
  list(index = grep(species, colnames(Species_means_0)))
})

names(species_col_indices) <- species_names


##Functions which calculate the overlap between two distributions based on their means (mu) and sds
min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

calculate_overlap <- function(mu1, sd1, mu2, sd2, n) {
  lower_limit <- min(mu1 - 3 * sd1, mu2 - 3 * sd2)
  upper_limit <- max(mu1 + 3 * sd1, mu2 + 3 * sd2)
  
  integrate(
    min.f1f2, lower_limit, upper_limit, 
    mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2
  )$value
}

for (i in 1:nrow(fixed_effects_posterior)) {
  for (j in seq_along(combo_names)) {
    Species1 <- spname_combinations[j, 1]
    Species2 <- spname_combinations[j, 2]
    
    # Indices for Species1 and Species2
    index_Sp1 <- species_col_indices[[Species1]]$index
    index_Sp2 <- species_col_indices[[Species2]]$index
    
    
    # Means and sd for each species in the pair in year 0
    mu1 <- Species_means_0[i, index_Sp1]
    sd1 <- Species_stdev_0[i, index_Sp1]
    mu2 <- Species_means_0[i, index_Sp2]
    sd2 <- Species_stdev_0[i, index_Sp2]
    #Calculate the overlap between the two distributions in year 0
    overlap_0 <- calculate_overlap(mu1, sd1, mu2, sd2, n)
    
    # Means and sd for each species in the pair in year 31
    mu3 <- Species_means_31[i, index_Sp1]
    sd3 <- Species_stdev_31[i, index_Sp1]
    mu4 <- Species_means_31[i, index_Sp2]
    sd4 <- Species_stdev_31[i, index_Sp2]
    #Calculate the overlap between the two distributions in year 31
    overlap_31 <- calculate_overlap(mu3, sd3, mu4, sd4, n)
    
    # Calculate the difference
    Differences[[combo_names[j]]][i] <- overlap_0 - overlap_31
    
    Initial_overlap[[combo_names[j]]][i] <- overlap_31
  }
}

#saveRDS(Differences, "Overlap_differences.RDS")
#saveRDS(Initial_overlap, "Initial_overlap.RDS")
