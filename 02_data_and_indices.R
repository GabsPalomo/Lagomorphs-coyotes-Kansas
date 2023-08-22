## New occupancy model 
## Interaction between lagomorphs and carnivores 
## Gabby Palomo and Mason Fidino 
## July 18 2022

# Packages 
library(readr)
library(magrittr)
library(tidyverse)

source("01_data_org.R")

###########################################
# Data that needs to be supplied to model #
###########################################
#
# 1. pred_psi_design_matrix: nsite by n_parm_pred, first column should be 1's for intercept.
# 2. prey_psi_design_matrix: nsite by n_parm_prey, first column should be 1's for intercept,
#      currently assumes that every covariate you put in here has an interaction with
#      predator presence. Could be modified to change that pretty easily.
# 3. pred_y: nsite by npred by nseason array. Each element is the number of secondary sample units the
#              predator was observed at a site during a seasonal survey. If a survey
#              did not happen at that site for a given season it should be NA.
# 4. prey_y: nsite by nprey by nseason array. Each element is the number of secondary sample units the
#              prey was observed at a site during a seasonal survey. If a survey
#              did not happen at that site for a given season it should be NA.
# 5. pred_rho_design_matrix: nsite by n_parm_pred, first column should be 1's for intercept.
# 6. prey_rho_design_matrix: nsite by n_parm_prey, first column should be 1's for intercept.
# 7. J: nsite by nseason matrix. Each element is the number of secondary sample units during a primary
#         sampling period at a given site. If no sampling occured, that cell element is 0.
#
#################################################
# Indices that need to be supplied to the model #
#################################################
# 1. npred: scalar. The number of predator.
# 2. nsite: scalar. The number of unique sites across the study.
# 3. nseason: The number of primary sampling periods.
# 4. nparm_prey_psi: The number of prey parameters for psi (i.e., columns in the prey_psi_design_matrix).
# 5. nparm_prey_rho: The number of prey parameters for rho (i.e., columns in the prey_rho_design_matrix).
# 6. nprey: The number of prey species.
# 7. nparm_pred_psi: The number of predator parameters for psi (i.e., columns in the pred_psi_design_matrix).
# 8. nparm_pred_rho: The number of predator parameters for rho (i.e., columns in the pred_rho_design_matrix).
# Latent state dominant species (ds) model

# Data -------------------------------------------------------------------------
## Predator response variable --------------------------------------------------
## pred_y: nsite by npred by nseason array.
## Not the most elegant way to code this but it works. 

# make a FULL data.frame of sites & survey year, this is just
# in case not all sites are sampled each year. 

full_site <- expand.grid(
  Site = unique(sites$Site),
  SurveyYear = unique(sites$SurveyYear)
)


coydf <- data.frame(coy)
coydf$Site <- sites$Site
coydf$SurveyYear <- sites$SurveyYear
coydf <- dplyr::left_join(
  full_site,
  coydf,
  by = c("Site","SurveyYear")
)

badgerdf <- data.frame(badger)
badgerdf$Site <- sites$Site
badgerdf$SurveyYear <- sites$SurveyYear
badgerdf <- dplyr::left_join(
  full_site,
  badgerdf,
  by = c("Site","SurveyYear")
)

sfoxdf <- data.frame(sfox)
sfoxdf$Site <- sites$Site
sfoxdf$SurveyYear <- sites$SurveyYear
sfoxdf <- dplyr::left_join(
  full_site,
  sfoxdf,
  by = c("Site","SurveyYear")
)

# calculate the number of weeks of observations per year
J_vec <- length(grep("X", colnames(coydf))) -
  rowSums(is.na(coydf[,grep("^X", colnames(coydf))]))

# compare to another specie to make sure its the same

J_vec2 <- length(grep("X", colnames(sfoxdf))) -
  rowSums(is.na(sfoxdf[,grep("^X", colnames(sfoxdf))]))

if(all.equal(J_vec, J_vec2)){
  rm(J_vec2)
} else {
  stop("Different amounts of sampling between species which should be impossible")
}

# make the J matrix
J <- matrix(
  J_vec,
  ncol = length(unique(coydf$SurveyYear)),
  nrow = length(unique(coydf$Site))
)

## IM REMOVING BOBCAT SO I START HERE
# nsite by npred by nsurveyyear
pred_y <- array(
  NA,
  dim = c(nrow(J), 
          3, # CHANGED FROM 4 TO 3 BECAUSE IM REMOVING BOBCAT
          ncol(J))
)

dim(pred_y) # sites, predator species, years

# coyote
pred_y[,1,] <- rowSums(coydf[,grep("^X", colnames(coydf))], na.rm=TRUE)
# badger
pred_y[,2,] <- rowSums(badgerdf[,grep("^X", colnames(badgerdf))], na.rm=TRUE)
# sfox
pred_y[,3,] <- rowSums(sfoxdf[,grep("^X", colnames(sfoxdf))], na.rm=TRUE)

# add NA if no sampling, and do a little check to make sure
#  we are not writing over data
to_na <- which(J == 0, arr.ind = TRUE)
for(i in 1:3){ #I MODIFIED THIS TOO
  for(j in 1:nrow(to_na)){
    be_0 <-   pred_y[
      to_na[j,1],
      i,
      to_na[j,2]
    ]
    if(be_0 == 0 | is.na(be_0)){
      pred_y[
        to_na[j,1],
        i,
        to_na[j,2]
      ] <- NA
    } else{
      stop("species detection where there should be none")
    }
  }
}


## Prey response variable -----------------------------------------------------
## prey_y: nsite by nprey by nseason array.

btjrdf <- data.frame(btjr)
btjrdf$Site <- sites$Site
btjrdf$SurveyYear <- sites$SurveyYear
btjrdf <- dplyr::left_join(
  full_site,
  btjrdf,
  by = c("Site","SurveyYear")
)

ectrdf <- data.frame(ectr)
ectrdf$Site <- sites$Site
ectrdf$SurveyYear <- sites$SurveyYear
ectrdf <- dplyr::left_join(
  full_site,
  ectrdf,
  by = c("Site","SurveyYear")
)

# We have already made the J matrix, so we just need to make prey_y
prey_y <- array(
  NA,
  dim = c(nrow(J),2,ncol(J))
)
# btjr
prey_y[,1,] <- rowSums(btjrdf[,grep("^X", colnames(btjrdf))], na.rm=TRUE)
#ectrdf
prey_y[,2,] <- rowSums(ectrdf[,grep("^X", colnames(ectrdf))], na.rm=TRUE)

dim(prey_y) #sites, prey species, years

# and do the same NA stuff with prey_y. Looks like a
#  btjr
to_na <- which(J == 0, arr.ind = TRUE)
for(i in 1:2){
  for(j in 1:nrow(to_na)){
    be_0 <-   prey_y[
      to_na[j,1],
      i,
      to_na[j,2]
    ]
    if(be_0 == 0 | is.na(be_0)){
      prey_y[
        to_na[j,1],
        i,
        to_na[j,2]
      ] <- NA
    } else{
      warning("species detection where there should be none")
      print(to_na[j,])
    }
  }
}


## 1 Occupancy design matrix ------------------------------------------------------
#├ Predators -------------------------------------------------------------------
pred_psi_design_matrix <- matrix(nrow = nrow(pred_y), ncol =1)
pred_psi_design_matrix[, 1] <- 1 # intercept 
dim(pred_psi_design_matrix) #383 sites and intercept

#├ Prey -----------------------------------------------------------------------
prey_psi_design_matrix <- matrix(nrow = nrow(prey_y), ncol=5)
prey_psi_design_matrix[, 1] <- 1 # intercept
#prey_psi_design_matrix[, 2] <- as.numeric(scale(landscape$RowcropPrp_1k))
prey_psi_design_matrix[, 2] <- as.numeric(scale(landscape$PrairieSum_1k))
prey_psi_design_matrix[, 3] <- as.numeric(scale(cov_mean$VegHeight_mean))
prey_psi_design_matrix[, 4] <- as.numeric(scale(cov_mean$GrassPrp_mean))
prey_psi_design_matrix[, 5] <- as.numeric(scale(cov_mean$ForbPrp_mean))
#prey_psi_design_matrix[, 7] <- as.numeric(scale(landscape$RowcropTE))
#prey_psi_design_matrix[, 8] <- as.numeric(scale(landscape$PrairieTE_2k))
head(prey_psi_design_matrix)
dim(prey_psi_design_matrix) #383 sites and 8 covariates 

# 5. pred_rho_design_matrix: nsite by n_parm_pred, first column should be 1's for intercept.

## Occupancy design matrix ------------------------------------------------------
#├ Predators -------------------------------------------------------------------
# scent_y <- scent %>% 
#   filter(SurveyYear == 2018) %>% 
#   select(starts_with("Scent"))

pred_rho_design_matrix <- matrix(nrow = nrow(pred_y), ncol = 3)
pred_rho_design_matrix[, 1] <- 1 # intercept
pred_rho_design_matrix[, 2] <- as.numeric(scale(cov_mean$VegHeight_mean))
pred_rho_design_matrix[, 3] <- as.numeric(scale(landscape$PrairieSum_1k))
head(pred_rho_design_matrix)
dim(pred_rho_design_matrix) # 383 sites and 3 covariates for rho for each predator

# 6. prey_rho_design_matrix: nsite by n_parm_prey, first column should be 1's for intercept.
#├ Prey -----------------------------------------------------------------------
prey_rho_design_matrix <- matrix(nrow = nrow(prey_y), ncol=3)
prey_rho_design_matrix[, 1] <- 1 # intercept
prey_rho_design_matrix[, 2] <- as.numeric(scale(cov_mean$VegHeight_mean))
prey_rho_design_matrix[, 3] <- as.numeric(scale(landscape$PrairieSum_1k))
head(prey_rho_design_matrix)
dim(prey_rho_design_matrix) # 383 sites and 3 cov for rho

# 7. J: nsite by nseason matrix. Each element is the number of secondary sample units during a primary
#         sampling period at a given site. If no sampling occured, that cell element is 0.

# Mason: We have calculated this above now, not every one is a 4 because
#  there are times when no sampling occured! See above for how I calculated
# this.

# Indices ----------------------------------------------------------------------
#├ 1. npred: scalar. The number of predator ------------------------------------
npred <- dim(pred_y)[2]
#├ 2. nsite: scalar. The number of unique sites across the study ---------------
nsite <- dim(pred_y)[1] # 383 sites sampled each year, total sites across years is nrow(coy)
#├ 3. nseason: The number of primary sampling periods --------------------------
nseason <- dim(pred_y)[3] # 3 years
#├ 4. nparm_prey_psi: The number of prey parameters for psi (i.e., columns in the prey_psi_design_matrix).
nparm_prey_psi <- ncol(prey_psi_design_matrix)
# 5. nparm_prey_rho: The number of prey parameters for rho (i.e., columns in the prey_rho_design_matrix).
nparm_prey_rho <- ncol(prey_rho_design_matrix)
# 6. nprey: The number of prey species.
nprey <- dim(prey_y)[2]
# 7. nparm_pred_psi: The number of predator parameters for psi (i.e., columns in the pred_psi_design_matrix).
nparm_pred_psi <- dim(pred_psi_design_matrix)[2]
# 8. nparm_pred_rho: The number of predator parameters for rho (i.e., columns in the pred_rho_design_matrix).
nparm_pred_rho <- dim(pred_rho_design_matrix)[2]
# Latent state dominant species (ds) model
