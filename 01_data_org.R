## New occupancy model 
## Interaction between lagomorphs and carnivores 
## Gabby Palomo and Mason Fidino 
## July 18 2022

# Packages 
library(readr)
library(magrittr)
library(tidyverse)

# Functions 
source("functions/functions.R")

# Detection histories  -------------------------------------------------------------------------
# detection histories obtained from the MasterCarnivoreStudyOccAndCovariates 
#├ Lagomorphs ------------------------------------------------------------------
btjr <- read_csv("data/BTJR_dh.csv") %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE)
ectr <- read_csv("data/ECTR_dh.csv") %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE)
#├ Carnivores ------------------------------------------------------------------
coy <- read_csv("data/coyote_dh_.csv") %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE)
badger <- read_csv("data/badger_dh.csv") %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE)
sfox <- read_csv("data/sfox_dh.csv") %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE)

# btjr and ectr don't have all the sites because they were NAs so we need to 
# input those sites and set the values as NA 
# sites tell us that we have 1149 rows so all our detection histories have to have
# 1149 rows. 
sites <- read_csv("data/sites.csv")
empty <-as.data.frame(lapply(coy, function(x) rep.int(NA, length(x))))
empty$Site <- coy$Site
empty$SurveyYear <- coy$SurveyYear
# Remove columns that start with'Coyote' and leave only Site and survey year
empty <- empty[,-grep("Coyote", colnames(empty))]

# Join the empty df with the other species detection histories 
btjr <- left_join(empty, btjr, by = c("SurveyYear", "Site"))
ectr <- left_join(empty, ectr, by = c("SurveyYear", "Site"))

#├ Shrink dh to occasions ---------------------------------------------------------
btjr <- shrink(btjr[, grep("Occ", colnames(btjr))], 7) # each occasion is composed of 7 days
ectr <- shrink(ectr[, grep("Occ", colnames(ectr))], 7)
coy <- shrink(coy[,grep("Occ", colnames(coy))], 7)
badger <- shrink(badger[,grep("Occ", colnames(badger))], 7)
sfox <- shrink(sfox[,grep("Occ", colnames(sfox))], 7)

# Detection covariates -----------------------------------------------------
doy <- read_csv("data/doy.csv")
scent <- read_csv("data/scent.csv")
days_active <- read_csv("data/days_active.csv")

doy <- doy %>% 
  select(Site, SurveyYear, DOY_7, DOY_14, DOY_21, DOY_28) %>% 
  mutate(across(3:6, round, 0)) # doy rounded so no decimals but just number of days

scent <- scent %>% 
  select(Site, SurveyYear, Scent_7, Scent_14, Scent_21, Scent_28)

## year matrix
site_year <- doy %>% 
  group_by(SurveyYear) %>% 
  summarise(total = n())

# make it a factor 
yeardf <- data.frame(
  year = factor(
    rep(
      c("2018", "2019", "2020"),
      times = c(383,383,383)
    )
  )
)
# make dummy response variable
yeardf$y <- 1

# Now we create the dummy variables for years. 2018 will be the reference year
# so if I have Year2019 = 0 and Year2020 = 0 then it must be 2018 (0 0); 
# If I have Year2019 =1 and Year2020 = 0 then it must be 2019 (1 0);
# If I have Year2019 = 0 and Year2020 = 1 then it must be 2020 (0 1)
year <- model.matrix(y~year, data = yeardf) 
# We need to remove the intercept column (for 2018)
year <- year[,grep("year", colnames(year))]

# Occupancy covariates -------------------------------------------------------------
#├ Fine scale covariates ------------------------------------------------------
cov <- read_csv("data/fine_scale_habitat.csv") %>% 
  relocate(SurveyYear, .before=Site) %>% 
  distinct(SurveyYear, Site, .keep_all = TRUE) %>% 
  mutate(across(GrassPct:OpenGroundPct, ~.x/100)) %>% 
  rename_all(~(str_replace(., "Pct", "Prp"))) 

#days_active <- left_join(empty, days_active, by = c("SurveyYear", "Site"))

# Cov has less rows because not all sites were surveyed all years 
cov <- left_join(empty, cov, by = c("SurveyYear", "Site")) %>% 
  select(1:8)

#cov$DaysActive <- days_active$DaysActive

# mean imputing by site if there are any NA's
cov <- split(cov, factor(cov$Site))
vcols <- colnames(cov[[1]])[-c(1:2)]

cat("Imputting missing covariate data...\n")
pb <- txtProgressBar(max = length(cov))
for(i in 1:length(cov)){
  setTxtProgressBar(pb, i)
  tmp <- cov[[i]]
  if(any(is.na(tmp))){
    #stop()
    for(j in 1:length(vcols)){
      tmp[is.na(tmp[,vcols[j]]),vcols[j]] <- mean(tmp[,vcols[j]], na.rm = TRUE)
    }
  }
  cov[[i]] <- tmp
}
cov <- dplyr::bind_rows(cov) # we unlist the elements and create a df
dim(cov)
 
#Now let's estimate the mean of each covariate across years. 
cov_mean <- cov %>% 
  group_by(Site) %>% 
  select(VegHeight, GrassPrp, ForbPrp) %>% 
  summarise(across(everything(), .f = list(mean = ~mean(.[!is.na(.) | . !=0]))))

## Site 294 (row 295) and site 339 (row 340) have NA across all years for all cov, so just to test the model
## I will place a 0 here but will have to figure how to deal with the data at some point. 
cov_mean[295, 2:4] <- 0
cov_mean[340, 2:4] <- 0

#├ Landscape covariates ------------------------------------------------------
landscape <- read_csv("data/landcover_scales.csv")

landscape <- landscape %>%  
  ## sum across all prairies 
  mutate(PrairieSum_750 = SGPPrp_750+MGPPrp_750+TGPPrp_750+SandsagePrp_750) %>% 
  mutate(PrairieSum_1k = SGPPrp_1k+MGPPrp_1k+TGPPrp_1k+SandsagePrp_1k) %>% 
  mutate(PrairieSum_1.5k = SGPPrp_1.5k+MGPPrp_1.5k+TGPPrp_1.5k+SandsagePrp_1.5k) %>% 
  mutate(PrairieSum_2k = SGPPrp_2k+MGPPrp_2k+TGPPrp_2k+SandsagePrp_2k) %>% 
  mutate(PrairieTE_2k = SGPTE+MGPTE+TGPTE+SandsageTE) %>% 
  filter(SurveyYear == 2018) 

# END------------------------------------------------------------------------
