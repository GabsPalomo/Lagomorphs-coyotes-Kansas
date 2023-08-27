## Packages 
library(jagsUI)
source("02_data_and_indices.R")

## Put data in a list ---------------------------------------------------------
data_list <- list(
  pred_psi_design_matrix = pred_psi_design_matrix,
  prey_psi_design_matrix = prey_psi_design_matrix,
  pred_y = pred_y, 
  prey_y = prey_y, 
  pred_rho_design_matrix = pred_rho_design_matrix, 
  prey_rho_design_matrix = prey_rho_design_matrix, 
  J = J, 
  npred = npred,
  nsite = nsite, 
  nseason = nseason, 
  nparm_prey_psi = nparm_prey_psi,
  nparm_prey_rho = nparm_prey_rho,
  nprey = nprey, 
  nparm_pred_psi = nparm_pred_psi, 
  nparm_pred_rho = nparm_pred_rho
)

# get parameters to track, that is the helpful part with providing initial
#  values to all of your parameters!
params <- names(my_inits(1))
# don't track the z array (unless you really want to)
params <- params[-grep("_z$|RNG", params)]

## These parms correspond with what is in the model 

## initial values -----------------------------------------------------------

inits1 <- my_inits(3)
inits2 <- my_inits(3)
inits3 <- my_inits(3)
inits4 <- my_inits(3)

# MCMC settings ----------------------------------------------------------------
nc <- 3 # number of chains
na <- 1000 # number of adaptations to run in the adaptive phase 
ni <- 530000  # number of iterations
nb <- 30000 # number of burning 
nt <- 10 # number of thinning 

iterations <- ((ni-nb)/nt)*nc
iterations #60,000 
(iterations/nc) #20,000 per chain (3)

##jagsui will produce (niter-nburnin)/nthin as the number of steps per chain as an end result

# fitting the model
fit <- jags(data = data_list, 
            inits = list(inits1, inits2, inits3), # I need to provide initial values per chain
            parameters.to.save = params,
            model.file = 'big_pred_model.R', 
            n.chains = nc, 
            n.adapt = na,
            n.iter = ni,
            n.burnin = nb, 
            n.thin = nt, 
            parallel = T)

# Results MCMC -------------------------------------------------------------

# saveRDS(fit, "2022_09_15_fit.rds")
fit <- readRDS('./rds/2022_09_15_fit.rds')

fit
summary(fit)
# Summary for model 'big_pred_model.R' 
# Saved parameters: pred_beta pred_theta pred_alpha inxs_b0 inxs_beta prey_theta prey_beta prey_alpha deviance 
# MCMC ran in parallel for 1079.712 minutes at time 2022-09-09 16:11:03.
# 
# For each of 3 chains:
#   Adaptation:            1000 iterations (sufficient)
# Burn-in:               30000 iterations
# Thin rate:             10 iterations
# Total chain length:    531000 iterations
# Posterior sample size: 50000 draws
# 
# Successful convergence based on Rhat values (all < 1.1). 
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 4354.2 and DIC = 10177.93 


names(fit)
# plot(fit)
# traceplot(fit)
## Let's use MCMCvis to visualize the traceplots 
library(MCMCvis)
MCMCplot(fit,
         params = params,
         horiz = FALSE, 
         rank = T,
         ref_ovl = FALSE)

MCMCsummary(fit)
# Put summary in a table and export
results <- MCMCsummary(fit, round = 2)
library(kableExtra)
results %>%
  kbl() %>%
  kable_styling(full_width = FALSE) 
# %>%
#   save_kable(file = "2022_09_15_summary.pdf", self_contained = TRUE)

# Export pdf with trace plots and density plots for each parameter
# MCMCtrace(fit,
#           iter = 10000,
#           filename = "2022_09_15_traceplot.pdf",
#           wd = getwd())

# Overall Prey and Predator occupancy results across all sites ----------------------------------------------------------
# You don't want to make predictions from the mean of each posterior.
# Rather, you want to do the predictions across the posterior and then calculate
# the mean (or median). 
mcmc_list <- bbsBayes::get_mcmc_list(fit) # to see objects in mcmc 
# we see in mcmc_list that there is a list called sims.list which is what we want
tmp <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(tmp)
# "pred_beta"  "pred_theta" "pred_alpha" "inxs_b0"    "inxs_beta"  "prey_theta"
# "prey_beta"  "prey_alpha" "deviance" 

# Average occupancy across study for prey ---------------------------------------
# Summary and traceplots
print(fit_sum <- MCMCsummary(fit, round = 2))
# write.csv(fit_sum, "2022_09_07_fit.csv")
mcmc_list <- bbsBayes::get_mcmc_list(fit) # to see objects in mcmc 
# we see in mcmc_list that there is a list called sims.list which is what we want
tmp <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(tmp)
# "pred_beta"  "pred_theta" "pred_alpha" "inxs_b0"    "inxs_beta"  "prey_theta"
# "prey_beta"  "prey_alpha" "deviance" 

prey_order <- c("btjr", "ectr")

#Prey occupancy estimate overall sites
psi1prey <- plogis(tmp$prey_beta[,,1])
psi2prey <- plogis(tmp$prey_beta[,,1] + tmp$prey_theta)
eocprey <- psi1prey / (psi1prey + (1-psi2prey))
head(eocprey)

# Use apply to estimate eoc for all species (columns in the tmp$pred array)
prey_occ1 <- t(
  apply(
    eocprey, # data vector 
    2, # it will be applied to columns, 1 means it will be applied to rows.
    quantile, # function to apply to eoc data
    probs = c(0.025,0.5,0.975) #argument of quantile function 
  )
)
prey_occ1

# Predator influence on prey occupancy overall sites
# Above: This would be prey occupancy given the absence of all the other
#  predator species. If we want to get the average we should also 
#  add in the expected occupancy of the predator species
# predator order 
pred_order <- c('coy', 'badger', 'sfox')
# Predator estimate across all sites 
pred_occ <- plogis(tmp$pred_beta)
pred_occ2 <- plogis(tmp$pred_beta[,,1] + tmp$pred_theta)
eocpred <- pred_occ[,,1] / (pred_occ[,,1] + (1 - pred_occ2))
eocpred <- apply(eocpred, 2, median)

# get predator 'influence' on prey given their expected occupancy
# rows are iterations, columns are prey, arrays are predators [iterations, prey, predators]
pred_inf_prey1 <- sweep(tmp$inxs_b0[,1,],2, eocpred,"*") #jackrabbit
pred_inf_prey2 <- sweep(tmp$inxs_b0[,2,],2, eocpred,"*") #cottontail
pred_inf_prey1 <- rowSums(pred_inf_prey1)
pred_inf_prey2 <- rowSums(pred_inf_prey2)
pred_inf_prey <- cbind(pred_inf_prey1, pred_inf_prey2)

# prey given eoc(predator) and the influence they have on prey
psi1prey_gp <- plogis(
  tmp$prey_beta[,,1] + pred_inf_prey
)
psi2prey_gp <- plogis(tmp$prey_beta[,,1] + tmp$prey_theta + pred_inf_prey)
eocprey_gp <- psi1prey_gp / (psi1prey_gp + (1-psi2prey_gp))

colMeans(eocprey_gp)

# Use apply to estimate eoc for all species (columns in the tmp$pred array)
prey_occ2 <- t(
  apply(
    eocprey_gp, # data vector 
    2, # it will be applied to columns, 1 means it will be applied to rows.
    quantile, # function to apply to eoc data
    probs = c(0.025,0.5,0.975) #argument of quantile function 
  )
)

prey_occ2 <- as.data.frame(prey_occ2, )
rownames(prey_occ2) <- NULL

prey_occ <- rbind(prey_occ1, prey_occ2)
prey_occ

prey_occ <- prey_occ %>% 
  mutate(species = c('btjr', 'ectr', 'btjr', 'ectr')) %>% 
  mutate(influence = c('no predator', 'no predator', 'predator', 'predator')) %>% 
  relocate(species, influence)

prey_occ

prey_occ[1, 1] <- "jackrabbit"
prey_occ[2, 1] <- "cottontail"
prey_occ[3, 1] <- "jackrabbit"
prey_occ[4, 1] <- "cottontail"

# Reorder so they are all organized by species 
prey_occ <- prey_occ[c(1, 3, 2, 4),]
row.names(prey_occ) <- NULL

# Put in a nice table to export 
library(kableExtra)
prey_occ %>% 
  select(species, influence, '2.5%', '50%', '97.5%') %>% 
  mutate(visualization = "") %>% 
  relocate(species, influence, "50%", visualization) %>% 
  rename("occupancy" = "50%") %>% 
  kbl(digits = 2) %>% 
  kable_styling(font_size = 14) %>% 
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  add_header_above(c("" , "", "", "", "Confidence intervals" = 2), bold = TRUE) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(4, 
              image = spec_pointrange(
                x = prey_occ[,4],
                xmin = prey_occ[,3], 
                xmax = prey_occ[,5],
                vline = 0.5
              ))-> prey_occ_g
prey_occ_g

#├ Average occupancy across study for predators.--------------------------------
# Summary and traceplots
print(fit_sum <- MCMCsummary(fit, round = 2))
# write.csv(fit_sum, "2022_09_07_fit.csv")
mcmc_list <- bbsBayes::get_mcmc_list(fit) # to see objects in mcmc 
# we see in mcmc_list that there is a list called sims.list which is what we want
tmp <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(tmp)
# "pred_beta"  "pred_theta" "pred_alpha" "inxs_b0"    "inxs_beta"  "prey_theta"
# "prey_beta"  "prey_alpha" "deviance" 

predator_order <- c("coyote", "badger", "swift fox")

## Predator occupancy estimate overall sites
psi1 <- plogis(tmp$pred_beta[,,1]) #vector
psi2 <- plogis(tmp$pred_beta[,,1] + tmp$pred_theta) #vector
eoc <- psi1 / (psi1 + (1 - psi2)) #vector
# Use apply to estimate eoc for all species (columns in the tmp$pred array)
pred_occ <- t(
  apply(
    eoc, # data vector 
    2, # columns
    quantile, # function to apply to eoc data
    probs = c(0.025,0.5,0.975) #argument of quantile function 
  )
)

prey_occ %>% 
  as.data.frame() %>% 
  mutate(visualization = "") %>% 
  relocate(species, influence, "50%", visualization) %>% 
  rename("Estimate" = "50%") %>% 
  rename("Species" = "species") %>% 
  rename("Influence" = "influence") %>% 
  rename("Visualization" = "visualization")-> t.prey.occ

t.prey.occ

pred_occ %>% 
  as.data.frame() %>% 
  mutate(species = predator_order,
         visualization = "", 
         influence = "-") %>% 
  relocate(species, influence, "50%", visualization) %>% 
  rename("Estimate" = "50%") %>% 
  rename("Species" = "species") %>% 
  rename("Influence" = "influence") %>% 
  rename("Visualization" = "visualization")-> t.pred.occ

t.pred.occ


t.all.occ <- rbind(t.prey.occ, t.pred.occ)
t.all.occ

# occ_g %>% 
#   save_kable("figures_06/occ_table_nobob.png", 
#              self_contained = TRUE,
#              density = 500)


# Overall prey detection estimate across all sites. -------------------------------------------
# same thing for detection, dont need to do it piecewise
prey_det <- t(
  apply(
    plogis(tmp$prey_alpha[,,1]), #This is what I want to apply to all the columns in my array
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)

prey_order <- c("btjr", "ectr")
prey_det
prey_det %>% 
  as.data.frame() %>% 
  mutate(species = prey_order,
         visualization = "",
         influence = "-") %>% 
  relocate(species, influence, "50%", visualization) %>% 
  rename("Estimate" = "50%") %>% 
  rename("Species" = "species") %>% 
  rename("Influence" = "influence") %>% 
  rename("Visualization" = "visualization") -> prey.det
prey.det[1, 1] <- "jackrabbit"
prey.det[2, 1] <- "cottontail"


occ.det <- rbind(t.all.occ, prey.det)
occ.det


# Overall predator detection estimate across all sites. -------------------------------------------
# same thing for detection, dont need to do it piecewise
pred_det <- t(
  apply(
    plogis(tmp$pred_alpha[,,1]), #This is what I want to apply to all the columns in my array
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)
pred_det
pred_det %>% 
  as.data.frame() %>% 
  mutate(species = predator_order,
         visualization = "",
         influence = "-") %>% 
  relocate(species, influence, "50%", visualization) %>% 
  rename("Estimate" = "50%") %>% 
  rename("Species" = "species") %>% 
  rename("Influence" = "influence") %>% 
  rename("Visualization" = "visualization") -> pred_det

pred_det

all.t <- rbind(occ.det, pred_det)

all.t

all.t %>% 
  kbl(digits = 2) %>% 
  kable_styling(font_size = 14) %>% 
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  add_header_above(c("" , "", "", "", "Credible intervals" = 2), bold = TRUE) %>% 
  row_spec(0, bold = TRUE) %>% 
  pack_rows("Occupancy", 1,7) %>% 
  pack_rows("Prey", 1, 4) %>%
  pack_rows("Predator", 5, 7) %>% 
  pack_rows("Detection", 8, 12) %>%  
  pack_rows("Prey", 8, 9) %>% 
  pack_rows("Predator", 10, 12) %>% 
  row_spec(1:7, align = 'c') %>% 
  column_spec(4, 
              image = spec_pointrange(
                x = all.t[,3],
                xmin = all.t[,5], 
                xmax = all.t[,6],
                vline = 0.5)) -> det_g
det_g

# Export the final table as a png file
# det_g %>%
#   save_kable("tables/occ_det_table_final.png",
#              self_contained = TRUE,
#              density = 700)

# END -------------------------------------------------------------------------
