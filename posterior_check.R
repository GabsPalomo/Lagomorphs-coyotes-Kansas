(.packages())

fit
pp.check(fit, observed = )

library(coda)
library(MCMCvis)
## The output of a coda is a list of matrices where each matrix contains the output
## of a single chain for all parameters being estimated. 
## Parameter values are stored in the columns of the matrix. 
## Values for one iteration of the chain are stored in each row. 
mcmc_list <- bbsBayes::get_mcmc_list(fit) # this is a coda object
mcmc_list[[1]] # First chain 
## The mean of each parameter
MCMCvis::MCMCsummary(mcmc_list$mcmc_list)
MCMCpstr(mcmc_list$mcmc_list, params = c("pred_beta"), func = mean)


## Gelman diagram is a function of package coda that calculates the 
## potential scale reduction factor for each variable, as well as upper and lower CI
## When the upper limit is close to 1 then chains have converged. 
gelman.diag(mcmc_list$mcmc_list)
gelman.plot(mcmc_list$mcmc_list)

?jagsUI::update
params
t <- jags(data = data_list, 
          inits = list(inits1, inits2, inits3), # I need to provide initial values per chain
          parameters.to.save = params,
          model.file = 'big_pred_model.R', 
          n.chains = 3, 
          #n.adapt = na,
          n.iter = 1000,
          # n.burnin = nb, 
          # n.thin = nt, 
          parallel = T)
tu <- update(t,parameters.to.save = params, n.adapt = 1000,n.iter = 5,000)

library(tidybayes)
## All my simulations from the model fit organized in a data frame
## only for parameter pred_beta and organized by predator (1, 2, 3)
fit %>% 
  spread_draws(pred_beta[predator]) %>% 
  mutate(predator = factor(predator, levels = c(1, 2, 3), 
                          labels = c("coyote", "badger", "swift fox")))->t

t <- t %>% 
  summarise_draws()

t %>%
  # spread_draws(pred_beta[predator]) %>%
  ggplot(aes(y = predator, x = pred_beta)) +
  stat_pointinterval()+
  theme_classic()


test <- fit$sims.list
test$prey_beta[,1,1]

