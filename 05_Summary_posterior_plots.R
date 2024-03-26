## Summary table 
## Gabby Palomo and Mason Fidino 

# packages ------------------------------------------------------------
library(jagsUI)
library(MCMCvis)
library(ggthemes)
library(rphylopic)
library(bayesplot)
library(cowplot)
library(kableExtra)

# Source functions and data ----------------------------------------------
source("functions/silhouettes.R")
source("02_data_and_indices.R")
fit <- readRDS("./rds/2022_09_15_fit.rds")


## Now a summary table ---------------------------------------------------------
# The indexes are as follows 
# predators order is coyote, badger, swift fox
# prey order is jackrabbit, cottontail 
# psi covariates are intercept, prairie, vegetation height, grass, forbs
# det covariates are intercept, vegetation height, prairie 

# pred_beta[predator, intercept]
# pred_theta[predator]
# pred_alpha[predator, det covar]
# inxs_b0[prey, predator]
# inxs_beta[prey, predator, psi covar]
# prey_theta[prey]
# prey_beta[prey, psi covar]
# prey_alpha[prey, det covar]

# vectors for each column 
prey_l <- c(rep(' ', 15), rep(c('jackrabbit', 'cottontail'), 24))
pred_l <- c(rep(c('coyote', 'badger', 'swift fox'), 5), 
            rep(c('coyote', 'coyote', 'badger', 'badger', 'swift fox', 'swift fox'), 5), 
            rep(' ', 2), 
            rep('no predators', 10), 
            rep(' ', 6))
covariate_l <- c(rep(' ', 9), 
                 rep('vegHeight', 3), 
                 rep('prairie', 3), 
                 rep(' ', 6),
                 rep('prairie', 6), 
                 rep('vegHeight', 6), 
                 rep('grass', 6),
                 rep('forbs', 6), 
                 rep(' ', 4),
                 rep('prairie', 2), 
                 rep('vegHeight', 2), 
                 rep('grass', 2),
                 rep('forbs', 2),
                 rep(' ', 2), 
                 rep('vegHeight', 2), 
                 rep('prairie', 2))

params_l <- c(rep('$\\beta_0$', 3), rep('$\\theta$', 3), 
              rep('$\\alpha_0$', 3), 
              rep('$\\alpha_1$', 3), 
              rep('$\\alpha_2$', 3), 
              rep('$\\beta_0$ ind x site', 6), 
              rep('$\\beta_1$ ind x site', 6), 
              rep('$\\beta_2$ ind x site', 6), 
              rep('$\\beta_3$ ind x site', 6), 
              rep('$\\beta_4$ ind x site', 6), 
              rep('$\\theta$',2), 
              rep('$\\beta_0$ ind x site', 2), 
              rep('$\\beta_1$ ind x site', 2), 
              rep('$\\beta_2$ ind x site', 2), 
              rep('$\\beta_3$ ind x site', 2), 
              rep('$\\beta_4$ ind x site', 2), 
              rep('$\\alpha_0$', 2), 
              rep('$\\alpha_1$', 2), 
              rep('$\\alpha_2$', 2))

param_list <- data.frame(prey_l, pred_l, covariate_l, params_l)

fit_summary <- MCMCsummary(fit, digits = 3)
fit_summary <- tibble::rownames_to_column(fit_summary, 
                                          "Parameter") %>% 
  rename('Lower 2.5' = 4,
         'Upper 97.5' = 6,
         Mean = 2,
         Median = 5) %>% 
  relocate(Median, .after = Mean) %>% 
  filter(!Parameter == 'deviance') %>% 
  select(!c(Rhat, sd, n.eff))

param_list <- cbind(param_list, fit_summary) %>% 
  select(!Parameter)
colnames(param_list) <- c('Prey', 'Predator', 'Covariate', 'Parameter', 'Mean', 'Median', 'Lower 2.5', 'Upper 97.5')

saveRDS(param_list, './tables/params_table.rds')


## Check that the table looks good 
## Knit it to pdf using an rmarkdown or quarto document
param_list %>% 
  kbl(digits = 2, align = c('l', 'l', 'l', 'c', 'c', 'c', 'c', 'c'), 
      booktabs = TRUE, 
      # format = 'latex',
      longtable = TRUE, 
      # make sure to FALSE so it will show latex notation
      escape = FALSE) %>% 
  kable_classic(full_width=FALSE) %>% 
  pack_rows('Predator occupancy model', 1, 6) %>% 
  pack_rows('Predator detection model', 7, 15) %>% 
  pack_rows('Predator-prey co-occurrence model', 16, 57) %>% 
  pack_rows('Prey detection model', 58, 63) %>% 
  row_spec(0, bold=TRUE, font_size = 18, align='c') %>% 
  add_header_above(c(" " = 6,
                     "95% CI" = 2),
                   font_size=18, 
                   bold=TRUE, 
                   align = 'c')

# 
#   pack_rows('Prey model', 7, 48) %>% 
#   pack_rows('jackrabbit', 7, 27, indent = FALSE) %>% 
#   pack_rows('cottontail', 28, 48, indent = FALSE) %>% 
#   pack_rows('Detection model', 49, 63) %>% 
#   pack_rows('jackrabbit', 49, 51) %>% 
#   pack_rows('cottontail', 52, 54) %>% 
#   pack_rows('coyote', 55, 57) %>% 
#   pack_rows('badger', 58, 60) %>% 
#   pack_rows('swift fox', 61, 63) %>% 

## END ----------------------------------------------------------
