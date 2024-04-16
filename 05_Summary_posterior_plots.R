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
# source("02_data_and_indices.R")
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


## Estimate proportion of a given effect

fit_effect <- do.call('rbind', fit$samples)
nms <- colnames(fit_effect)[1:63]

fit_effect_df <- fit_effect |> 
  data.frame() |> 
  pivot_longer(everything(), 
               names_to = 'parameter', 
               values_to = 'samples') |> 
  filter(parameter != 'deviance')

nms_df <- unique(fit_effect_df$parameter)

# merge with original names list 
nms_df <- cbind(param_list, nms_df)

# estimate proportion of effect based on direction 
proportion_effect <- fit_effect_df |> 
  group_by(parameter) |> 
  summarise(median = round(median(samples), 4), 
            prop_greater = round(sum(samples>0, na.rm = TRUE)/n(), 2), 
            prop_less = round(sum(samples<0, na.rm = TRUE)/n(), 2)) |>  
  mutate(direction_effect = case_when(median > 0 ~ 'positive', 
                                      median < 0 ~ 'negative', 
                                      median == 0 ~ 'zero')) |> 
  mutate('Pr(direction_effect)' = case_when(median > 0 ~ prop_greater, 
                                median < 0 ~ prop_less, 
                                median == 0 ~ 0)) |> 
  arrange(factor(parameter, levels = nms))

# merge names and proportions 
proportion_effect <- left_join(nms_df, proportion_effect, by = join_by('nms_df' == 'parameter'))


## Now we add the values of each parameter 
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

param_list <- cbind(proportion_effect, fit_summary)

# reorder list 
param_list <- param_list |> 
  select(!c(nms_df, median, prop_greater, prop_less, Parameter)) |> 
  relocate(direction_effect, .after = 'Median') |> 
  relocate('Pr(direction_effect)', .after = direction_effect) |> 
  relocate(params_l)

colnames(param_list) <- c('Parameter','Prey', 'Predator', 'Covariate',  'Mean', 'Median', 'Direction effect', 'Pr(direction effect)', 'Lower 2.5', 'Upper 97.5')

param_list <- rowid_to_column(param_list, var = 'id')

# saveRDS(param_list, './tables/params_table.rds')

# Let's reorganize some rows so the info is more clear and understandable 
param_list_test1 <- param_list |> # 1
  slice(1:6)

param_list_test2 <- param_list |> #2
  slice(46:47)

param_list_test3 <- param_list |> #3
  slice(16:21)

param_list_test31 <- param_list |> #3.1
  slice(48:49)

param_list_test4 <- param_list |> #4
  slice(50:57)

param_list_test5 <- param_list |> #5
  slice(22:45)

param_list_test6 <- param_list |> #6
  slice(58:63)

param_list_test7 <- param_list |> #7
  slice(7:15)

param_list_beta <- rbind(param_list_test3, param_list_test31)
param_list_beta <- param_list_beta |> 
  mutate(Prey = factor(Prey, levels = c('jackrabbit', 'cottontail')), 
         Predator = factor(Predator, levels = c('no predators', 'coyote', 'badger', 'swift fox'))) |> 
  arrange(Prey, Predator)
  
param_list_coo <- rbind(param_list_test4, param_list_test5)
param_list_coo <- param_list_coo |> 
  mutate(Prey = factor(Prey, levels = c('jackrabbit', 'cottontail')), 
         Covariate = factor(Covariate, levels = c('prairie', 'vegHeight', 'grass', 'forbs'))) |> 
  arrange(Prey, Covariate)

param_list_final <- rbind(param_list_test1, 
                          param_list_test2, 
                          param_list_beta,
                          param_list_coo,
                          param_list_test6, 
                          param_list_test7)

param_list_final <- select(param_list_final, -id)

## Check that the table looks good 
## Knit it to pdf using an rmarkdown or quarto document
param_list_final %>% 
  kbl(digits = 2, 
      align = c('l', 'l', 'l', 'c', 'c', 'c', 'c', 'c', 'c', 'c'), 
      booktabs = TRUE,
      row.names = FALSE,
      # format = 'latex',
      longtable = TRUE, 
      # make sure to FALSE so it will show latex notation
      escape = FALSE) %>% 
  kable_classic(full_width=FALSE) %>% 
  pack_rows('Predator occupancy model', 1, 3) %>% 
  pack_rows('Autologistic term', 4, 8) %>% 
  pack_rows('Predator-prey co-occurrence model', 9, 48) %>% 
  row_spec(12, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(16, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(20, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(24, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(28, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(32, extra_css = "border-bottom: 2px solid")  %>% 
  row_spec(36, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(40, extra_css = "border-bottom: 1px solid")  %>% 
  row_spec(44, extra_css = "border-bottom: 1px solid")  %>% 
  pack_rows('Prey detection model', 49, 54) %>% 
  pack_rows('Predator detection model', 55, 63) %>% 
  row_spec(0, bold=TRUE, font_size = 18, align='c') %>% 
  add_header_above(c(" " = 8,
                     "95% CI" = 2),
                   font_size=18, 
                   bold=TRUE, 
                   align = 'c') -> param.table

param.table

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
