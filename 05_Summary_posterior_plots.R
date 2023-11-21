## Final plots 
## Gabby Palomo and Mason Fidino 
## These plots include 4 covariates for psi and 2 for detection
## 3 predators (coyote, badger, swift fox) and 2 prey (btjr and ectr)

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
fit <- readRDS("./2022_09_15_fit.rds")

## Plot posterior distributions ----------------------------------------------------
posterior_sims <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(posterior_sims)

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

# Let's convert each matrix into a df 
post_samples <- lapply(posterior_sims, function(x) as.data.frame(x))
post_samples <- post_samples[-9] #Remove deviance 
# Let's change the names of each column 
preds <- c('coyote', 'badger', 'sfox')
preys <- c('jackrabbit', 'cottontail')
preds_det_covar <- c('coyote-intercept', 'coyote-vegH', 'coyote-prairie', 
                     'badger-intercept', 'badger-vegH', 'badger-prairie', 
                     'sfox-intercept', 'sfox-vegH', 'sfox-prairie')
preys_det_covar <- c('jackrabbit-intercept', 'jackrabbit-vegH', 'jackrabbit-prairie', 
                     'cottontail-intercept', 'cottontail-vegH', 'cottontail-prairie')
b0prey_pred <- c('jackrabbit-coyote', 'jackrabbit-badger', 'jackrabbit-sfox',
                 'cottontail-coyote', 'cottontail-badger', 'cottontail-sfox')
prey_psi_covar <- c('jackrabbit-intercept', 'jackrabbit-prairie', 'jackrabbit-vegH', 'jackrabbit-grass', 'jackrabbit-forbs',
                    'cottontail-intercept', 'cottontail-prairie', 'cottontail-vegH', 'cottontail-grass', 'cottontail-forbs')

b_is <- c('jackrabbit-coyote-forbs','jackrabbit-badger-forbs', 'jackrabbit-sfox-forbs',
          'jackrabbit-coyote-grass', 'jackrabbit-badger-grass', 'jackrabbit-sfox-grass',
          'jackrabbit-coyote-vegH', 'jackrabbit-badger-vegH', 'jackrabbit-sfox-vegH',
          'jackrabbit-coyote-prairie', 'jackrabbit-badger-prairie', 'jackrabbit-sfox-prairie',
          'cottontail-coyote-forbs','cottontail-badger-forbs', 'cottontail-sfox-forbs',
          'cottontail-coyote-grass', 'cottontail-badger-grass', 'cottontail-sfox-grass',
          'cottontail-coyote-vegH', 'cottontail-badger-vegH', 'cottontail-sfox-vegH',
          'cottontail-coyote-prairie', 'cottontail-badger-prairie', 'cottontail-sfox-prairie')

names(post_samples$pred_beta) <- preds
names(post_samples$pred_theta) <- preds
names(post_samples$pred_alpha) <- preds_det_covar
names(post_samples$prey_alpha) <- preys_det_covar
names(post_samples$prey_theta) <- preys
names(post_samples$inxs_b0) <- b0prey_pred
names(post_samples$inxs_beta) <- b_is
names(post_samples$prey_beta) <- prey_psi_covar

## The df has the data in wider format, we have to make it longer so we can 
## plot everything 
post_samples <- lapply(post_samples, function(x)
  x %>% 
    pivot_longer(everything(), 
                 names_to = 'observation', 
                 values_to = 'values')
)

# We want to create separate sets of graphs for each of the components of the model 
## Let's separate the predator and prey models and the interaction term
post_sampleS_pred <- post_samples[c(1,2,3)]
post_sampleS_prey <- post_samples[c(6,7,8)]
post_inxs <- post_samples[c(4,5)]

## Predator model -------------------------------------------------------------
## Let's start with empty lists to fill in with plots 
mypredplots <- list() # predators 
mypreyplots <- list() # preys 

# Let's set the color palettes to use 
col_pal1 <- c('#587b7c', '#713e86', '#f57600', 
              '#c9d382', '#7fbed8', '#d8b0cc',
              '#38875f', '#7f7ead', '#886a63',
              '#234c5b', '#49608e', '#64447c')
col_pal2 <- c('#c9d382', '#7fbed8', '#d8b0cc',
              '#38875f', '#7f7ead', '#886a63',
              '#234c5b', '#89d4d1', '#64447c')

## Now we create a for loop to create the plots for each part of the predator model
for(i in 1:length(post_sampleS_pred)){
  mypredplots[[i]] <- ggplot(post_sampleS_pred[[i]], 
                             aes(x=values, fill=observation))+
    geom_density(aes(y = after_stat(density - 0.01 * group)), 
                 alpha=0.3)+
    theme_classic()+
    theme(legend.position = 'none', 
          legend.title = element_blank(), 
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16),
          plot.title = element_text(hjust = 0, size = 18, face='bold'))+
    ylab('')+
    xlab(names(post_sampleS_pred)[[i]])+
    scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values=col_pal1)
}

# Extract the legend from the first plot 
leg1 <- mypredplots[[1]]+
  theme(legend.position = 'right')+
  ggtitle('Predator model')
lgnd <- get_legend(leg1)
plot(lgnd)
# Extract the title 
titleplot <- get_title(leg1)

leg2 <- mypredplots[[3]]+
  scale_fill_manual(values=col_pal2)+
  theme(legend.position = 'right')
tp_rw <- plot_grid(titleplot)
mp_rw <- plot_grid(mypredplots[[1]], mypredplots[[2]], lgnd, ncol=3, rel_widths=c(1, 1, 0.2))
bt_rw <- plot_grid(leg2)  
plot_grid(tp_rw, mp_rw, bt_rw, nrow=3, rel_heights = c(0.1, 1, 1))->pred_plots
pred_plots

## Prey model -----------------------------------------------------------------
for(i in 1:length(post_sampleS_prey)){
  mypreyplots[[i]] <- ggplot(post_sampleS_prey[[i]], 
                             aes(x=values, fill=observation))+
    geom_density(aes(y = after_stat(density - 0.001 * group)), 
                 alpha=0.3)+
    theme_classic()+
    theme(legend.position = 'right', 
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16),
          plot.title = element_text(hjust = 0, size = 18, face='bold'))+
    ylab('')+
    xlab(names(post_sampleS_prey)[[i]])+
    scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values=col_pal1)
}


# Extract the legend from the first plot 
lega <- mypreyplots[[1]]+
  theme(legend.position = 'right')+
  ggtitle('Prey model')
# Extract the title 
titleplt <- get_title(lega)

## We just want [[1]] and [[3]]
tp_rwy <- plot_grid(titleplt)
mp_rwy <- plot_grid(mypreyplots[[1]], mypreyplots[[3]],ncol=2, rel_widths=c(0.4, 0.6))
bt_rwy <- plot_grid(mypreyplots[[2]])  
cowplot::plot_grid(tp_rwy, mp_rwy, bt_rwy, nrow = 3, rel_heights = c(0.1, 1, 1))->prey_plots
prey_plots

## Now we do the interaction term (inxs) part -----------------------------------------------
inxs_list1 <- post_samples[4]
inxs_list2 <- post_samples[5]
inxs_list2_split <- lapply(inxs_list2, function(x) split(x, x$observation))
inxs_list2_split <- unlist(inxs_list2_split, recursive = FALSE)
nmt <- names(inxs_list2_split)
nmt <- str_remove_all(nmt, 'inxs_beta.')
names(inxs_list2_split) <- nmt

names_gr <- c("cottontail-badger", "cottontail-coyote", "cottontail-sfox", 
              "jackrabbit-badger", "jackrabbit-coyote", "jackrabbit-sfox")
full_list <- vector("list", length = length(names_gr))
names(full_list) <- names_gr

for(i in 1:length(names_gr)){
  full_list[[i]] <- names(inxs_list2_split) %>% 
    str_detect(names_gr[i]) %>% 
    keep(inxs_list2_split,.)
}

full_list2 <-  vector("list", length = length(names_gr))
names(full_list2) <- names_gr
for(i in 1:length(full_list)){
  full_list2[[i]] <- do.call(rbind, c(full_list[[i]], make.row.names=FALSE))
}


myinxsplots <- list()
for(i in 1:length(full_list2)){
  myinxsplots[[i]] <- ggplot(full_list2[[i]], 
                             aes(x=values, fill=observation))+
    geom_density(aes(y = after_stat(density - 0.01 * group)), 
                 alpha=0.3)+
    theme_classic()+
    theme(legend.position = 'none', 
          legend.title = element_blank(), 
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16),
          plot.title = element_text(hjust = 0, 
                                    size = 18, 
                                    face='bold'), 
          plot.margin=margin(0, 1, 0,0, 'cm'))+
    ylab('')+
    xlab(paste0('ind x site: ', names_gr[[i]]))+
    scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values=col_pal1, labels=c('Forbs', 'Grass', 'Prairie', 'Veg Height'))
}

plot1 <- myinxsplots[[1]]+
  theme(legend.position = 'top')+
  ggtitle('Individuals x site model')
titleplot <- get_title(plot1)
lgnd <- get_legend(plot1)
trow <- plot_grid(titleplot)
mrow <- plot_grid(plotlist = myinxsplots)
brow <- plot_grid(lgnd)
plot_grid(trow, mrow, brow, nrow = 3, 
          rel_heights = c(0.1, 1, 0.1)) -> plots_inxs
plots_inxs
## Now let's put all plots together -------------------------------------------
col1_plots <- align_plots(pred_plots, plots_inxs, align = 'v', axis = 1)
top_row <- plot_grid(top_plots[[1]], prey_plots)
plot_grid(top_row, col1_plots[[2]], 
          nrow = 2, ncol = 1) -> all_plots 
all_plots

ggsave(filename = 'plots/posteriors_all_plots2.jpg', 
       plot=all_plots, 
       dpi=600, 
       width = 16, 
       height = 12, 
       units = c('in'))

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
  select(!c(Rhat, sd))

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
