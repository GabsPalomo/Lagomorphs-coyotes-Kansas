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

# Source functions and data ----------------------------------------------
source("functions/silhouettes.R")
source("02_data_and_indices.R")
fit <- readRDS("./2022_09_15_fit.rds")

## Plot posterior distributions ----------------------------------------------------
posterior_sims <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(posterior_sims)

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
prey_psi_covar <- c('jackrabbit-intercept', 'jackrabbit-forbs', 'jackrabbit-grass', 'jackrabbit-vegH', 'jackrabbit-prairie',
                    'cottontail-intercept', 'cottontail-forbs', 'cottontail-grass', 'cottontail-vegH', 'cottontail-prairie')

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

tp_rwy <- plot_grid(titleplt)
mp_rwy <- plot_grid(mypreyplots[[1]], mypreyplots[[3]],ncol=2, rel_widths=c(0.4, 0.6))
bt_rwy <- plot_grid(mypreyplots[[2]])  
cowplot::plot_grid(tp_rwy, mp_rwy, bt_rwy, nrow = 3, rel_heights = c(0.1, 1, 1))->prey_plots
prey_plots

## Now we do the interaction term (inxs) part -----------------------------------------------
inxs_list1 <- post_samples[4]
inxs_list2 <- post_samples[5]
inxs_list2_split <- lapply(inxs_list2, function(x) split(x, x$observation))
inxs_list2_split <- unlist(test, recursive = FALSE)
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
plot_grid(pred_plots, prey_plots, plots_inxs, nrow = 1) -> all_plots 
ggsave(filename = 'plots/posteriors_all_plots.jpg', 
       plot=all_plots, 
       dpi=600, 
       width = 20, 
       height = 8, 
       units = c('in'))

## Now a summary ---------------------------------------------------------
fit_summary <- MCMCsummary(fit, digits = 3)
fit_summary <- tibble::rownames_to_column(fit_summary, 
                                          "Parameter") %>% 
  rename( lower2.5 = 4,
          upper97.5 = 6,
          median = 5)

