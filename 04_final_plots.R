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

post_samples <- lapply(post_samples, function(x)
  x %>% 
    pivot_longer(everything(), 
                 names_to = 'observation', 
                 values_to = 'values')
  )

## Let's separate the predator and prey models 
post_sampleS_pred <- post_samples[c(1,2,3)]
post_sampleS_prey <- post_samples[c(6,7,8)]
post_inxs <- post_samples[c(4,5)]

mypredplots <- list()
mypreyplots <- list()

## Predator model 
col_pal1 <- c('#587b7c', '#713e86', '#f57600', 
              '#c9d382', '#7fbed8', '#d8b0cc',
              '#38875f', '#7f7ead', '#886a63',
              '#234c5b', '#49608e', '#64447c')
col_pal2 <- c('#c9d382', '#7fbed8', '#d8b0cc',
              '#38875f', '#7f7ead', '#886a63',
              '#234c5b', '#89d4d1', '#64447c')
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

leg1 <- mypredplots[[1]]+
  theme(legend.position = 'right')+
  ggtitle('Predator model')
titleplot <- get_title(leg1)
leg1 <- get_legend(leg1)
plot(leg1)

leg2 <- mypredplots[[3]]+
  scale_fill_manual(values=col_pal2)+
  theme(legend.position = 'right')
tp_rw <- plot_grid(titleplot)
mp_rw <- plot_grid(mypredplots[[1]], mypredplots[[2]], leg1, ncol=3, rel_widths=c(1, 1, 0.2))
bt_rw <- plot_grid(leg2)  
plot_grid(tp_rw, mp_rw, bt_rw, nrow=3, rel_heights = c(0.1, 1, 1))->pred_plots

## Prey model 
for(i in 1:length(post_sampleS_prey)){
  mypreyplots[[i]] <- ggplot(post_sampleS_prey[[i]], 
                             aes(x=values, fill=observation))+
    geom_density(aes(y = after_stat(density - 0.01 * group)), 
                 alpha=0.3)+
    theme_classic()+
    theme(legend.position = 'right', 
          legend.title = element_blank())+
    ylab('')+
    xlab(names(post_sampleS_prey)[[i]])+
    scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values=col_pal1)
}

tp_rwy <- plot_grid(mypreyplots[[1]], mypreyplots[[3]],ncol=2, rel_widths=c(0.4, 0.6))
bt_rwy <- plot_grid(mypreyplots[[2]])  
cowplot::plot_grid(tp_rwy, bt_rwy, nrow = 2)->prey_plots

## Now we do the inxs part 
inxs_list1 <- post_samples[4]
inxs_list2 <- post_samples[5]
test <- lapply(inxs_list2, function(x) split(x, x$observation))
test <- unlist(test, recursive = FALSE)
nmt <- names(test)
nmt <- str_remove_all(nmt, 'inxs_beta.')
names(test) <- nmt

names_gr <- c("cottontail-badger", "cottontail-coyote", "cottontail-sfox", 
              "jackrabbit-badger", "jackrabbit-coyote", "jackrabbit-sfox")
full_list <- list()
names(full_list) <- names_gr
for(i in 1:length(names_gr)){
  full_list[[i]] <- names(test) %>% 
    str_detect(names_gr[i]) %>% 
    keep(test,.)
}

full_list2 <- list()
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
                                    face='bold'))+
    ylab('')+
    xlab(paste0('individual x site: ', names_gr[[i]]))+
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
          rel_heights = c(0.1, 1, 0.1))

# Get posterior data --------------------------------------------------------
fit_summary <- MCMCsummary(fit, digits = 3)
fit_summary <- tibble::rownames_to_column(fit_summary, 
                                          "Parameter") %>% 
  rename( lower2.5 = 4,
          upper97.5 = 6,
          median = 5)

fit_summary %>% 
  filter(!Parameter == 'deviance') %>% 
  ggplot(aes(x= median, y=Parameter))+
  geom_errorbar(aes(xmin = lower2.5, xmax=upper97.5),
                width = 0, linewidth = 2.5, color = "#92b6b1")+
  geom_point(color= "#666a86", size = 3)+
  geom_vline(xintercept = 0, color = '#b2c9ab')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size=18))+
  theme_hc() -> params.median

# ggsave(plot = params.median,
#        filename = "figures/params.median.png",
#        width = 15, 
#        height = 20, 
#        units = "in",
#        dpi =300
#        )

mcmc_list <- bbsBayes::get_mcmc_list(fit) # to see objects in mcmc 
# we see in mcmc_list that there is a list called sims.list which is what we want
tmp <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(tmp)

## Plot occupancy --------------------------------------------------------------

# named vector, elements match column names
cov_order <- c(
  # "Row crop (proportion)" = "RowcropPrp_1k",
  "Prairie (proportion)" = "PrairieSum_1k",
  "Vegetation height (cm)" = "VegHeight",
  "Grass (proportion)" = "GrassPrp",
  "Forb (proportion)" = "ForbPrp"
  # ,
  # "Rowcrop total edge" = "RowcropTE"
)

predator_order <- c("coyote", "badger", "swift fox") # removed bobcat

my_ranges <- matrix(
  NA,
  ncol = 2, # 2 cols because min and max for the ranges 
  nrow = length(cov_order) # one row for each covariate 
)

# get nice ranges for each
for(i in 1:length(cov_order)){
  if(cov_order[i] %in% colnames(landscape)){
    my_ranges[i,] <- range(landscape[,cov_order[i]],na.rm = TRUE)  
  } else {
    my_ranges[i,] <- range(cov[,cov_order[i]],na.rm = TRUE)  
  }
}
# look at it
my_ranges 

# make 'pretty' ranges for plotting
# This is going to be the x axis so they don't need to be scaled. 
my_ranges <- matrix(
  c(
    0,1,      #cov1
    # 0,1,      
    0,150,
    0,1,
    0,0.7
    # ,
    # 0, 111000
  ),
  ncol = 2,
  nrow = length(cov_order),
  byrow = TRUE
)

# function to scale data
scale_dat <- function(x,y){
  (x - mean(y, na.rm = TRUE)) / sd(y,na.rm = TRUE)
}

# make a list object to store in predictions for each prey species with and
# without each predator. This is a nested list that goes
# 2 prey species > 5 covariates > prey without predators and prey with each predator
prey_ests <- vector("list", length = 2) # this is a list with 2 objects, one for each prey species
names(prey_ests) <- c("btjr", "ectr") # we rename the lists to match each prey
# add in sub-list for btjr with 5 objects which are the number of covariates for psi
prey_ests$btjr <- vector("list", length = length(cov_order))
# rename the sub-lists to match the names of covariates
names(prey_ests$btjr) <- cov_order
#Now with second prey
prey_ests$ectr <- vector("list", length = length(cov_order))
names(prey_ests$ectr) <- cov_order

# Create sublists for the combination of each prey species, 
# covariate, and predators (including no predators)
# I use a for loop to create all the sublists including all the combinations 
# so that I don't have to do them individually
for(i in 1:length(cov_order)){
  prey_ests$btjr[[i]] <- vector(
    "list",
    length = npred + 1 # +1 to include "no predators" 
  )
  names(prey_ests$btjr[[i]]) <- c("no_predators", predator_order)
  prey_ests$ectr[[i]] <- vector(
    "list",
    length = npred + 1
  )
  names(prey_ests$ectr[[i]]) <- c("no_predators", predator_order)
}

# start filling in all of this! 
# This takes time
# I use a for loop so I don't have to fill these lists and sublists individually for all covariates 
# but sort of sequentially and in fewer lines. 
pb <- txtProgressBar(max = length(cov_order))
for(i in 1:length(cov_order)){
  setTxtProgressBar(pb, i) # progress bar for each covariate 
  rc.pred <- seq( # prediction for row crop
    my_ranges[i,1], # min 
    my_ranges[i,2], # max 
    length.out = 200 # length or number of predictions for plot (x)
  )
  # 
  if(cov_order[i] %in% colnames(landscape)){
    rc.pred.scaled <- scale_dat(
      rc.pred,
      landscape[,cov_order[i], drop = TRUE]
    )
  }else{
    rc.pred.scaled <- scale_dat(
      rc.pred,
      cov[,cov_order[i], drop = TRUE]
    )
  }
  pred_mat <- matrix(
    0, # fill them with zeros at first
    ncol = nparm_prey_psi, # one column per covariate (5 + 1 intercept)
    nrow = length(rc.pred)
  )
  # Now we fill pred_mat with data 
  # 1 for intercept
  pred_mat[,1] <- 1
  # add in the covariate at the correct location + 1 because intercept
  pred_mat[,i+1] <- rc.pred.scaled
  # make predictions for each species
  for(j in 1:nprey){
    psi1prey <- plogis(
      tmp$prey_beta[,j,] %*% t(pred_mat)
      # remember tmp$prey_beta is an array with rows=sims, cols=2 prey, matrix=6 covariates (1 is intercept)
      # AND pred_mat has rows = 2 prey and cols = 6 covariates (1 is intercept)
      # so the t before pred_mat transpose it so it has 2 columns just like tmp$prey_beta
    )
    # season 2
    psi2prey <- tmp$prey_beta[,j,] %*% t(pred_mat) 
    # add on theta
    psi2prey <- psi2prey + rep(
      tmp$prey_theta[,j],
      ncol(psi2prey)
    )
    psi2prey <- plogis(
      psi2prey
    )
    eoc <- psi1prey / (psi1prey + (1 - psi2prey))
    prey_ests[[j]][[i]]$no_predators <- data.frame(
      covariate = rc.pred,
      t(
        apply(
          eoc,
          2,
          quantile,
          probs = c(0.025,0.5,0.975)
        )
      )
    )
    colnames(prey_ests[[j]][[i]]$no_predators) <- c(
      "covariate",
      "lower95",
      "estimate",
      "upper95"
    )
    # and now do the predators too!
    for(s in 1:npred){
      # get inxs stuff
      tmp_mcmc <- cbind(
        tmp$inxs_b0[,j,s],
        tmp$inxs_beta[,j,s,]
      )
      psi1prey <- plogis(
        tmp$prey_beta[,j,] %*% t(pred_mat) +
          tmp_mcmc %*% t(pred_mat)
      )
      # season 2
      psi2prey <- tmp$prey_beta[,j,] %*% t(pred_mat) +
        tmp_mcmc %*% t(pred_mat)
      # add on theta
      psi2prey <- psi2prey + rep(
        tmp$prey_theta[,j],
        ncol(psi2prey)
      )
      psi2prey <- plogis(
        psi2prey
      )
      eoc <- psi1prey / (psi1prey + (1 - psi2prey))
      # +1 because first one is no predators
      prey_ests[[j]][[i]][[s+1]] <- 
        data.frame(
          covariate = rc.pred,
          t(
            apply(
              eoc,
              2,
              quantile,
              probs = c(0.025,0.5,0.975)
            )
          )
        )
      colnames(prey_ests[[j]][[i]][[s+1]]) <- c(
        "covariate",
        "lower95",
        "estimate",
        "upper95"
      )
      
    }
  }
}



# plot them all out! One plot for each covariate.


## Plot results using ggplot2 -----------------------------------------------------
## Get data ready to plot 
## Extract data frames from list 
library(tidyverse)
prey_df <- map_dfr(prey_ests, ~as.data.frame(.x), .id = "id")
nmtst <- names(prey_df)

prey_df %>% 
  rename_with(~str_replace(., 'swift.fox', 'swiftfox')) %>%
  rename_with(~str_replace(., 'no_predators', 'nopredators')) %>%
  pivot_longer(cols = !id, 
               names_to = c("name", "predator", "variable"),
               names_pattern = "([A-Za-z0-9_]+)\\.([a-zA-Z0–9]+)\\.([a-zA-Z\\d{2}]+)") %>% 
  group_by(predator, variable) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = variable, 
              values_from = value) %>% 
  select(!row) %>% 
  mutate_at(vars(predator),
            list(~ factor(., levels = c("badger", "coyote", "swiftfox", "nopredators"))))->all.prey.pred.g

colors_pred <- c(
  '#4d194d', '#312244', '#065a82', '#1c7293', '#9eb3c2'
)

colors_pred2 <- c(
  "#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51"
)

colors_pred3 <- c("#e07a5f","#3d405b","#81b29a","#f2cc8f")
colors_pred3 <- c("#f18805","#9631AA","#86bbbd","#8a817c")

# Forbs has a different x scale so I plot it separate 
frbplot <- all.prey.pred.g %>% 
  filter(name == 'ForbPrp') 

ggplot(frbplot) +
  geom_line(aes(x=covariate, y=estimate, color = predator), size = 2)+
  geom_ribbon(aes(x = covariate, 
                  ymin = lower95, 
                  ymax = upper95, 
                  fill=predator), alpha = 0.2)+
  facet_grid(rows = vars(id),
             cols = vars(predator))+
  xlim(0, 0.3)+
  scale_color_manual(values = colors_pred3)+
  scale_fill_manual(values = colors_pred3)+
  labs(x = 'Proportion of Forbs',
       y = 'Occupancy')+
  theme_classic()+
  theme(legend.position = "none",
        panel.spacing = unit(0.5, 'cm'),
        strip.text = element_blank(),
        axis.title = element_text(size = 18, face = 'bold'),
        axis.text = element_text(size=16, color = 'black'),
        panel.grid.major.y = element_line(color = "gray80"),
        strip.background = element_blank())-> pplot
        # strip.text = element_textbox_highlight(
        #   size = 12, 
        #   face = "bold",
        #   r = unit(5, "pt"),
        #   box.color = "gray50",
        #   color = "black",
        #   linetype = 1, 
        #   width = unit(1, "npc"),
        #   halign = 0.5, 
        #   padding = margin(5, 0, 5, 0)
        
pplot

ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotation_custom(badger.s, xmin = 0.05, xmax=0.2)+
  annotation_custom(coyote.s, xmin = 0.3, xmax = 0.5)+
  annotation_custom(sfox.s, xmin = 0.57, xmax=0.7) +
  annotation_custom(np.s, xmin=0.78, xmax=1.1)-> pred.sil 
pred.sil
  
ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotation_custom(jack.s, ymin=0.7, ymax=0.85) +
  annotation_custom(cotton.s, ymin=0.25, ymax=0.4) -> prey.sil

gridExtra::grid.arrange(pred.sil, pplot, heights=c(0.1, 0.9))->forbs.sil
gridExtra::grid.arrange(forbs.sil, prey.sil, widths=c(0.9, 0.1))->forbs.sil

forbs.sil

ggsave(filename = paste0("./plots/",
                         "ForbsPrp_scale",
                         "_prey_occ_plot.png"),
       plot = forbs.sil,
       width = 15,
       height = 10,
       units = "in")


## Create function for sequential graphing of data by covariate (e.g., RowcropPrp_1k). 
make_plot <- function(data){
  # a vector of names of covariates to loop over 
  covar <- unique(all.prey.pred.g$name)
  names <- as.factor(c(
    # "Proportion of Rowcrop (buffer = 1k)", 
    "Proportion of Prairie (buffer = 1k)", 
    "Vegetation height (cm)", 
    "Proportion of Grass", 
    "Proportion of Forbs"
    # ,
    # "Rowcrop total edge"
  ))
  
  # a loop to produce the rest of the graphs 
  for (i in seq_along(covar)){
    plot <- data %>% 
      ggplot()+
      geom_line(data = filter(data, name == covar[i]),
                aes(x=covariate, y=estimate, color = predator), 
                linewidth = 2)+
      geom_ribbon(data = filter(data, name == covar[i]),
                  aes(x = covariate, 
                      ymin = lower95, 
                      ymax = upper95, 
                      fill=predator), 
                  alpha = 0.2)+
      facet_grid(
        # data = filter(data, name == covar[i]),
        rows = vars(id),
        cols = vars(predator))+
      scale_color_manual(values = colors_pred3)+
      scale_fill_manual(values = colors_pred3)+
      labs(
        x = names[i],
        y = 'Occupancy')+
      ylim(0,1)+
      theme_classic()+
      theme(legend.position = "none",
            panel.spacing = unit(0.5, 'cm'),
            strip.text = element_blank(),
            axis.title = element_text(size = 18, face = 'bold'),
            axis.text = element_text(size=16, color='black'),
            panel.grid.major.y = element_line(color = "gray80"),
            strip.background = element_blank()
            # strip.text = element_textbox_highlight(
            #   size = 12, 
            #   face = "bold",
            #   r = unit(5, "pt"),
            #   box.color = "gray50",
            #   color = "black",
            #   linetype = 1, 
            #   width = unit(1, "npc"),
            #   halign = 0.5, 
            #   padding = margin(5, 0, 5, 0)
            #)
      )
    gridExtra::grid.arrange(pred.sil, plot, heights=c(0.1, 0.9))->plot.sil
    gridExtra::grid.arrange(plot.sil, prey.sil, widths=c(0.9, 0.1))->plot.sil
    
        # create folder to save the plots to
    if (dir.exists("plots")) { } 
    else {dir.create("plots")}
    
    # save plots to the 'figures' folder
    ggsave(filename = paste0("./plots/",
                             covar[i],
                             "_prey_occ_plot.png"),
           plot = plot.sil,
           width = 18,
           height = 12,
           units = "in")
    
    # print each plot to screen
    print(plot)
  }
  
}

make_plot(all.prey.pred.g)-> test


## Plot coefficients of each covariate and CI ------------------------------
library(ggthemes)
all.prey.pred.g %>% 
  group_by(id, name, predator) %>% 
  summarise(occupancy = mean(estimate),
            lower = quantile(estimate, c(0.05)),
            upper = quantile(estimate, c(0.95))) %>% 
  mutate(name = case_when(name == 'ForbPrp' ~ 'Prp Forbs', 
                          name == 'GrassPrp' ~ 'Prp Grass', 
                          name == 'PrairieSum_1k' ~ 'Prp Prairie 1k', 
                          name == 'VegHeight' ~ 'Vegetation H'),
         id = case_when(id == 'btjr' ~ 'jackrabbit', 
                          id == 'ectr' ~ 'cottontail'))-> data.mean

colnames(data.mean)
unique(data.mean$name)

colors_var <- c(
  "#4d7c8a","#8C2155","#9B9BC8","#EF8354"
)

data.mean %>% 
  ggplot(aes(x = occupancy))+
  facet_grid(vars(name), vars(id))+
  tidybayes::geom_pointinterval(aes(x = occupancy, 
                                    y = predator, 
                                    xmin = lower, 
                                    xmax = upper,
                         group = id, color = name),
                     fatten_point = 3,
                     linewidth = 3)+
  scale_color_manual(values = colors_var)+
  theme_hc()+
  ylab('')+
  theme(legend.position = 'none',
        strip.background =element_blank(), 
        strip.text = element_text(face = 'bold', size=16), 
        axis.title.x = element_text(face = 'bold', size=16), 
        axis.text = element_text(size=14))
ggsave(filename = './plots/occ_coef.png', 
       last_plot(),
       width = 8,
       height = 8,
       units = "in")


# Plot Detection --------------------------------------------------------------------
# named vector, elements match column names
cov_order <- c(
  "Vegetation Height (cm)" = "VegHeight",
  "Prairie (proportion)" = "PrairieSum_1k"
)

my_ranges <- matrix(
  NA,
  ncol = 2, # 2 cols because min and max for the ranges 
  nrow = length(cov_order) # one row for each covariate 
)

# get nice ranges for each
for(i in 1:length(cov_order)){
  if(cov_order[i] %in% colnames(landscape)){
    my_ranges[i,] <- range(landscape[,cov_order[i]],na.rm = TRUE)  
  } else {
    my_ranges[i,] <- range(cov[,cov_order[i]],na.rm = TRUE)  
  }
}
# look at it
my_ranges 

# make 'pretty' ranges for plotting
# This is going to be the x axis so they don't need to be scaled. 
my_ranges <- matrix(
  c(
    0,150,    # Rowcrop_1k
    0,1     # Prairie_1k
  ),
  ncol = 2,
  nrow = length(cov_order),
  byrow = TRUE
)

# function to scale data
scale_dat <- function(x,y){
  (x - mean(y, na.rm = TRUE)) / sd(y,na.rm = TRUE)
}

# make a list object to store in predictions for each PREY species & covariate. 
# This is a nested list that goes
# 2 prey species > 2 covariates 
prey_ests <- vector("list", length = 2) # this is a list with 2 objects, one for each prey species
names(prey_ests) <- c("btjr", "ectr") # we rename the lists to match each prey
# add in sub-list for btjr with 2 objects which are the number of covariates for p
prey_ests$btjr <- vector("list", length = length(cov_order))
# rename the sub-lists to match the names of covariates
names(prey_ests$btjr) <- cov_order
#Now with second prey
prey_ests$ectr <- vector("list", length = length(cov_order))
names(prey_ests$ectr) <- cov_order


# start filling in all of this! 
# This takes time
# I use a for loop so I don't have to fill these lists individually for all covariates 
# but sort of sequentially and in fewer lines. 
pb <- txtProgressBar(max = length(cov_order))
for(i in 1:length(cov_order)){
  setTxtProgressBar(pb, i) # progress bar for each covariate 
  rc.pred <- seq( # prediction for row crop
    my_ranges[i,1], # min 
    my_ranges[i,2], # max 
    length.out = 200 # length or number of predictions for plot (x)
  )
  # Check to see if the covariates I'm using are scaled or not
  if(cov_order[i] %in% colnames(landscape)){
    rc.pred.scaled <- scale_dat(
      rc.pred,
      landscape[,cov_order[i], drop = TRUE]
    )
  }else{
    rc.pred.scaled <- scale_dat(
      rc.pred,
      cov[,cov_order[i], drop = TRUE]
    )
  }
  
 
  # This is the prediction matrix that we will fill with estimate values
  # and place in each sublist. 
  pred_mat <- matrix(
    0, # fill them with zeros at first
    ncol = nparm_prey_rho, # one column per covariate (2 + 1 intercept)
    nrow = length(rc.pred)
  )
  # Now we fill pred_mat with data 
  # 1 for intercept
  pred_mat[,1] <- 1
  # add in the covariate at the correct location + 1 because intercept
  pred_mat[,i+1] <- rc.pred.scaled
  # make predictions for each prey species
  for(j in 1:nprey){
    det1prey <- plogis(
      tmp$prey_alpha[,j,] %*% t(pred_mat)
      # remember tmp$prey_alpha is an array with rows=sims, cols=2 prey, matrix=3 covariates (1 is intercept)
      # AND pred_mat has rows = sims and cols = 3 covariates (1 is intercept)
      # so the t before pred_mat transpose it so it has 2 columns just like tmp$prey_alpha
    )
    # +1 because first one is no predators
    prey_ests[[j]][[i]] <- 
      data.frame(
        covariate = rc.pred,
        t(
          apply(
            det1prey,
            2,
            quantile,
            probs = c(0.025,0.5,0.975)
          )
        )
      )
    colnames(prey_ests[[j]][[i]]) <- c(
      "covariate",
      "lower95",
      "estimate",
      "upper95"
    )
  }
  
}


# Now we plot them all out! One plot for each covariate.

## Plot results using ggplot2 -----------------------------------------------------
#├ Prey detection plots --------------------------------------------------------
## Get data ready to plot 
## Extract data frames from list 
library(tidyverse)
map_dfr(prey_ests, ~as.data.frame(.x), .id = "id") -> prey_df
nm_preydf <- names(prey_df)

prey_df %>% 
  pivot_longer(cols = !id, 
               names_to = c("name", "variable"),
               names_pattern = "([A-Za-z0-9_]+)\\.([a-zA-Z\\d{2}]+)") %>% 
  group_by(variable) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = variable, 
              values_from = value) %>% 
  select(!row) %>% 
  mutate_at(vars(id),
            list(~ factor(., levels = c("btjr", "ectr")))) %>% 
  mutate_at(vars(name),
            list(~ factor(., levels = c("VegHeight", "PrairieSum_1k"))))->det_prey_df

# View(det_prey_df)



colors_prey <- c(
  '#001219', '#005f73' 
)

## Create function for sequential graphing of data by covariate (e.g., RowcropPrp_1k). 
make_plot_det <- function(data){
  # a vector of names of covariates to loop over 
  covar <- unique(det_prey_df$name)
  names <- as.factor(c("Vegetation Height (cm)", 
                       "Proportion of prairie (buffer = 1k)"
  ))
  # a loop to produce the rest of the graphs 
  for (i in seq_along(covar)){
    plot <- data %>% 
      ggplot()+
      geom_line(data = filter(data, name == covar[i]),
                aes(x=covariate, y=estimate, color = id), size = 1.5)+
      geom_ribbon(data = filter(data, name == covar[i]),
                  aes(x = covariate, ymin = lower95, 
                      ymax = upper95, fill=id), alpha = 0.2)+
      scale_color_manual(values = colors_prey)+
      scale_fill_manual(values = colors_prey)+
      facet_grid(
        # data = filter(data, name == covar[i]),
        cols = vars(id))+
      labs(x = names[i],
           y = 'Detection')+
      theme_classic()+
      theme(legend.position = "none",
            panel.spacing = unit(0.5, 'cm'),
            strip.text = element_blank(),
            axis.title = element_text(size = 18, face = 'bold'),
            axis.text = element_text(size=16, color = 'black'),
            panel.grid.major.y = element_line(color = "gray80"),
            strip.background = element_blank())
    
    
    # create folder to save the plots to
    if (dir.exists("figures_06")) { } 
    else {dir.create("figures_06")}
    
    # save plots to the 'figures' folder
    ggsave(filename = paste0("figures_06/",
                             covar[i],
                             "_prey_det_plot.png"),
           plot = plot,
           width = 15, 
           height = 10, 
           units = "in")
    
    # print each plot to screen
    print(plot)
  }
}

make_plot_det(det_prey_df)
names(det_prey_df)
head(det_prey_df)

#├ Predator detection plots --------------------------------------------------------
##  CURRENTLY I DIDN'T INCLUDE ANY DETECTION COVARIATES FOR THE PREDATORS
## SO THIS CANNOT BE COMPUTED. I CAN ALWAYS CHANGE IT. 
# named vector, elements match column names
cov_order <- c(
  "Vegetation Height (cm)" = "VegHeight",
  "Prairie (proportion)" = "PrairieSum_1k"
)


my_ranges <- matrix(
  NA,
  ncol = 2, # 2 cols because min and max for the ranges 
  nrow = length(cov_order) # one row for each covariate 
)

# get nice ranges for each
for(i in 1:length(cov_order)){
  if(cov_order[i] %in% colnames(landscape)){
    my_ranges[i,] <- range(landscape[,cov_order[i]],na.rm = TRUE)  
  } else {
    my_ranges[i,] <- range(cov[,cov_order[i]],na.rm = TRUE)  
  }
}
# look at it
my_ranges 

# make 'pretty' ranges for plotting
# This is going to be the x axis so they don't need to be scaled. 
my_ranges <- matrix(
  c(
    0,150,    # Veg Height
    0,1     # Prairie_1k
  ),
  ncol = 2,
  nrow = length(cov_order),
  byrow = TRUE
)

# function to scale data
scale_dat <- function(x,y){
  (x - mean(y, na.rm = TRUE)) / sd(y,na.rm = TRUE)
}

# make a list object to store in predictions for each PREDATOR species & covariate. 
# This is a nested list that goes
# 2 PREDATOR species > 2 covariates 
pred_ests <- vector("list", length = 3) # this is a list with 2 objects, one for each predator species
names(pred_ests) <- predator_order # we rename the lists to match each prey
# add in sub-list for btjr with 2 objects which are the number of covariates for p
pred_ests$coyote <- vector("list", length = length(cov_order))
# rename the sub-lists to match the names of covariates
names(pred_ests$coyote) <- cov_order
#Now with second predator
pred_ests$badger <- vector("list", length = length(cov_order))
names(pred_ests$badger) <- cov_order
#Now with third predator
pred_ests$'swift fox' <- vector("list", length = length(cov_order))
names(pred_ests$'swift fox') <- cov_order

# start filling in all of this! 
# This takes time
# I use a for loop so I don't have to fill these lists individually for all covariates 
# but sort of sequentially and in fewer lines. 
pb <- txtProgressBar(max = length(cov_order))
for(i in 1:length(cov_order)){
  setTxtProgressBar(pb, i) # progress bar for each covariate 
  rc.pred <- seq( # prediction for row crop
    my_ranges[i,1], # min 
    my_ranges[i,2], # max 
    length.out = 200 # length or number of predictions for plot (x)
  )
  # Check to see if the covariates I'm using are scaled or not
  if(cov_order[i] %in% colnames(landscape)){
    rc.pred.scaled <- scale_dat(
      rc.pred,
      landscape[,cov_order[i], drop = TRUE]
    )
  }else{
    rc.pred.scaled <- scale_dat(
      rc.pred,
      cov[,cov_order[i], drop = TRUE]
    )
  }
  # This is the prediction matrix that we will fill with estimate values
  # and place in each sublist. 
  pred_mat <- matrix(
    0, # fill them with zeros at first
    ncol = nparm_pred_rho, # one column per covariate (2 + 1 intercept)
    nrow = length(rc.pred)
  )
  # Now we fill pred_mat with data 
  # 1 for intercept
  pred_mat[,1] <- 1
  # add in the covariate at the correct location + 1 because intercept
  pred_mat[,i+1] <- rc.pred.scaled
  # make predictions for each prey species
  for(j in 1:npred){
    det1pred <- plogis(
      tmp$pred_alpha[,j,] %*% t(pred_mat)
      # remember tmp$prey_alpha is an array with rows=sims, cols=2 prey, matrix=3 covariates (1 is intercept)
      # AND pred_mat has rows = sims and cols = 3 covariates (1 is intercept)
      # so the t before pred_mat transpose it so it has 2 columns just like tmp$prey_alpha
    )
    # +1 because first one is no predators
    pred_ests[[j]][[i]] <- 
      data.frame(
        covariate = rc.pred,
        t(
          apply(
            det1pred,
            2,
            quantile,
            probs = c(0.025,0.5,0.975)
          )
        )
      )
    colnames(pred_ests[[j]][[i]]) <- c(
      "covariate",
      "lower95",
      "estimate",
      "upper95"
    )
  }
  
}


# Now we plot them all out! One plot for each covariate.

## Get data ready to plot 
## Extract data frames from list 
map_dfr(pred_ests, ~as.data.frame(.x), .id = "id") -> pred_df
nm_preddf <- names(pred_df)

pred_df %>% 
  pivot_longer(cols = !id, 
               names_to = c("name", "variable"),
               names_pattern = "([A-Za-z0-9_]+)\\.([a-zA-Z\\d{2}]+)") %>% 
  group_by(variable) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = variable, 
              values_from = value) %>% 
  select(!row) %>% 
  mutate_at(vars(id),
            list(~ factor(., levels = predator_order))) %>% 
  mutate_at(vars(name),
            list(~ factor(., levels = c("VegHeight", "PrairieSum_1k"))))->det_predator_df



# View(det_prey_df)
head(det_predator_df)

colors_pred3 <- c("#ee9b00","#0a9396","#bb3e03")
cov_order <- c(
  "Vegetation Height (cm)" = "VegHeight",
  "Prairie (proportion)" = "PrairieSum_1k"
)



## Create function for sequential graphing of data by covariate (e.g., RowcropPrp_1k). 
make_plot_det <- function(data){
  # a vector of names of covariates to loop over 
  covar <- unique(det_prey_df$name)
  names <- as.factor(c("Vegetation Height (cm)" ,
                       "Proportion of prairie (buffer = 1k)"
  ))
  # a loop to produce the rest of the graphs 
  for (i in seq_along(covar)){
    plot <- data %>% 
      ggplot()+
      geom_line(data = filter(data, name == covar[i]),
                aes(x=covariate, y=estimate, color = id), size = 1.5)+
      geom_ribbon(data = filter(data, name == covar[i]),
                  aes(x = covariate, 
                      ymin = lower95, 
                      ymax = upper95, 
                      fill=id), alpha = 0.2)+
      scale_color_manual(values = colors_pred3)+
      scale_fill_manual(values = colors_pred3)+
      facet_grid(
        # data = filter(data, name == covar[i]),
        cols = vars(id))+
      labs(x = names[i],
           y = 'detection')+
      theme_classic()+
      theme(legend.position = "none",
            panel.spacing = unit(0.5, 'cm'),
            strip.text = element_blank(),
            axis.title = element_text(size = 18, face = 'bold'),
            axis.text = element_text(size=16, color = 'black'),
            panel.grid.major.y = element_line(color = "gray80"),
            strip.background = element_blank())
    
    
    # create folder to save the plots to
    if (dir.exists("plots")) { } 
    else {dir.create("plots")}
    
    # save plots to the 'figures' folder
    # ggsave(filename = paste0("./plots/",
    #                          covar[i],
    #                          "_predator_det_plot.png"),
    #        plot = plot,
    #        width = 15, 
    #        height = 10, 
    #        units = "in")
    
    # print each plot to screen
    print(plot)
  }
}

make_plot_det(det_predator_df)

## Let's plot predators and prey detection together in one plot per covar. 
head(det_prey_df, 2)
head(det_predator_df, 2)
det_all_df <- rbind(det_prey_df, det_predator_df)


## Silhouettes for detection plots all together 
ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotation_custom(jack.s, xmin=0.02, xmax=0.18) +
  annotation_custom(cotton.s, xmin=0.25, xmax=0.32)+
  annotation_custom(coyote.s, xmin = 0.4, xmax=0.6)+
  annotation_custom(badger.s, xmin = 0.65, xmax = 0.8)+
  annotation_custom(sfox.s, xmin = 0.85, xmax=1) -> det.sil 
det.sil

# Color palette
colors_pp <- c('#5F0F40', '#1e6091', colors_pred3)

## Create function for sequential graphing of data by covariate (e.g., RowcropPrp_1k). 
make_plot_det <- function(data){
  # a vector of names of covariates to loop over 
  covar <- unique(det_prey_df$name)
  names <- as.factor(c("Vegetation Height (cm)", 
                       "Proportion of prairie (buffer = 1k)"
  ))
  # a loop to produce the rest of the graphs 
  for (i in seq_along(covar)){
    #Function to only use 2 decimal places on the y axis 
    scaleFUN <- function(x) sprintf("%.1f", x)
    # FUnction to delete leading zero in the x axis
    dropLeadingZero <- function(l){
      lnew <- c()
      for(i in l){
        if(i==0){ #zeros stay zero
          lnew <- c(lnew,"0")
        } else if (i>1){ #above one stays the same
          lnew <- c(lnew, as.character(i))
        } else
          lnew <- c(lnew, gsub("(?<![0-9])0+", "", i, perl = TRUE))
      }
      as.character(lnew)
    }
    
    plot <- data %>% 
      ggplot()+
      geom_line(data = filter(data, name == covar[i]),
                aes(x=covariate, y=estimate, color = id), size = 1.5)+
      geom_ribbon(data = filter(data, name == covar[i]),
                  aes(x = covariate, ymin = lower95, ymax = upper95, fill=id), alpha = 0.2)+
      scale_color_manual(values = colors_pp)+
      scale_fill_manual(values = colors_pp)+
      facet_grid(
        # data = filter(data, name == covar[i]),
        cols = vars(id))+
      labs(x = names[i],
           y = 'Detection')+
      scale_y_continuous(breaks=seq(0,1, 0.2), limits = c(0, 1),
                         labels = dropLeadingZero)+
      scale_x_continuous(labels = dropLeadingZero)+
      theme_classic()+
      theme(legend.position = "none",
            axis.title = element_text(size = 26, face = 'bold'),
            axis.text = element_text(size=20, color = 'black'),
            panel.grid.major.y = element_line(color = "gray80"),
            strip.background = element_blank(),
            strip.text = element_blank()
      )
    
    gridExtra::grid.arrange(det.sil, plot, heights=c(0.1, 0.9))->plot
    
    # create folder to save the plots to
    if (dir.exists("plots")) { } 
    else {dir.create("plots")}
    
    # save plots to the 'figures' folder
    ggsave(filename = paste0("./plots/",
                             covar[i],
                             "_pp_det_plot.png"),
           plot = plot,
           width = 20, 
           height = 10, 
           units = "in")
    
    # print each plot to screen
    print(plot)
  }
}

make_plot_det(det_all_df)

# END ---------------------------------------------------------------------