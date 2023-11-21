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
library(ggh4x)

# Source functions and data ----------------------------------------------
source("functions/silhouettes.R")
source("02_data_and_indices.R")
fit <- readRDS("./2022_09_15_fit.rds")

# Get posterior data --------------------------------------------------------
mcmc_list <- bbsBayes::get_mcmc_list(fit) # to see objects in mcmc 
# we see in mcmc_list that there is a list called sims.list which is what we want
tmp <- fit$sims.list # sims.list contains vectorized posterior samples produced by jagsUI
names(tmp)

## Plot occupancy --------------------------------------------------------------
# named vector, elements match column names
cov_order <- c(
  "Prairie (proportion)" = "PrairieSum_1k",
  "Vegetation height (cm)" = "VegHeight",
  "Grass (proportion)" = "GrassPrp",
  "Forb (proportion)" = "ForbPrp"
)

predator_order <- c("coyote", "badger", "swift fox") 

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
  rc.pred <- seq( # prediction for covariate i
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
    # Estimate expected occupancy 
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
prey_df <- map_dfr(prey_ests, ~as.data.frame(.x), .id = "id")

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
  mutate_at(vars(id),
            list(~ factor(., levels = c("btjr", "ectr"))))%>% 
  mutate_at(vars(predator),
            list(~ factor(., levels = c("coyote", "badger", "swiftfox", "nopredators"))))->all.prey.pred.g

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
  geom_line(aes(x=covariate, y=estimate, color = predator), linewidth = 2)+
  geom_ribbon(aes(x = covariate, 
                  ymin = lower95, 
                  ymax = upper95, 
                  fill=predator), alpha = 0.2)+
  ggh4x::facet_grid2(id~predator, 
                     axes = 'x')+
  # facet_grid(id~predator,
  #            scales='free')+
  xlim(0, 0.3)+
  scale_color_manual(values = colors_pred3)+
  scale_fill_manual(values = colors_pred3)+
  labs(x = 'Proportion of Forbs',
       y = 'Occupancy probability')+
  theme_classic()+
  theme(legend.position = "none",
        panel.spacing = unit(0.5, 'cm'),
        strip.text = element_blank(),
        axis.title = element_text(size = 18, face = 'bold'),
        axis.title.x=element_text(vjust = -2),
        axis.title.y=element_text(vjust=2.5),
        axis.text = element_text(size=16, color = 'black'),
        plot.margin =(margin(0,0,0.5,0.5, unit = 'cm')),
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
  annotate('text', x=c(0.1,0.21,0.42,0.61,0.85,0.9), y=1, 
           label=c(' ','coyote', 'badger', 'swift fox', 'no predators', ' '), 
           fontface=2, 
           size=6)->label_sils
  
ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotation_custom(coyote.s, xmin = 0.06, xmax=0.2)+
  annotation_custom(badger.s, xmin = 0.27, xmax = 0.5)+
  annotation_custom(sfox.s, xmin = 0.57, xmax=0.7) +
  annotation_custom(np.s, xmin=0.78, xmax=1.1)-> pred.sil 
pred.sil

ggplot(mapping = aes(x=1, y=0:1))+
  theme_void()+
  annotate('text',x=1, y=c(0.10, 0.2, 0.52, 0.8),  
           label=c(' ', 'cottontail', 'jackrabbit',' '), 
           fontface=2, 
           size=6)+
  annotation_custom(jack.s, ymin=0.51, ymax=0.63) +
  annotation_custom(cotton.s, ymin=0.2, ymax=0.3) -> prey.sil

gridExtra::grid.arrange(pred.sil, label_sils, pplot, heights=c(0.1, 0.05, 0.9))->forbs.sil

gridExtra::grid.arrange(forbs.sil,prey.sil,  widths=c(0.9, 0.10))->forbs.sil

forbs.sil

ggsave(filename = paste0("./plots2/",
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
  covar <- covar[1:3]
  names <- as.factor(c(
    "Proportion of Prairie (buffer = 1k)", 
    "Vegetation height (cm)", 
    "Proportion of Grass", 
    "Proportion of Forbs"
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
      ggh4x::facet_grid2(id~predator, 
                         axes = 'x')+
      # facet_grid(
      #   # data = filter(data, name == covar[i]),
      #   rows = vars(id),
      #   cols = vars(predator))+
      scale_color_manual(values = colors_pred3)+
      scale_fill_manual(values = colors_pred3)+
      labs(
        x = names[i],
        y = 'Occupancy probability')+
      ylim(0,1)+
      theme_classic()+
      theme(legend.position = "none",
            panel.spacing = unit(0.5, 'cm'),
            strip.text = element_blank(),
            axis.title = element_text(size = 18, face = 'bold'),
            axis.title.x=element_text(vjust = -2),
            axis.title.y=element_text(vjust=2.5),
            axis.text = element_text(size=16, color = 'black'),
            plot.margin =(margin(0,0,0.5,0.5, unit = 'cm')),
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
    
    ggplot(mapping = aes(x=0:1, y=1))+
      theme_void()+
      annotate('text', x=c(0.1,0.21,0.42,0.61,0.85,0.9), y=1, 
               label=c(' ','coyote', 'badger', 'swift fox', 'no predators', ' '), 
               fontface=2, 
               size=6)->label_sils
    
    ggplot(mapping = aes(x=0:1, y=1))+
      theme_void()+
      annotation_custom(coyote.s, xmin = 0.06, xmax=0.2)+
      annotation_custom(badger.s, xmin = 0.27, xmax = 0.5)+
      annotation_custom(sfox.s, xmin = 0.57, xmax=0.7) +
      annotation_custom(np.s, xmin=0.78, xmax=1.1)-> pred.sil 
    pred.sil
    
    ggplot(mapping = aes(x=1, y=0:1))+
      theme_void()+
      annotate('text',x=1, y=c(0.10, 0.2, 0.52, 0.8),  
               label=c(' ', 'cottontail', 'jackrabbit',' '), 
               fontface=2, 
               size=6)+
      annotation_custom(jack.s, ymin=0.49, ymax=0.63) +
      annotation_custom(cotton.s, ymin=0.2, ymax=0.3) -> prey.sil
    
    gridExtra::grid.arrange(pred.sil, label_sils, plot, heights=c(0.1, 0.05, 0.9))->plot.sil
    
    gridExtra::grid.arrange(plot.sil,prey.sil,  widths=c(0.9, 0.1))->plot.sil
    
    plot.sil
    
    
    # gridExtra::grid.arrange(pred.sil, plot, heights=c(0.1, 0.9))->plot.sil
    # gridExtra::grid.arrange(plot.sil, prey.sil, widths=c(0.9, 0.1))->plot.sil
    
        # create folder to save the plots to
    if (dir.exists("plots2")) { } 
    else {dir.create("plots2")}
    
    # save plots to the 'figures' folder
    ggsave(filename = paste0("./plots2/",
                             covar[i],
                             "_prey_occ_plot.png"),
           plot = plot.sil,
           width = 18,
           height = 12,
           units = "in")
    
    # print each plot to screen
    print(plot.sil)
  }
  
}

make_plot(all.prey.pred.g)

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
map_dfr(prey_ests, ~as.data.frame(.x), .id = "id") -> prey_df

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
            list(~ factor(., levels = c("VegHeight", "PrairieSum_1k")))) ->det_prey_df
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
                aes(x=covariate, y=estimate, color = id), linewidth = 1.5)+
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

    # print each plot to screen
    print(plot)
  }
}

make_plot_det(det_prey_df)

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
            list(~ factor(., levels = c("coyote", "badger", "swiftfox", "nopredator")))) %>% 
  mutate_at(vars(name),
            list(~ factor(., levels = c("VegHeight", "PrairieSum_1k"))))->det_predator_df

# View(det_prey_df)
# head(det_predator_df)

colors_pred3 <- c("#ee9b00","#0a9396","#bb3e03")
cov_order <- c(
  "Vegetation Height (cm)" = "VegHeight",
  "Prairie (proportion)" = "PrairieSum_1k"
)

## Create function for sequential graphing of data by covariate (e.g., RowcropPrp_1k). 
make_plot_det <- function(data){
  # a vector of names of covariates to loop over 
  covar <- unique(det_predator_df$name)
  names <- as.factor(c("Vegetation Height (cm)" ,
                       "Proportion of prairie (buffer = 1k)"
  ))
  # a loop to produce the rest of the graphs 
  for (i in seq_along(covar)){
    plot <- data %>% 
      ggplot()+
      geom_line(data = filter(data, name == covar[i]),
                aes(x=covariate, y=estimate, color = id), linewidth = 1.5)+
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
    # print each plot to screen
    print(plot)
  }
}

make_plot_det(det_predator_df)

## Let's plot predators and prey detection together in one plot per covar. 
det_all_df <- rbind(det_predator_df, det_prey_df)

## Silhouettes for detection plots all together 
ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotate('text', x=c(0.02, 0.11, 0.29,0.47,0.66,0.86,0.9), y=1, 
           label=c(' ', 'coyote', 'badger', 'swift fox', 'jackrabbit', 'cottontail', ' '), 
           fontface=2, 
           size=6)->sil_lab
sil_lab

ggplot(mapping = aes(x=0:1, y=1))+
  theme_void()+
  annotation_custom(coyote.s, xmin=0.01, xmax=0.19) +
  annotation_custom(badger.s, xmin=0.2, xmax=0.4)+
  annotation_custom(sfox.s, xmin = 0.41, xmax=0.6)+
  annotation_custom(jack.s, xmin = 0.65, xmax = 0.8)+
  annotation_custom(cotton.s, xmin = 0.9, xmax=1) -> det.sil 
det.sil


# Color palette
colors_pp <- c('#5F0F40', '#1e6091', colors_pred3)

## Det Plots Together 
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

vplot <- det_all_df %>%
  filter(name == "VegHeight") %>% 
  ggplot()+
  geom_line(aes(x=covariate, y=estimate, color = id), size = 1.5)+
  geom_ribbon(aes(x = covariate, ymin = lower95, ymax = upper95, fill=id), alpha = 0.2)+
  scale_color_manual(values = colors_pp)+
  scale_fill_manual(values = colors_pp)+
  facet_grid(
    # data = filter(data, name == covar[i]),
    cols = vars(id))+
  labs(x = "Vegetation Height (cm)",
       y = 'Detection probability')+
  scale_y_continuous(breaks=seq(0,1, 0.2), limits = c(0, 1),
                     labels = dropLeadingZero)+
  scale_x_continuous(labels = dropLeadingZero)+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size=20, color = 'black'),
        axis.title.x = element_text(margin = margin(0.5, 0,0.5,0, 'cm')), 
        axis.title.y = element_text(margin=margin(0, 0.5, 0,0.5, 'cm')),
        panel.grid.major.y = element_line(color = "gray80"),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

pplot <- det_all_df %>%
  filter(name == "PrairieSum_1k") %>% 
  ggplot()+
  geom_line(aes(x=covariate, y=estimate, color = id), size = 1.5)+
  geom_ribbon(aes(x = covariate, ymin = lower95, ymax = upper95, fill=id), alpha = 0.2)+
  scale_color_manual(values = colors_pp)+
  scale_fill_manual(values = colors_pp)+
  facet_grid(
    # data = filter(data, name == covar[i]),
    cols = vars(id))+
  labs(x = "Proportion of prairie (buffer = 1k)",
       y = 'Detection probability')+
  scale_y_continuous(breaks=seq(0,1, 0.2), limits = c(0, 1),
                     labels = dropLeadingZero)+
  scale_x_continuous(labels = dropLeadingZero)+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size=20, color = 'black'),
        axis.title.x = element_text(margin = margin(0.5, 0,0.5,0, 'cm')), 
        axis.title.y = element_text(margin=margin(0, 0.5, 0,0.5, 'cm')),
        panel.grid.major.y = element_line(color = "gray80"),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

gridExtra::grid.arrange(det.sil,sil_lab, vplot, heights=c(0.1, 0.05, 0.9)) ->plot01
gridExtra::grid.arrange(det.sil,sil_lab, pplot, heights=c(0.1, 0.05, 0.9)) ->plot02    

gridExtra::grid.arrange(plot01, pplot, ncol=1) -> det_plots

ggsave("./plots2/det_plot.png",
       plot = det_plots,
       width = 20,
       height = 10,
       units = "in")

# END ---------------------------------------------------------------------