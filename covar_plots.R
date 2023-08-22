## Gabriela Palomo
## 2022/02/25

## Create areas to estimate temporal activity per area

library(tidyverse)
library(ggdist)
library(ggthemes)
source("01_data_org.R")
## Site needs to be a character 

##Comparing all years per variable

landscape %>% 
  select(Site, RowcropPrp_1k, PrairieSum_1k) %>% 
  pivot_longer(cols = ends_with("1k"), 
               names_to = "covariate",
               values_to = "proportion") %>%  
  mutate(covariate= factor(covariate)) %>% 
  ggplot(aes(x=covariate, 
             y =proportion,
             color = covariate,
             fill = covariate))  +
  ggdist::stat_halfeye(
    height = -0.5, 
    width = .3, 
    justification = -.7, 
    point_colour = NA) + 
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_point(
    size = 1.3,
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  scale_fill_manual(values = c("#156064", "#fb8f67", "#c9ada7"))+
  scale_color_manual(values = c("#156064", "#fb8f67", "#c9ada7"))+
  xlab("")+
  theme_hc() +
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 45),
        axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold", size = 14))-> g1

g1

cov %>% 
  group_by(SurveyYear) %>% 
  select(SurveyYear, VegHeight) %>% 
  pivot_longer(cols = VegHeight,
               names_to = "covariate",
               values_to = "Height_cm") %>%  
  mutate(covariate= factor(covariate)) %>% 
  mutate(SurveyYear = factor(SurveyYear)) %>% 
  ggplot(aes(x=SurveyYear, 
             y = Height_cm,
             color = SurveyYear,
             fill = SurveyYear)) + 
  ggdist::stat_halfeye(
    height = -0.5, 
    width = .3, 
    justification = -.7, 
    point_colour = NA) + 
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_point(
    size = 1.3,
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  scale_fill_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  scale_color_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  xlab("")+
  theme_hc() +
  theme(legend.position = "none",
        axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold", size = 14))-> g2.1
g2.1


cov %>% 
  group_by(SurveyYear) %>% 
  select(SurveyYear, GrassPrp) %>% 
  pivot_longer(cols = GrassPrp,
               names_to = "covariate",
               values_to = "proportion") %>%  
  mutate(covariate= factor(covariate)) %>% 
  mutate(SurveyYear = factor(SurveyYear)) %>% 
  ggplot(aes(x=SurveyYear, 
             y = proportion,
             color = SurveyYear,
             fill = SurveyYear)) + 
  ggdist::stat_halfeye(
    height = -0.5, 
    width = .3, 
    justification = -.7, 
    point_colour = NA) + 
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_point(
    size = 1.3,
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  ylab("Proportion of grass")+
  scale_fill_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  scale_color_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  xlab("")+
  theme_hc() +
  theme(legend.position = "none",
        axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold", size = 14))-> g2.2
g2.2

cov %>% 
  group_by(SurveyYear) %>% 
  select(SurveyYear, ForbPrp) %>% 
  pivot_longer(cols = ForbPrp,
               names_to = "covariate",
               values_to = "proportion") %>%  
  mutate(covariate= factor(covariate)) %>% 
  mutate(SurveyYear = factor(SurveyYear)) %>% 
  ggplot(aes(x=SurveyYear, 
             y = proportion,
             color = SurveyYear,
             fill = SurveyYear)) + 
  ggdist::stat_halfeye(
    height = -0.5, 
    width = .3, 
    justification = -.7, 
    point_colour = NA) + 
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_point(
    size = 1.3,
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  ylab("Proportion of forbs")+
  scale_fill_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  scale_color_manual(values = c("#fb8b24", "#5f0f40", "#14213d"))+
  xlab("")+
  ylim(0,1)+
  theme_hc() +
  theme(legend.position = "none",
        axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold", size = 14))-> g2.3
g2.3
library(cowplot)
plot_grid(g1, g2.1, g2.2, g2.3, ncol = 4)

ggsave(plot = last_plot(),
       path = "figures_05",
       filename = "covariates_dist.png",
       dpi = 300, 
       width = 20,
       height = 10, 
       units = "in")
#####################################

landcor <- landscape %>% 
  select(RowcropPrp_1k, PrairieSum_1k, RowcropTE)

# Correlation plot
library(corrplot)
cor <- cor(landcor)
corrplot(cor, type = 'lower', method = 'color', addCoef.col = 'black', number.cex = 0.7)

covcor <- cov %>% 
  select(VegHeight, GrassPrp, ForbPrp)
cor <- cor(covcor, use = "complete.obs")
corrplot(cor, type = 'lower', method = 'color', addCoef.col = 'black', number.cex = 0.7)

cor(landcor, covcor)


png(height=1800, width=1800, file="figs/correlation.png", type = "cairo")
corrplot(cor, type = 'lower', method = 'color', addCoef.col = 'black', tl.cex = 2.5, tl.col = 'black', number.cex = 2.5, col = COL2('PuOr', 10))
dev.off()

