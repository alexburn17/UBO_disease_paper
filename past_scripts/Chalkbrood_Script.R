# UBO Joint Project
# P. Alexander Burnham
# 31 January 2023


# set directory:
setwd("~/Documents/GitHub/UBO_disease_paper")


# install libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(tidyr)
library(viridis)
library(car)
library(imputeTS)
library(cowplot)
library(scales)
library(ggmosaic)


# read in chalkbrood data
ds <- read.csv("data/Chalkbrood_Aus.csv", header = TRUE, stringsAsFactors = FALSE)

# select relavant variables
ds <- select(ds, Test_No, Total_Brood, Brood_with_Chalk, White_Chalk, Spore_Dark_Chalk, Percentage_UBO)

# calculate percent chalkbrood
ds$percent_chalk <- ds$Brood_with_Chalk/ds$Total_Brood

# calculate chalk binary
ds$chalk_binary <- ifelse(ds$percent_chalk > 0, 1, 0)

# ubo binary
ds$UBO_binary <- ifelse(ds$Percentage_UBO < 20, "low UBO", "high UBO")

# make long form
ds_long <- gather(ds, chalk_type, chalkbrood, White_Chalk:Spore_Dark_Chalk, factor_key=TRUE)

# rename chalk type
ds_long$chalk_type <- ifelse(ds_long$chalk_type == "White_Chalk", "White Chalk", "Black Chalk Spores")

# plot chalk brood
ggplot(ds_long, aes(x=Percentage_UBO, y=chalkbrood, 
               color=as.character(chalk_type))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1) +
  geom_point(size=2.3) +
  ylab("Chalkbrood (cells/frame)") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3"), name=" ") + # color pallets option = A-H
  guides(color = guide_legend(override.aes = list(label = '')))


# plot chalk brood
ggplot(ds_long, aes(x=Percentage_UBO, y=percent_chalk)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1) +
  geom_point(size=2.3) +
  ylab("% frames with chalkbrood") + # y axis label
  xlab("Proportion UBO Response") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) 

ds

varroaPrevSum <- ds %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(varroa_binary, na.rm=T), # mean
    n = length(varroa_binary),
    a = sum(varroa_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 






mod <- aov(data = ds_long, chalkbrood ~ UBO_binary * chalk_type)
summary(mod)


mod <- lm(data = ds_long, chalkbrood ~ Percentage_UBO * chalk_type)
Anova(mod)

cor.test(x = ds_long$Percentage_UBO, y = ds_long$chalkbrood, method = 'spearman', alternative = "greater")


mod <- lm(data = ds, Total_Brood ~ Percentage_UBO)
Anova(mod)


cor.test(x = ds_long$Percentage_UBO, y = ds_long$chalkbrood, method = 'spearman', alternative = "greater")







