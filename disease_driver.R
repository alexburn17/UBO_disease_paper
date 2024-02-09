# UBO Joint Disease Project
# P. Alexander Burnham
# 31 January 2023
# Last Modified: 13 November 2023

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
library(wesanderson)

# set directory:
setwd("~/Documents/GitHub/UBO_disease_paper")

# read in data
#ds <- read.csv("data/Chalkbrood_Aus.csv", header = TRUE, stringsAsFactors = FALSE)
ds <- read.csv("data/Chalkbrood_Aus.csv", header = TRUE, stringsAsFactors = FALSE)

ds23 <- read.csv("data/UBO_Data_2023.csv", header = TRUE, stringsAsFactors = FALSE)

# 2022 data
cleanDS_22 <- ds23[ds23$year == 2022,]




#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Chalk Brood
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################


# select relavant variables
ds <- select(ds, Test_No, Total_Brood, Brood_with_Chalk, White_Chalk, Spore_Dark_Chalk, Empty_Cells, Total_Cells, Cleared_Cells)

ds$Percentage_UBO <- ((ds$Cleared_Cells - ds$Empty_Cells) / ds$Total_Cells) * 100

# calculate percent chalkbrood
ds$percent_chalk <- ds$Brood_with_Chalk/ds$Total_Brood

# ubo binary
ds$UBO_binary <- ifelse(ds$Percentage_UBO < 20, "low UBO", "high UBO")

# make long form
ds_long <- gather(ds, chalk_type, chalkbrood, White_Chalk:Spore_Dark_Chalk, factor_key=TRUE)

# rename chalk type
ds_long$chalk_type <- ifelse(ds_long$chalk_type == "White_Chalk", "White Chalk", "Black Chalk Spores")

# make chalkbrood binary
ds_long$chalk_binary <- ifelse(ds_long$chalkbrood > 0, 1, 0)

# short data
ds_short <- ds_long[ds_long$chalk_type == "White Chalk",]

# total chalk
ds$TotalChalk = ds$White_Chalk + ds$Spore_Dark_Chalk


# ubo and chalk binary
UBO_chalk_binary <- ds_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(chalk_type, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(chalk_binary, na.rm=T), # mean
    n = length(chalk_binary),
    a = sum(chalk_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 



# ubo and chalk CONTINUOUS
UBO_chalk_brood <- ds_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(chalk_type, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(chalkbrood, na.rm=T), # mean
    n = length(chalkbrood),
    
  ) 



# PREVALENCE OF CHALKBROOD BY UBO SCORE
####################################################################################################

# model for binary chalk by ubo score
mod <- glm(data = ds_long, chalk_binary ~ Percentage_UBO * chalk_type, binomial(link = "logit"))
Anova(mod)

# plot chalk brood binary
ggplot(ds_long, aes(x=Percentage_UBO, y=chalk_binary, color=as.character(chalk_type))) +
  #geom_point(size=0) + 
  geom_smooth(method="glm", se=FALSE, fullrange=TRUE, size = 1, method.args = list(family=binomial)) +
  geom_point(size=2.3) +
  ylab("Chalkbrood Binary") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) +
  scale_color_manual(values = c("darkturquoise", "tomato3"), name=" ")  # color pallets option = A-H




# PERCENT INFECTED FRAMES BY UBO SCORE
####################################################################################################

# plot chalk brood
ggplot(ds_short, aes(x=Percentage_UBO, y=percent_chalk)) +
  #geom_point(size=0) + 
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE, size = 1) +
  geom_point(size=2.3) +
  ylab("% frames with chalkbrood") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) 
  

# model for percentage ubo on frames
mod <- lm(data = ds_short, percent_chalk ~ Percentage_UBO)
Anova(mod)

# original model
cor.test(x = ds_short$Percentage_UBO, y = ds_short$percent_chalk, method = 'spearman', alternative = "greater")






# CHALKBROOD CELLS BY UBO SCORE
####################################################################################################
# plot chalk brood
ggplot(ds_long, aes(x=Percentage_UBO, y= (chalkbrood + 1), 
                    color=as.character(chalk_type))) +
  #geom_point(size=0) + 
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE, size = 1) +
  geom_point(size=2.3) +
  ylab("Chalkbrood (cells/frame)") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3"), name=" ") + # color pallets option = A-H
  guides(color = guide_legend(override.aes = list(label = ''))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# glm for chalkbrood count data by ubo score and chalk type (poisson)
mod <- glm(data = ds_long, chalkbrood ~ Percentage_UBO * chalk_type, family = poisson(link = "log"))
Anova(mod)


mod <- lm(data = ds_long, chalkbrood ~ Percentage_UBO * chalk_type)
summary(mod)

# original model
cor.test(x = ds$Percentage_UBO, y = ds$TotalChalk, method = 'spearman', alternative = "greater")



# CHALKBROOD type regression
####################################################################################################
plot(y = log10(ds$White_Chalk+1), x = log10(ds$Spore_Dark_Chalk + 1))

x <- lm(data = ds, log10(Spore_Dark_Chalk + 1) ~ log10(White_Chalk+1))
summary(x)




#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2022 Nosema
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################

cleanDS_22 <- dplyr::select(cleanDS_22, year, beekeeper, sampling_event, lab_ID, assay_score, sampling_date, varroa_count, 
                            nosema_count, FKB_percentile)


# select ubo testing periods
x <- cleanDS_22[cleanDS_22$assay_score >= 0,]
uboMap <- select(x, lab_ID, assay_score)[complete.cases(select(x, lab_ID, assay_score)),]
cleanDS_22$assay_score <- NULL # delete assay score
cleanDS_22 <- merge(x = cleanDS_22, y = uboMap, by = "lab_ID", all.x=T) # merge UBO scores

# add month 
cleanDS_22$monthNum <- substr(cleanDS_22$sampling_date, 1, 1)
cleanDS_22$Month <- ifelse(cleanDS_22$monthNum == 5, "May", 
                           ifelse(cleanDS_22$monthNum == 6, "June", 
                                  ifelse(cleanDS_22$monthNum == 7, "July", 
                                         ifelse(cleanDS_22$monthNum == 8, "August",
                                                ifelse(cleanDS_22$monthNum == 9, "Sept.", NA)))))

# add binary variable
cleanDS_22$UBO_binary <- ifelse(cleanDS_22$assay_score >= 0.6, 1, 0) #"hygienic", "nonhygienic")
cleanDS_22$Nosema_binary <- ifelse(cleanDS_22$nosema_count > 0, 1, 0) 
cleanDS_22$Varroa_binary <- ifelse(cleanDS_22$varroa_count > 0, 1, 0)


# create anonymous beekeeper names
cleanDS_22$anonBeek <- ifelse(cleanDS_22$beekeeper == "Andrew Munkres", "beekeeper 1",
                              ifelse(cleanDS_22$beekeeper == "Mike Palmer", "beekeeper 2", 
                                     ifelse(cleanDS_22$beekeeper == "Jack Rath", "beekeeper 3", NA
                                     )))


##############################################################################
# Nosema 22
##############################################################################


# nosema summary
nosePrevSum <- cleanDS_22 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(Nosema_binary, na.rm=T), # mean
    n = length(Nosema_binary),
    a = sum(Nosema_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 

nosePrevSum <- nosePrevSum[complete.cases(nosePrevSum),]


nosemaLoad_Sum <- cleanDS_22 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(((nosema_count*4000000)/80), na.rm=T), # mean\
    n = length(((nosema_count*4000000)/80)),
    sd = sd(((nosema_count*4000000)/80), na.rm = TRUE),
    se = sd / sqrt(n)
    
  ) 

nosemaLoad_Sum <- nosemaLoad_Sum[complete.cases(nosemaLoad_Sum),]


# add factor data and make ubo a char
nosePrevSum <- nosePrevSum[!is.na(nosePrevSum$UBO_binary),]
nosePrevSum$Month <- factor(nosePrevSum$Month, levels = c("June", "July", "August", "Sept."))
nosePrevSum$UBO_Char <- ifelse(nosePrevSum$UBO_binary==1, "UBO High", "UBO Low")

# plot prevalence
nosPrev <- ggplot(nosePrevSum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.7,.95)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Nosema Prevalence", color=" ") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  scale_y_continuous(labels = scales::percent)



# add factor data and make ubo a char
nosemaLoad_Sum <- nosemaLoad_Sum[!is.na(nosemaLoad_Sum$UBO_binary),]
nosemaLoad_Sum$Month <- factor(nosemaLoad_Sum$Month, levels = c("June", "July", "August", "Sept."))
nosemaLoad_Sum$UBO_Char <- ifelse(nosemaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

contNos <-ggplot(nosemaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5) +
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(3,9)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Date", y="Nosema Load (spores/bee)", color=" ") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c( "tomato3", "darkturquoise")) 


# make a multi panel plot
plot_grid(nosPrev, contNos,
          labels = "AUTO", 
          label_size = 20)

mod2 <- glmer(data = cleanDS_22, Nosema_binary ~ UBO_binary * Month + (1 | lab_ID), family = binomial(link="logit"))
Anova(mod2)

mod3 <- glmer(data = cleanDS_22, (1+nosema_count) ~ UBO_binary * Month + (1 | lab_ID), family = Gamma(link=identity))
Anova(mod3)




#####################################################################################
############################ VIRUS DATA 2021 ########################################
#####################################################################################

# READ DATA
virus <- read.csv("data/UBO_VirusData_2021.csv", header = TRUE, stringsAsFactors = FALSE)
virus$RPS5..Cq. <- NULL # remove un-needed col

# make long form
virus_long <- gather(virus, virus, virus_load, DWV.ACopies.µl:IAPVCopies.µl, factor_key=TRUE)

# remove "Copies.ul" from virus name
virus_long$virus <- gsub('Copies.µl', '', virus_long$virus)

# ubo scores and binary virus and mites and log load
virus_long$virus_binary <- ifelse(virus_long$virus_load > 0, 1, 0)
virus_long$ubo_binary_june <- ifelse(virus_long$June.UBO > 60, "UBO +", "UBO -")
virus_long$ubo_binary_august <- ifelse(virus_long$August.UBO > 60, "UBO +", "UBO -")
virus_long$mite_binary_june <-  ifelse(virus_long$June.Mite > 0, 1, 0) 
virus_long$mite_binary_august <- ifelse(virus_long$August.Mite > 0, 1, 0)
virus_long$logLoad <- log10(virus_long$virus_load + 1)

virus_long <- virus_long[!is.na(virus_long$ubo_binary_august),]

# grouped boxplot
ggplot(data = virus_long, aes(x = virus, y = virus_load, fill = ubo_binary_august)) + 
  geom_boxplot() +
  theme_minimal(base_size = 20) +
  theme(legend.position = "top") +
  labs(x="Virus", y="Virus Load (copies/ul)", fill = "UBO Status:") +
  scale_fill_manual(values = c("#5071A0", "#E77624")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# virus binary summary
virusPrevSum <- virus_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(virus, ubo_binary_august) %>% # pick variables to group by
  dplyr::summarise(
    mean = mean(virus_binary, na.rm=T), # mean
    n = length(virus_binary),
    a = sum(virus_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
  ) 

virusPrevSum <- virusPrevSum[complete.cases(virusPrevSum),]


ggplot(virusPrevSum, aes(x=virus, y=mean, fill=ubo_binary_august)) + 
  theme_minimal(base_size = 20) +
  theme(legend.position = "top") +
  labs(x="Virus", y="Virus Prevalence", fill = "UBO Status:") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("#5071A0", "#E77624")) +
  scale_y_continuous(labels = scales::percent)

# make data wide
virus_wide <- select(virus_long, -virus_binary, -ubo_binary_june, -ubo_binary_august,-mite_binary_june, -virus_load)
virus_wide <- virus_wide %>% pivot_wider(names_from = virus, values_from = logLoad)
virus_wide$ubo_binary_june <- ifelse(virus_wide$June.UBO > 60, "UBO +", "UBO -")
virus_wide$ubo_binary_august <- ifelse(virus_wide$August.UBO > 60, "UBO +", "UBO -")

res.man <- manova(cbind(June.Mite,August.Mite,DWV.A,DWV.B,LSV,SBV,BQCV,IAPV) ~ ubo_binary_august, data = virus_wide)
summary.aov(res.man)







