# UBO Bee Breeder Project Data Analysis
# P. Alexander Burnham, S. Miller, C. McKay, and Samantha Alger
# 30 November 2022


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
library(wesanderson)



# load in data
#ds <- read.csv("SARE_Field_database2022.csv", header = TRUE, stringsAsFactors = FALSE)
ds22 <- read.csv("data/UBO_Data_2022.csv", header = TRUE, stringsAsFactors = FALSE)
ds23 <- read.csv("data/UBO_Data_2023.csv", header = TRUE, stringsAsFactors = FALSE)

# 2022 data
cleanDS_22 <- ds23[ds23$year == 2022,]

#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2023 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################

# simplify dataset
cleanDS <- dplyr::select(ds23, year, beekeeper, sampling_event, lab_ID, assay_score, sampling_date, varroa_count, 
                   nosema_count, FKB_percentile, project, breeder_colony_UBO_status, breeder_mating_system, queen_mating_system)

cleanDS <- cleanDS[cleanDS$year == 2023,]

# create heritibility dataset
herit <- cleanDS[cleanDS$project=="heritability",]

# create main 2023 dataset
cleanDS <- cleanDS[!cleanDS$project=="heritability",]

# select ubo testing periods
x <- cleanDS[cleanDS$assay_score >= 0,]
uboMap <- select(x, lab_ID, assay_score)[complete.cases(select(x, lab_ID, assay_score)),]
cleanDS$assay_score <- NULL # delete assay score
cleanDS <- merge(x = cleanDS, y = uboMap, by = "lab_ID", all.x=T) # merge UBO scores

# add month 
cleanDS$monthNum <- substr(cleanDS$sampling_date, 1, 1)
cleanDS$Month <- ifelse(cleanDS$monthNum == 5, "May", 
                       ifelse(cleanDS$monthNum == 6, "June", 
                              ifelse(cleanDS$monthNum == 7, "July", 
                                     ifelse(cleanDS$monthNum == 8, "August", NA))))

# add binary variable
cleanDS$UBO_binary <- ifelse(cleanDS$assay_score >= 0.6, 1, 0) #"hygienic", "nonhygienic")
cleanDS$Nosema_binary <- ifelse(cleanDS$nosema_count > 0, 1, 0) 


##############################################################################
# Heritabilty 23
##############################################################################

herit$breeder <- paste0(herit$breeder_mating_system, " + " ,herit$queen_mating_system)

herit$breeder <- ifelse(herit$breeder_colony_UBO_status == "low", "OM + OM (low)", herit$breeder)

herit <- herit[!herit$breeder=="OM + II",]

herit$UBO_binary <- ifelse(herit$assay_score >= 0.6, 1, 0)
thresh <- mean(ifelse(cleanDS_22$assay_score >= .6, 1, 0), na.rm=T)

# ubo summary
uboPrevSum <- herit %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(breeder) %>% # pick variables to group by
  summarise(
    
    mean = mean(UBO_binary, na.rm=T), # mean
    n = length(UBO_binary),
    a = sum(UBO_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 




# Hygienic behavior by breeding type
ggplot(herit, aes(x=breeder, y=assay_score)) + 
  geom_boxplot(size=1, outlier.shape = NA, aes(color=breeder)) +
  geom_jitter(size=3, aes(color=breeder)) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("UBO Score") + # y axis label
  xlab("Breeding System") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "orange", "firebrick", "forestgreen")) +
  geom_hline(yintercept=c(.6), linetype="dashed",
             color = c("black"), size=1) +
  stat_summary(fun=mean, geom="point", shape=18, size=8, color="black") +
  scale_y_continuous(labels = scales::percent)

# percent UBO by breeding type
ggplot(uboPrevSum, aes(x=breeder, y=mean, color = breeder)) + 
  geom_point(size = 10, shape = 18) +
  theme_classic(base_size = 17) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, .77)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.2, color=breeder))+
  labs(x="Breeding System", y="Percent High UBO", color="UBO Status:") +
  scale_color_manual(values = c("darkturquoise", "orange", "firebrick", "forestgreen")) +
  geom_hline(yintercept=thresh, linetype="dashed",
             color = c("black"), size=1) +
  scale_y_continuous(labels = scales::percent)


mod4 <- aov(data=herit, assay_score~breeder)
summary(mod4)

mod5 <- glm(data = herit, UBO_binary~breeder, family = binomial(link="logit"))
Anova(mod5)






##############################################################################
# Nosema 23
##############################################################################


# nosema summary
nosePrevSum <- cleanDS %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
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


nosemaLoad_Sum <- cleanDS %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
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
nosePrevSum$Month <- factor(nosePrevSum$Month, levels = c("May", "June", "July"))
nosePrevSum$UBO_Char <- ifelse(nosePrevSum$UBO_binary==1, "UBO High", "UBO Low")

# plot prevalence
nosPrev <- ggplot(nosePrevSum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Nosema Prevalence", color="UBO Status:") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  scale_y_continuous(labels = scales::percent)



# add factor data and make ubo a char
nosemaLoad_Sum <- nosemaLoad_Sum[!is.na(nosemaLoad_Sum$UBO_binary),]
nosemaLoad_Sum$Month <- factor(nosemaLoad_Sum$Month, levels = c("May", "June", "July"))
nosemaLoad_Sum$UBO_Char <- ifelse(nosemaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

contNos <-ggplot(nosemaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.6,.9)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Date", y="Nosema Load (spores/bee)", color=" ") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c( "tomato3", "darkturquoise")) +
  geom_hline(yintercept=1000000, linetype="dashed",
             color = c("black"), size=1)


# make a multi panel plot
plot_grid(nosPrev, contNos,
          labels = "AUTO", 
          label_size = 20)

mod2 <- glmer(data = cleanDS, Nosema_binary ~ UBO_binary * Month + (1 | lab_ID), family = binomial(link="logit"))
Anova(mod2)

mod3 <- glmer(data = cleanDS, nosema_count ~ UBO_binary * Month + (1 | lab_ID), family = "gaussian")
Anova(mod3)



ggplot(cleanDS_no0, aes(x = assay_score, y = nosema_count, color = Month)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 17) +
  ylab("Nosema Load") + # y axis label
  xlab("UBO") +# x axis label
  guides(color=guide_legend(title="Month:")) +
  geom_smooth(method='lm', formula= y~x, se=F, size = 2) +
  theme(legend.position="top") +
  scale_color_manual(values = c("tomato3", "darkturquoise", "goldenrod"))




##############################################################################
# Varroa 23
##############################################################################


# remove 0s
cleanDS$varroa_binary <- ifelse(cleanDS$varroa_count > 0, 1, 0)
Varroa_no0 <-cleanDS[!cleanDS$varroa_binary==0,]
Varroa_no0$varroa_load <- as.numeric(Varroa_no0$varroa_count)

varroaPrevSum <- cleanDS %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(varroa_binary, na.rm=T), # mean
    n = length(varroa_binary),
    a = sum(varroa_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 
varroaPrevSum <- varroaPrevSum[complete.cases(varroaPrevSum),]



# add factor data and make ubo a char
varroaPrevSum <- varroaPrevSum[!is.na(varroaPrevSum$UBO_binary),]
varroaPrevSum$Month <- factor(varroaPrevSum$Month, levels = c("May", "June", "July"))
varroaPrevSum$UBO_Char <- ifelse(varroaPrevSum$UBO_binary==1, "UBO High", "UBO low")


# plot prevalence
varPrev <- ggplot(varroaPrevSum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(4,5)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Prevalence", color="") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  scale_y_continuous(labels = scales::percent)



varroaLoad_Sum <- cleanDS %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(100*(varroa_count/300), na.rm=T), # mean\
    n = length(100*(varroa_count/300)),
    sd = sd(100*(varroa_count/300), na.rm = TRUE),
    se = sd / sqrt(n)
    
  ) 
varroaLoad_Sum <- varroaLoad_Sum[complete.cases(varroaLoad_Sum),]

# add factor data and make ubo a char
varroaLoad_Sum <- varroaLoad_Sum[!is.na(varroaLoad_Sum$UBO_binary),]
varroaLoad_Sum$Month <- factor(varroaLoad_Sum$Month, levels = c("May", "June", "July"))
varroaLoad_Sum$UBO_Char <- ifelse(varroaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

varLoad <-ggplot(varroaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.4,.9)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Load (mites/100 bees)", color=" ") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  geom_hline(yintercept=c(2), linetype="dashed",
             color = c("black"), size=1)
varLoad


# make a multi panel plot
plot_grid(varPrev, varLoad,
          labels = "AUTO", 
          label_size = 20)


mod1 <- glmer(data = cleanDS, (varroa_count+1) ~ UBO_binary * Month + (1 | lab_ID), family = "Gamma")
Anova(mod1)

mod2 <- glmer(data = cleanDS, varroa_binary ~ UBO_binary * Month + (1 | lab_ID), family = binomial(link="logit"))
Anova(mod2)


ggplot(Varroa_no0, aes(x = assay_score, y = varroa_load, color = Month)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 17) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("UBO") +# x axis label
  guides(color=guide_legend(title="Month:")) +
  geom_smooth(method='lm', formula= y~x, se=F, size = 2) +
  theme(legend.position="top") +
  scale_color_manual(values = c("tomato3", "darkturquoise", "goldenrod"))

mod <- glmer(data = Varroa_no0[Varroa_no0$Month=="May",], varroa_load ~ assay_score + (1 | beekeeper), family = "Gamma")
Anova(mod)


#############################################################################
# Beekeeper 23
#############################################################################

colors <- wes_palette("Zissou1", 5, "discrete")[c(1,2,4,5)]

cleanDS$anonBeek <- ifelse(cleanDS$beekeeper=="Andrew Munkres", "Beekeeper 1", 
                          ifelse(cleanDS$beekeeper=="Mike Palmer", "Beekeeper 2",
                                 ifelse(cleanDS$beekeeper=="Jack Rath", "Beekeeper 3", "Beekeeper 4")))


# ubo summary
uboPrevSumBeek <- cleanDS %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(anonBeek) %>% # pick variables to group by
  summarise(
    
    mean = mean(UBO_binary, na.rm=T), # mean
    n = length(UBO_binary),
    a = sum(UBO_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 


# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(cleanDS, aes(x=anonBeek, y=assay_score, color=beekeeper)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  geom_jitter(size=3) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("UBO Score") + # y axis label
  xlab("Beekeeper") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_manual(values = c("#3B9AB2", "#78B7C5", "#F21A00","#E1AF00")) +
  geom_hline(yintercept=.6, linetype="dashed",
             color = "black", size=1) +
  scale_y_continuous(labels = scales::percent) +
  stat_summary(fun=mean, geom="point", shape=18, size=8, color="black")



# percent UBO by breeding type
ggplot(uboPrevSumBeek, aes(x=anonBeek, y=mean, color = anonBeek)) + 
  geom_point(size = 10, shape = 18) +
  theme_classic(base_size = 17) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, .77)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.2, color=anonBeek))+
  labs(x="Beekeeper", y="Percent High UBO", color="UBO Status:") +
  scale_color_manual(values = c("#3B9AB2", "#E1AF00", "#F21A00", "#78B7C5")) +
  geom_hline(yintercept=thresh, linetype="dashed",
             color = c("black"), size=1) +
  scale_y_continuous(labels = scales::percent)


#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2022
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


##############################################################################
# Varroa 22
##############################################################################


varroaPrevSum <- cleanDS_22 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(Varroa_binary, na.rm=T), # mean
    n = length(Varroa_binary),
    a = sum(Varroa_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 
varroaPrevSum <- varroaPrevSum[complete.cases(varroaPrevSum),]



# add factor data and make ubo a char
varroaPrevSum <- varroaPrevSum[!is.na(varroaPrevSum$UBO_binary),]
varroaPrevSum$Month <- factor(varroaPrevSum$Month, levels = c("June", "July", "August", "Sept."))
varroaPrevSum$UBO_Char <- ifelse(varroaPrevSum$UBO_binary==1, "UBO High", "UBO low")


# plot prevalence
varPrev <- ggplot(varroaPrevSum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(4,5)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Prevalence", color="") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  scale_y_continuous(labels = scales::percent)



varroaLoad_Sum <- cleanDS_22 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(Month, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(100*(varroa_count/300), na.rm=T), # mean\
    n = length(100*(varroa_count/300)),
    sd = sd(100*(varroa_count/300), na.rm = TRUE),
    se = sd / sqrt(n)
    
  ) 
varroaLoad_Sum <- varroaLoad_Sum[complete.cases(varroaLoad_Sum),]

# add factor data and make ubo a char
varroaLoad_Sum <- varroaLoad_Sum[!is.na(varroaLoad_Sum$UBO_binary),]
varroaLoad_Sum$Month <- factor(varroaLoad_Sum$Month, levels = c("June", "July", "August", "Sept."))
varroaLoad_Sum$UBO_Char <- ifelse(varroaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

varLoad <-ggplot(varroaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.3,.82)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Load (mites/100 bees)", color=" ") +
  scale_color_manual(values = c("tomato3", "darkturquoise")) +
  geom_hline(yintercept=c(2), linetype="dashed",
             color = c("black"), size=1)
varLoad


# make a multi panel plot
plot_grid(varPrev, varLoad,
          labels = "AUTO", 
          label_size = 20)


mod1 <- glmer(data = cleanDS_22, (varroa_count+1) ~ UBO_binary * Month + (1 | lab_ID), family = Gamma(link="identity"))
Anova(mod1)

mod2 <- glmer(data = cleanDS_22, Varroa_binary ~ UBO_binary * Month + (1 | lab_ID), family = binomial(link="logit"))
Anova(mod2)


#############################################################################
# Beekeeper 22
#############################################################################


cleanDS_22$anonBeek <- ifelse(cleanDS_22$beekeeper=="Andrew Munkres", "Beekeeper 1", 
                           ifelse(cleanDS_22$beekeeper=="Mike Palmer", "Beekeeper 2",
                                  ifelse(cleanDS_22$beekeeper=="Jack Rath", "Beekeeper 3", "Beekeeper 4")))


# ubo summary
uboPrevSumBeek <- cleanDS_22 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(anonBeek) %>% # pick variables to group by
  summarise(
    
    mean = mean(UBO_binary, na.rm=T), # mean
    n = length(UBO_binary),
    a = sum(UBO_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 


# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(cleanDS_22, aes(x=anonBeek, y=assay_score, color=beekeeper)) + 
  geom_boxplot(size=1, outlier.shape = NA) +
  geom_jitter(size=3) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("UBO Score") + # y axis label
  xlab("Beekeeper") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_manual(values = c("#3B9AB2", "#F21A00","#E1AF00")) +
  geom_hline(yintercept=.6, linetype="dashed",
             color = "black", size=1) +
  scale_y_continuous(labels = scales::percent) +
  stat_summary(fun=mean, geom="point", shape=18, size=8, color="black")



# percent UBO by breeding type
ggplot(uboPrevSumBeek, aes(x=anonBeek, y=mean, color = anonBeek)) + 
  geom_point(size = 10, shape = 18) +
  theme_classic(base_size = 17) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, .77)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.2, color=anonBeek))+
  labs(x="Beekeeper", y="Percent High UBO", color="UBO Status:") +
  scale_color_manual(values = c("#3B9AB2", "#E1AF00", "#F21A00")) +
  geom_hline(yintercept=thresh, linetype="dashed",
             color = c("black"), size=1) +
  scale_y_continuous(labels = scales::percent)



#############################################################################
# Survival 22
#############################################################################

surv23 <- cleanDS[cleanDS$beekeeper=="Mike Palmer",]
surv22 <- cleanDS_22[cleanDS_22$beekeeper=="Mike Palmer",]


# add survival data to 22
surv22$survival <- surv22$lab_ID %in% surv23$lab_ID
surv22_clean <- distinct(select(surv22, UBO_binary, assay_score, survival, lab_ID))

surv22_clean$surv_BINARY <- ifelse(surv22_clean$survival==TRUE, 1, 0)

x = aov(surv22_clean$assay_score~surv22_clean$surv_BINARY)
summary(x)

mod <- glm(data = surv22_clean, surv_BINARY ~ UBO_binary, family = binomial(link="logit"))
Anova(mod)





































# OLD
#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2022
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################


# UBO cont and binary 
# create binary variable for UBO 
ds22$UBO_binary <- ifelse(ds22$assay_score >= 0.6, 1, 0) #"hygienic", "nonhygienic")
mean(ds22$UBO_binary, na.rm=T) # get percentage of hygienic UBO

# create anonymous beekeeper names
ds22$anonBeek <- ifelse(ds$beekeeper == "Andrew Munkres", "beekeeper 1",
                        ifelse(ds$beekeeper == "Jack Rath", "beekeeper 2", "beekeeper 3"
                        ))



#################################################################################
# NOSEMA Analysis
#################################################################################

# create nosema data frame and make long form
NosemaDS <- select(ds22, beekeeper, yard, lab_ID, june_nosema_load_spores.bee, august_nosema_load_spores.bee, UBO_binary, assay_score)
NosemaDS_long <- gather(NosemaDS, time, nosmea_load, june_nosema_load_spores.bee:august_nosema_load_spores.bee, factor_key=TRUE)
NosemaDS_long$time <- ifelse(NosemaDS_long$time=="june_nosema_load_spores.bee", "June", "August")
NosemaDS_long$nosmea_load_log <- log10(NosemaDS_long$nosmea_load + 1)
NosemaDS_long$nosema_binary <- ifelse(NosemaDS_long$nosmea_load > 0, 1, 0)
NosemaDS_long$rescaledNosema <- NosemaDS_long$nosmea_load/sum(NosemaDS_long$nosmea_load, na.rm = TRUE)
NosemaDS_long$lab_ID <- as.character(NosemaDS_long$lab_ID)

# remove 0s
NosemaDS_long_no0 <- NosemaDS_long[!NosemaDS_long$nosema_binary==0,]

x=NosemaDS_long_no0[NosemaDS_long_no0$time=="August",]
x[x$UBO_binary==1,]

nosePrevSum <- NosemaDS_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(nosema_binary, na.rm=T), # mean
    n = length(nosema_binary),
    a = sum(nosema_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 

# add factor data and make ubo a char
nosePrevSum <- nosePrevSum[!is.na(nosePrevSum$UBO_binary),]
nosePrevSum$time <- factor(nosePrevSum$time, levels = c("June", "August"))
nosePrevSum$UBO_Char <- ifelse(nosePrevSum$UBO_binary==1, "UBO Pos.", "UBO Neg.")

# plot prevalence
nosPrev <- ggplot(nosePrevSum, aes(x=time, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Nosema Prevalence", color="UBO Status:") +
  scale_color_manual(values = c("tomato3", "darkturquoise"))


nosemaLoad_Sum <- NosemaDS_long_no0 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(nosmea_load, na.rm=T), # mean\
    n = length(nosmea_load),
    sd = sd(nosmea_load, na.rm = TRUE),
    se = sd / sqrt(n)

) 

# add factor data and make ubo a char
nosemaLoad_Sum <- nosemaLoad_Sum[!is.na(nosemaLoad_Sum$UBO_binary),]
nosemaLoad_Sum$time <- factor(nosemaLoad_Sum$time, levels = c("June", "August"))
nosemaLoad_Sum$UBO_Char <- ifelse(nosemaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

contNos <-ggplot(nosemaLoad_Sum, aes(x=time, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.2,.9)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Date", y="Nosema Load (spores/bee)", color=" ") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("tomato3", "darkturquoise"))
contNos


# make a multi panel plot
plot_grid(nosPrev, contNos,
          labels = "AUTO", 
          label_size = 20)



#glm with gamma distribution rescaled nosema assay score by time
mod <- glmer(data = NosemaDS_long_no0, rescaledNosema ~ assay_score * time + (1 | yard/beekeeper), family = "Gamma")
mod1 <- glmer(data = NosemaDS_long, nosema_binary ~ assay_score * time + (1 | yard/beekeeper), family = binomial(link="logit"))
Anova(mod)
Anova(mod1)

#glm with gamma distribution rescaled nosema ubo binary by time
mod2 <- glmer(data = NosemaDS_long, nosema_binary ~ UBO_binary * time + (1 | yard/beekeeper), family = binomial(link="logit"))
mod3 <- glmer(data = NosemaDS_long_no0, rescaledNosema ~ UBO_binary * time + (1 | yard/beekeeper), family = "Gamma")
Anova(mod2)
Anova(mod3)






#################################################################################
# VARROA Analysis
#################################################################################


VarroaDS <- select(ds, beekeeper, yard, lab_ID, june_varroa_load_mites.100.bees, august_varroa_load_mites.100.bees, UBO_binary, assay_score)
VarroaDS_long <- gather(VarroaDS, time, varroa_load, june_varroa_load_mites.100.bees:august_varroa_load_mites.100.bees, factor_key=TRUE)
VarroaDS_long$time <- ifelse(VarroaDS_long$time=="june_varroa_load_mites.100.bees", "June", "August")
VarroaDS_long$varroa_binary <- ifelse(VarroaDS_long$varroa_load > 0, 1, 0)
VarroaDS_long$rescaledVarroa <- VarroaDS_long$varroa_load/sum(VarroaDS_long$varroa_load, na.rm = TRUE)
VarroaDS_long$lab_ID <- as.character(VarroaDS_long$lab_ID)


# remove 0s
VarroaDS_long_no0 <- VarroaDS_long[!VarroaDS_long$varroa_binary==0,]


varroaPrevSum <- VarroaDS_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(varroa_binary, na.rm=T), # mean
    n = length(varroa_binary),
    a = sum(varroa_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
    
  ) 


# add factor data and make ubo a char
varroaPrevSum <- varroaPrevSum[!is.na(varroaPrevSum$UBO_binary),]
varroaPrevSum$time <- factor(varroaPrevSum$time, levels = c("June", "August"))
varroaPrevSum$UBO_Char <- ifelse(varroaPrevSum$UBO_binary==1, "UBO Pos.", "UBO Neg.")


# plot prevalence
nosPrev <- ggplot(varroaPrevSum, aes(x=time, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(8,8)) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Prevalence", color="UBO Status:") +
  scale_color_manual(values = c("tomato3", "darkturquoise"))
nosPrev





varroaLoad_Sum <- VarroaDS_long_no0 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, UBO_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(varroa_load, na.rm=T), # mean\
    n = length(varroa_load),
    sd = sd(varroa_load, na.rm = TRUE),
    se = sd / sqrt(n)
    
  ) 


# add factor data and make ubo a char
varroaLoad_Sum <- varroaLoad_Sum[!is.na(varroaLoad_Sum$UBO_binary),]
varroaLoad_Sum$time <- factor(varroaLoad_Sum$time, levels = c("June", "August"))
varroaLoad_Sum$UBO_Char <- ifelse(varroaLoad_Sum$UBO_binary==1, "UBO High", "UBO Low")

varLoad <-ggplot(varroaLoad_Sum, aes(x=time, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5)+
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.2,.9)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Month", y="Varroa Load", color=" ") +
  scale_color_manual(values = c("tomato3", "darkturquoise"))
varLoad


# make a multi panel plot
plot_grid(nosPrev, varLoad,
          labels = "AUTO", 
          label_size = 20)



#glm with gamma distribution rescaled nosema assay score by time
mod <- glmer(data = VarroaDS_long_no0, rescaledVarroa ~ assay_score * time + (1 | yard/beekeeper), family = "Gamma")
mod1 <- glmer(data = VarroaDS_long, varroa_binary ~ assay_score * time + (1 | yard/beekeeper), family = binomial(link="logit"))
Anova(mod)
Anova(mod1)

#glm with gamma distribution rescaled nosema ubo binary by time
mod2 <- glmer(data = VarroaDS_long, varroa_binary ~ UBO_binary * time + (1 | yard/beekeeper), family = binomial(link="logit"))
mod3 <- glmer(data = VarroaDS_long_no0, rescaledVarroa ~ UBO_binary * time + (1 | yard/beekeeper), family = "Gamma")
Anova(mod2)
Anova(mod3)


x=merge(NosemaDS_long, VarroaDS_long)
chisq.test(x$varroa_binary, x$nosema_binary)


x %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(time, varroa_binary) %>% # pick variables to group by
  summarise(
    
    mean = mean(nosema_binary, na.rm=T), # mean
    n = length(nosema_binary),
    sd = sd(nosema_binary, na.rm = TRUE),
    se = sd / sqrt(n)
    
  ) 

































# log transform data
ds22$log_june_varroa_load_mites.100.bees <- log10(ds22$june_varroa_load_mites.100.bees + 0.0001)
ds22$log_august_varroa_load_mites.100.bees <- log10(ds22$august_varroa_load_mites.100.bees + 0.0001)
ds22$log_june_nosema_load_spores.bee <- log10(ds22$june_nosema_load_spores.bee + 1)

# split the data by beekeeper
dsSplit <- split(ds22, ds22$beekeeper)

## UBO by June Varroa loads continuous 
# Add regression lines
juneBeek <- ggplot(ds22, aes(x=assay_score, y=june_varroa_load_mites.100.bees, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab(NULL) + # y axis label
  xlab(NULL) + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.8,.8)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
juneBeek

# discrete analysis
dsNo0 <- ds[!ds$june_varroa_load_mites.100.bees==0,]
boxplot(dsNo0$UBO_binary, dsNo0$june_varroa_load_mites.100.bees)
x = aov(dsNo0$june_varroa_load_mites.100.bees~dsNo0$UBO_binary)
summary(x)


# june varroa for just andrew by ubo
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$june_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$june_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$june_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO by August Varroa loads continuous
# Add regression lines
augBeek <- ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab(NULL) + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(5, 5)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
augBeek

# august varroa load by beekeeper
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$august_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$august_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$august_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO June composite continuous 
# Add regression lines
juneVar <- ggplot(ds22, aes(x=assay_score, y=june_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "darkturquoise") +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab("June Varroa Load (mites/100 bees)") + # y axis label
  xlab(NULL) + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text("aes(label=lab_ID)", ) +
  guides(color = guide_legend(override.aes = list(label = '')))
juneVar

# june varroa  by ubo
cor.test(ds$assay_score, ds$june_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO August composite continuous 
# Add regression lines
augVar <- ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "darkturquoise") +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab("August Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
augVar

# june varroa  by ubo
cor.test(ds$assay_score, ds$august_varroa_load_mites.100.bees, method="spearman", exact = F)



# make a multi panel plot
plot_grid(juneVar, juneBeek, augVar, augBeek,
          labels = "AUTO", 
          label_size = 17)





# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(ds, aes(x=anonBeek, y=assay_score, color=beekeeper)) + 
  geom_boxplot(size=1) +
  geom_point(size=3) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Beekeeper") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey")) +
  geom_hline(yintercept=.6, linetype="dashed",
             color = "black", size=1)



ds %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(beekeeper) %>% # pick variables to group by
  summarise(
    
    mean = mean(UBO_binary, na.rm=T), # mean
    sum = sum(UBO_binary, na.rm=T),
    N = length(UBO_binary), # sample size
    
  ) 



## UBO by June varroa load binary 
boxplot(ds$june_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$june_varroa_load_mites.100.bees~ds$UBO_binary))





## UBO by August Varroa loads binary
boxplot(ds$august_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$august_varroa_load_mites.100.bees~ds$UBO_binary))






#################################################################################
# Nosema Analysis
#################################################################################


## UBO by June Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=june_nosema_load_spores.bee, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  ylab("June Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$june_nosema_load_spores.bee, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$june_nosema_load_spores.bee, method="spearman", exact = F)
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$june_nosema_load_spores.bee, method="spearman", exact = F)



## UBO by June Nosema load continuous wit a log transform
# Add regression lines
ggplot(ds, aes(x=assay_score, y=log_june_nosema_load_spores.bee, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  ylab("June Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


summary(lm(dsSplit$`Andrew Munkres`$log_june_nosema_load_spores.bee ~ dsSplit$`Andrew Munkres`$assay_score))


# remove 0s
dsNos_no0 <- ds[!ds$june_nosema_load_spores.bee==0,] 

# remove NA
dsNos_no0 <- dsNos_no0[!is.na(dsNos_no0$june_nosema_load_spores.bee), ]
dsNos_no0 <- dsNos_no0[!is.na(dsNos_no0$assay_score), ]






dsNos_no0$rescaledNosema <- dsNos_no0$june_nosema_load_spores.bee/sum(dsNos_no0$june_nosema_load_spores.bee)










# june nosema by ubo
cor.test(ds$assay_score, ds$log_june_nosema_load_spores.bee, method="spearman", exact = F)


summary(lm(dsNos_no0$log_june_nosema_load_spores.bee ~ dsNos_no0$assay_score))
plot(dsNos_no0 $assay_score, dsNos_no0 $log_june_nosema_load_spores.bee)

## UBO by June Nosema load binary 
boxplot(dsNos_no0$log_june_nosema_load_spores.bee~dsNos_no0$UBO_binary)
summary(aov(ds$log_june_nosema_load_spores.bee~ds$UBO_binary))

t.test(ds$june_nosema_load_spores.bee~ds$UBO_binary, alternative="greater")

kruskal.test(ds$june_nosema_load_spores.bee~ds$UBO_binary)


 ## UBO by August Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=august_nosema_load_spores.bee, 
               color=as.character(beekeeper))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("August Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## UBO by June Nosema load binary 
boxplot(ds$august_nosema_load_spores.bee~ds$UBO_binary)
summary(aov(ds$august_nosema_load_spores.bee~ds$UBO_binary))



# here we can create a character variable below 10 is small greater than ten is large


# UBO by frames of brood 
ggplot(ds, aes(x=frames_of_brood, y=assay_score, 
               color=as.character(beekeeper))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("UBO Assay Score)") + # y axis label
  xlab("Colony Size") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))






