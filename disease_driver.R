# UBO Joint Disease Project
# P. Alexander Burnham
# 31 January 2023
# Last Modified: 1 May 2024

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
library(multcomp)
library(gridExtra)

# set directory:
setwd("~/Documents/GitHub/UBO_disease_paper")

# read in data
ds <- read.csv("data/Chalkbrood_Aus.csv", header = TRUE, stringsAsFactors = FALSE)
ds23 <- read.csv("data/UBO_Data_2023.csv", header = TRUE, stringsAsFactors = FALSE)

# 2022 data
cleanDS_22 <- ds23[ds23$year == 2022,]




#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Chalk Brood
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################


# select relevant variables
ds <- dplyr::select(ds, Test_No, Total_Brood, Brood_with_Chalk, White_Chalk, Spore_Dark_Chalk, Empty_Cells, Total_Cells, Cleared_Cells)

ds$Percentage_UBO <- ((ds$Cleared_Cells - ds$Empty_Cells) / ds$Total_Cells) * 100

# calculate percent chalkbrood
ds$percent_chalk <- ds$Brood_with_Chalk/ds$Total_Brood

# ubo binary
ds$UBO_binary <- ifelse(ds$Percentage_UBO < 20, "low UBeeO", "high UBeeO")

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
summary(mod)

# plot chalk brood binary
ggplot(ds_long, aes(x=Percentage_UBO, y=chalk_binary, color=as.character(chalk_type))) +
  geom_smooth(method="glm", se=FALSE, fullrange=TRUE, size = 2, method.args = list(family=binomial)) +
  geom_point(size=4) +
  ylab("Chalkbrood Binary") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 20) + # size of the text and label ticks
  theme(legend.position = "none") +
  scale_color_manual(values = c("#519E9A", "#9E519F"), name=" ")  # color pallets option = A-H




# PERCENT INFECTED FRAMES BY UBO SCORE
####################################################################################################

# plot chalk brood
ggplot(ds_short, aes(x=Percentage_UBO, y=percent_chalk)) +
  #geom_point(size=0) + 
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE, size = 2, color = c("#E77624")) +
  geom_point(size=2.8) +
  ylab("% frames with chalkbrood") + # y axis label
  xlab("Percent UBO Response") + # x axis label
  theme_minimal(base_size = 20) + # size of the text and label ticks
  theme(legend.position = c(.76, .9)) 
  

# model for percentage ubo on frames
mod <- lm(data = ds_short, percent_chalk ~ Percentage_UBO)
Anova(mod)
summary(mod)



# CHALKBROOD CELLS BY UBO SCORE
####################################################################################################
# plot chalk brood
ds_long$uboScore <- ds_long$Percentage_UBO/100
ds_long$logChalk <- log10(ds_long$chalkbrood+1)


bot <- ggplot(ds_long, aes(x=uboScore, y = (chalkbrood), 
                    linetype=as.character(chalk_type), shape = as.character(chalk_type))) +
  geom_point(size=6) +
  ylab(" ") + # y axis label
  xlab(" ") + # x axis label
  theme_bw(base_size = 20) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  coord_cartesian(xlim = c(0, .5), ylim = c(0, 50)) +
  scale_linetype_manual(values = c(1, 3), name=" ", guide = FALSE) + # color pallets option = A-H
  scale_shape_manual(values = c(20, 1), name=" ") + 
  guides(color = guide_legend(override.aes = list(label = ''))) +
  annotate("segment", x = 0, xend = (0.8833005), y = 5.374, yend = 0,
            colour = "darkturquoise", size = 1.2, linetype=1) +
  scale_x_continuous(labels = scales::percent) +
  geom_smooth(method="lm", se=F, fullrange=TRUE, size = 2, color = "black")

# calculate average slope
x <- lm(data = ds_long, chalkbrood ~ uboScore)
summary(x)


top <- ggplot(ds_long, aes(x=uboScore, y = (chalkbrood), 
                           linetype=as.character(chalk_type), shape = as.character(chalk_type))) +
  geom_point(size=6) +
  theme_bw(base_size = 20) + # size of the text and label ticks
  theme(legend.position = c(.75, .5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(0, .5)) +
  ylab(" ") + # y axis label
  scale_linetype_manual(values = c(1, 3), name=NULL, guide = FALSE) + # color pallets option = A-H
  scale_shape_manual(values = c(20, 1), name=NULL) + 
  guides(color = guide_legend(override.aes = list(label = ''))) +
  scale_y_continuous(limits = c(150, 200), breaks = seq(150, 200, by = 25))
top

plt <- plot_grid(top, bot, ncol = 1, rel_heights = c(.30, 1), align = "v")

plot_grid(plt + draw_label("UBeeO Score", x=0.55, y=  0, vjust=-.8, angle= 0, size = 20) +
            draw_label("Chalkbrood (cells/colony)", x=  0, y=0.55, vjust= 1.5, angle=90, size = 20))




# calculate average slope
x <- lm(data = ds_long, chalkbrood ~ uboScore)
summary(x)



xdf <- ds_long[ds_long$chalk_type == "White Chalk",]

ggplot(xdf, aes(x=uboScore, y = logChalk)) +
  geom_point(size=6) +
  ylab("Chalkbrood (cells/frame)") + # y axis label
  xlab("UBeeO Score") + # x axis label
  theme_minimal(base_size = 20) + # size of the text and label ticks
  theme(legend.position = c(.75, .9)) + # place the legend at the top
  coord_cartesian(xlim = c(0, .5)) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  #  annotate("segment", x = 0, xend = (229.2726/100), y = 2.239076, yend = 0,
  #           colour = "darkturquoise", size = 1.2, linetype=1) +
  scale_x_continuous(labels = scales::percent) +
  #  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                labels = trans_format("log10", math_format(10^.x))) + 
  geom_smooth(method="loess", se=F, fullrange=TRUE, size = 2, color = "black")



# glm for chalkbrood count data by ubo score and chalk type (poisson)\
ds_long$chalk_type <- factor(ds_long$chalk_type, levels = c("White Chalk", "Black Chalk Spores"))
ds_long$ID <- as.character(ds_long$Test_No)

# CHALKBROOD MODEL
mod <- glm(data = ds_long, chalkbrood ~ Percentage_UBO * chalk_type, family = poisson(link = "log"))
Anova(mod)
summary(mod)





#################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2022 Nosema
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#################################################################################

cleanDS_22 <- dplyr::select(cleanDS_22, year, beekeeper, yard, sampling_event, lab_ID, assay_score, sampling_date, varroa_count, 
                            nosema_count, FKB_percentile)


# select ubo testing periods
x <- cleanDS_22[cleanDS_22$assay_score >= 0,]
uboMap <- dplyr::select(x, lab_ID, assay_score)[complete.cases(dplyr::select(x, lab_ID, assay_score)),]
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
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1, color=UBO_Char))+
  labs(x="Sampling Month", y="Nosema Prevalence", color=" ", tag = "A") +
  scale_color_manual(values = c("#5071A0", "#E77624")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(guide = guide_axis(angle = 45))



# add factor data and make ubo a char
nosemaLoad_Sum <- nosemaLoad_Sum[!is.na(nosemaLoad_Sum$UBO_binary),]
nosemaLoad_Sum$Month <- factor(nosemaLoad_Sum$Month, levels = c("June", "July", "August", "Sept."))
nosemaLoad_Sum$UBO_Char <- ifelse(nosemaLoad_Sum$UBO_binary==1, "UBeeO High", "UBeeO Low")

contNos <-ggplot(nosemaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5) +
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Date", y="Nosema Load (spores/bee)", color=" ", tag = "B") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("#5071A0", "#E77624")) +
  scale_x_discrete(guide = guide_axis(angle = 45))

leg <-ggplot(nosemaLoad_Sum, aes(x=Month, y=mean, group=UBO_Char)) +
  geom_point(aes(color=UBO_Char), size=5) +
  geom_line(aes(color=UBO_Char), size=1.5) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "top") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = 0.1 ,color=UBO_Char))+
  labs(x="Sampling Date", y="Nosema Load (spores/bee)", color="UBeeO Status:", tag = "B") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("#5071A0", "#E77624"))
  

# make a multi panel plot
# get legend
grobs <- ggplotGrob(leg)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# plot figures
plt <- grid.arrange(nosPrev, contNos, ncol=2)
plot_grid(legend, plt, nrow = 2, rel_heights = c(.1, 1))

# prevalence model
mod2 <- glmer(data = cleanDS_22, Nosema_binary ~ UBO_binary * Month + (1 | yard), family = binomial(link="logit"))
Anova(mod2)
summary(mod2)

# load model
mod3 <- glmer(data = cleanDS_22, (1+nosema_count) ~ UBO_binary * Month + (1 | yard), family = Gamma(link="log"))
Anova(mod3)
summary(mod3)

# multcomp
cleanDS_22$inter <- interaction(cleanDS_22$UBO_binary, cleanDS_22$Month)
interMod <- glmer(data = cleanDS_22, (1+nosema_count) ~ inter + (1 | yard), family = Gamma(link="log"))
cht <- glht(interMod, linfct=mcp(inter = "Tukey"))
summary(cht, test = univariate())






##########################################################################
# Nosema Scatter plot

# make blocking var
cleanDS_22$block <- paste0(cleanDS_22$beekeeper, "_", cleanDS_22$yard)

# make month a factor
cleanDS_22_noMonth <- cleanDS_22[!is.na(cleanDS_22$Month),]
cleanDS_22_noMonth$Month <- factor(cleanDS_22_noMonth$Month, levels = c("June", "July", "August", "Sept."))

ggplot(cleanDS_22_noMonth, aes(x=assay_score, y=((nosema_count*4000000)/80), color=Month, shape=Month)) +
  geom_point(size=4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 1.9) +
  theme_minimal(base_size = 20) +
  theme(legend.position = c(.8,.8)) +
  labs(x="UBeeO Score", y="Nosema Load (spores/bee)", color="Month") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("#9E519F","#519E9A", "#5071A0", "#E77624")) 


nosPos <- cleanDS_22_noMonth[cleanDS_22_noMonth$nosema_count>0, ]


# facet wrap version
ggplot(nosPos, aes(x=assay_score, y=1+((nosema_count*4000000)/80))) +
  geom_point(size=4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "#619B50") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="UBeeO Score", y="Nosema Load (spores/bee)", color="Month") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 2),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_continuous(labels = scales::percent)



# month as a covariate
y <- glmer(data = nosPos, nosema_count ~ assay_score * Month + (1| block), family = Gamma(link = "log"))
Anova(y)
summary(y)

# split month
nosemaSplit <- split(nosPos, nosPos$Month)

# regression for nosema
mod <- glmer(data = nosemaSplit$June, nosema_count+1 ~ assay_score + (1|block), family = Gamma(link = "log"))
mod1 <- glmer(data = nosemaSplit$July, nosema_count+1 ~ assay_score + (1|block), family = Gamma(link = "log"))
mod2 <- glmer(data = nosemaSplit$August, nosema_count+1 ~ assay_score + (1|block), family = Gamma(link = "log"))
mod3 <- glmer(data = nosemaSplit$Sept., nosema_count+1 ~ assay_score + (1|block), family = Gamma(link = "log"))

Anova(mod)
Anova(mod1)
Anova(mod2)
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
virus_long$ubo_binary_june <- ifelse(virus_long$June.UBO > 60, "UBeeO High", "UBeeO Low")
virus_long$ubo_binary_august <- ifelse(virus_long$August.UBO > 60, "UBeeO High", "UBeeO Low")
virus_long$mite_binary_june <-  ifelse(virus_long$June.Mite > 0, 1, 0) 
virus_long$mite_binary_august <- ifelse(virus_long$August.Mite > 0, 1, 0)
virus_long$logLoad <- log10(virus_long$virus_load + 1)

virus_long <- virus_long[!is.na(virus_long$ubo_binary_august),]

# grouped boxplot
vld <- ggplot(data = virus_long, aes(x = virus, y = logLoad, fill = ubo_binary_june)) + 
  geom_boxplot() +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="Virus", y="Log10(load) (copies/ul)", fill = "UBeeO Status:", tag = "B") +
  scale_fill_manual(values = c("#5071A0", "#E77624")) +
  scale_x_discrete(guide = guide_axis(angle = 45))+
  annotate("text", x = c(1:6), y = c(9.6, 8,10,9,8.7,8.4), label = c("ns", "*", "**", "**", "ns", "ns"), size = 6)

vld


# virus binary summary
virusPrevSum <- virus_long %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(virus, ubo_binary_june) %>% # pick variables to group by
  dplyr::summarise(
    mean = mean(virus_binary, na.rm=T), # mean
    n = length(virus_binary),
    a = sum(virus_binary, na.rm = T)+1,
    b = n - a + 1,
    lower = qbeta(.025, shape1 = a, shape2 = b),
    upper = qbeta(.975, shape1 = a, shape2 = b),
  ) 

virusPrevSum <- virusPrevSum[complete.cases(virusPrevSum),]
virusPrevSum[virusPrevSum$mean==0,]$lower <- NA
virusPrevSum[virusPrevSum$mean==0,]$upper <- NA


vP <- ggplot(virusPrevSum, aes(x=virus, y=mean, fill=ubo_binary_june)) + 
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="Virus", y="Virus Prevalence", fill = "UBeeO Status:", tag = "A") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("#5071A0", "#E77624")) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  annotate("text", x = c(1:6), y = c(1.06, .9,.98,.88,1.05,1.06), label = c("ns", "ns", "**", "ns", "**", "ns"), size = 6)
vP

# extract legend
leg <- ggplot(virusPrevSum, aes(x=virus, y=mean, fill=ubo_binary_june)) + 
  theme_minimal(base_size = 20) +
  theme(legend.position = "top") +
  labs(x="Virus", y="Virus Prevalence", fill = "UBeeO Status:", tag = "A") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("#5071A0", "#E77624")) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(guide = guide_axis(angle = 45))

# make data wide
virus_wide <- dplyr::select(virus_long, -virus_binary, -ubo_binary_june, -ubo_binary_august,-mite_binary_june, -virus_load)
virus_wide <- virus_wide %>% pivot_wider(names_from = virus, values_from = logLoad)
virus_wide$ubo_binary_june <- ifelse(virus_wide$June.UBO > 60, "UBeeO +", "UBeeO -")
virus_wide$ubo_binary_august <- ifelse(virus_wide$August.UBO > 60, "UBeeO +", "UBeeO -")





# get legend
grobs <- ggplotGrob(leg)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# plot figures
plt <- grid.arrange(vP, vld, ncol=2)
plot_grid(legend, plt, nrow = 2, rel_heights = c(.1, 1))




ggplot(virus_long, aes(x=(June.UBO/100), y=(virus_load+1))) +
  geom_point(size=4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "#619B50") +
  theme_minimal(base_size = 20) +
  theme(legend.position = c(.8,.8)) +
  labs(x="UBeeO Score", y="Virus Load (copies/ul)", ) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(vars(virus)) +
  scale_x_continuous(labels = scales::percent, guide = guide_axis(angle = 45))



# split dataframe by virus and run regression on log virus load
splitVirus <- split(virus_long, virus_long$virus)


# VIRUS PREV - UBO BINARY
mod6 <- glm(data = splitVirus$BQCV, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))
mod7 <- glm(data = splitVirus$DWV.A, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))
mod8 <- glm(data = splitVirus$DWV.B, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))
mod9 <- glm(data = splitVirus$IAPV, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))
mod10 <- glm(data = splitVirus$LSV, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))
mod11 <- glm(data = splitVirus$SBV, virus_binary ~ ubo_binary_june, family = binomial(link = "logit"))

Anova(mod6)
Anova(mod7)
Anova(mod8)
Anova(mod9)
Anova(mod10)
Anova(mod11)

# VIRUS LOAD - UBO BINARY
mod12 <- glm(data = splitVirus$BQCV, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))
mod13 <- glm(data = splitVirus$DWV.A, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))
mod14 <- glm(data = splitVirus$DWV.B, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))
mod15 <- glm(data = splitVirus$IAPV, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))
mod16 <- glm(data = splitVirus$LSV, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))
mod17 <- glm(data = splitVirus$SBV, (virus_load+1) ~ ubo_binary_june, family = Gamma(link = "identity"))

Anova(mod12)
Anova(mod13)
Anova(mod14)
Anova(mod15)
Anova(mod16)
Anova(mod17)


# VIRUS LOAD - continuous UBO:
mod <- lm(data = splitVirus$BQCV, logLoad ~ June.UBO)
mod2 <- lm(data = splitVirus$DWV.A, logLoad ~ June.UBO)
mod1 <- lm(data = splitVirus$DWV.B, logLoad ~ June.UBO)
mod3 <- lm(data = splitVirus$IAPV, logLoad~ June.UBO)
mod4 <- lm(data = splitVirus$LSV, logLoad ~ June.UBO)
mod5 <- lm(data = splitVirus$SBV, logLoad ~ June.UBO)

Anova(mod)
Anova(mod2)
Anova(mod1)
Anova(mod3)
Anova(mod4)
Anova(mod5)

# NOTE: consider benjamini-hochberg to reduce FDR (have to think about hypothesis a bit more to check groupings)




