#####################################################################
#
#                              ANOVA
#
#####################################################################
library(car)
library(FSA)
library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("/Users/polina/Desktop/Frontiers/R") 
fungi.raw<-read.csv (file="Biomass+Copper.csv",header=TRUE)
fungi.raw<-na.omit(fungi.raw)
fungi.raw$Treatment<-factor(fungi.raw$Treatment,levels = c ("C", "G", "A"))

###########---Assumptions for ANOVA for Biomass---###########
assumptions<-matrix(nrow = 1, ncol = 2)
colnames(assumptions)<-c("p_normality of residuals","p_homoscedasticity")
assumptions[1,1]<-as.numeric(shapiro.test(residuals(aov(fungi.raw$Biomass~fungi.raw$Treatment)))[2])
assumptions[1,2]<-as.numeric(leveneTest(fungi.raw$Biomass~fungi.raw$Treatment)[1,3])
qqnorm(residuals(aov(fungi.raw$Biomass~fungi.raw$Treatment)))
qqline(residuals(aov(fungi.raw$Biomass~fungi.raw$Treatment)))

#  p_normality of residuals     p_homoscedasticity
#      3.377735e-07                 6.904303e-11

###########---ANOVA for Biomass---###########
B_Ctrl<-fungi.raw$Biomass[fungi.raw$Treatment=="C"] #subset for control group
B_Trt_G<-fungi.raw$Biomass[fungi.raw$Treatment=="G"] #subset for gradualgroup
B_Trt_A<-fungi.raw$Biomass[fungi.raw$Treatment=="A"] #subset for abrupt group

results_B<-matrix(nrow = 7, ncol = 1)
row.names(results_B)<-c("mean_C","mean_A","mean_G", "p_value ANOVA", "C-G", "C-A", "G-A")
results_B[1,1]<-as.numeric(mean(B_Ctrl))
results_B[2,1]<-as.numeric(mean(B_Trt_G))
results_B[3,1]<-as.numeric(mean(B_Trt_A))
results_B[4,1]<-as.numeric(summary(aov(fungi.raw$Biomass~fungi.raw$Treatment))[[1]]$"Pr(>F)"[1])
b<-(TukeyHSD(aov(fungi.raw$Biomass~fungi.raw$Treatment)))
results_B[5,1]<-as.numeric(b$`fungi.raw$Treatment`[1,4])
results_B[6,1]<-as.numeric(b$`fungi.raw$Treatment`[2,4])
results_B[7,1]<-as.numeric(b$`fungi.raw$Treatment`[3,4])
results_B


###########---Assumptions for ANOVA for Cu uptake---###########   
assumptions2<-matrix(nrow = 1, ncol = 2)
colnames(assumptions2)<-c("p_normality of residuals","p_homoscedasticity")
assumptions2[1,1]<-as.numeric(shapiro.test(residuals(aov(fungi.raw$Copper~fungi.raw$Treatment)))[2])
assumptions2[1,2]<-as.numeric(leveneTest(fungi.raw$Copper~fungi.raw$Treatment)[1,3])
qqnorm(residuals(aov(fungi.raw$Copper~fungi.raw$Treatment)))
qqline(residuals(aov(fungi.raw$Copper~fungi.raw$Treatment)))
assumptions2

###########--- ANOVA for Cu accumulation---###########
Cu_Ctrl<-fungi.raw$Copper[fungi.raw$Treatment=="C"] #subset for control group
Cu_Trt_G<-fungi.raw$Copper[fungi.raw$Treatment=="G"] #subset for gradual group
Cu_Trt_A<-fungi.raw$Copper[fungi.raw$Treatment=="A"] #subset for abrupt group

results_Cu<-matrix(nrow = 7, ncol = 1)
row.names(results_Cu)<-c("mean_C","mean_A","mean_G", "p_value ANOVA", "C-G", "C-A", "G-A") 
results_Cu[1,1]<-as.numeric(mean(Cu_Ctrl))
results_Cu[2,1]<-as.numeric(mean(Cu_Trt_G))
results_Cu[3,1]<-as.numeric(mean(Cu_Trt_A))
results_Cu[4,1]<-as.numeric(summary(aov(fungi.raw$Copper~fungi.raw$Treatment))[[1]]$"Pr(>F)"[1])
c<-(TukeyHSD(aov(fungi.raw$Copper~fungi.raw$Treatment)))
results_Cu[5,1]<-as.numeric(c$`fungi.raw$Treatment`[1,4])
results_Cu[6,1]<-as.numeric(c$`fungi.raw$Treatment`[2,4])
results_Cu[7,1]<-as.numeric(c$`fungi.raw$Treatment`[3,4])
results_Cu

