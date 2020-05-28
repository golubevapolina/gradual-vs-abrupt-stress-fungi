
library(car)
library(FSA)
library(ggplot2)
library(gridExtra)
library(ggpubr)
#setwd("/Users/polina/Desktop/R")
setwd("/Users/polina/Desktop/R.files/Liquid medium experiment/ANOVA")

fungi.raw<-read.csv (file="Biomass+Copper.csv",header=TRUE)
fungi.raw<-na.omit(fungi.raw)

# Check that each variable is in the rigth class
for (i in 1:length(colnames(fungi.raw))) print(c(colnames(fungi.raw[i]),class(fungi.raw[,i])))

#====Create subsets====
#This loop creates a subset of data for each species and saves it in a data frame with the name "biomass_name.of.sp."
#It also creates a list with the names of those subsets, named "subset_list"
subset_list<-vector()
for (i in 1:length(levels(fungi.raw$Sp.))) {
  assign(paste("biomass",levels(fungi.raw$Sp.)[i],sep="_"),
         subset(fungi.raw,fungi.raw$Sp.==levels(fungi.raw$Sp.)[i]))
  subset_list[i]<-paste("biomass",levels(fungi.raw$Sp.)[i],sep="_")
}
subset_list
#======Assumptions for ANOVA======#
assumptions<-matrix(NA,length(subset_list),4,
                    dimnames = list(levels(fungi.raw$Sp.)
                    ,c("p_normality of residuals","normality of residuals"
                    ,"p_homoscedasticity","homoscedasticity")))

for (i in 1:length(subset_list)) {
  a<-get(subset_list[i])
  assumptions[i,1]<- as.numeric(unlist(shapiro.test(residuals(aov(a$Biomass~a$Treatment))))[2])
  assumptions[i,3]<- as.numeric(unlist(leveneTest(a$Biomass~a$Treatment))[3])
  qqnorm(residuals(aov(a$Biomass~a$Treatment)))
  qqline(residuals(aov(a$Biomass~a$Treatment)))}

for (i in 1:length(subset_list)) {
  if (assumptions[i,1]<=0.05) {assumptions[i,2]<-"no"}
  else assumptions[i,2]<-"yes"}
for (i in 1:length(subset_list)) {
  if (assumptions[i,3]<=0.05) {assumptions[i,4]<-"no"}
  else assumptions[i,4]<-"yes"}

#====Analysis====
#Creates a matrix where we will store our results, names rows with species
results<-matrix(NA,length(subset_list),7,dimnames = list(levels(fungi.raw$Sp.),
                                                        c("mean_C","mean_A","mean_G", "p_value ANOVA",
                                                       "C-A","G-A","G-C")))
                                                          
for (i in 1:length(subset_list)) {
  a<-get(subset_list[i])
  results[i,1]<-as.numeric(mean(a$Biomass[a$Treatment=="C"]))
  results[i,2]<-as.numeric(mean(a$Biomass[a$Treatment=="A"]))
  results[i,3]<-as.numeric(mean(a$Biomass[a$Treatment=="G"]))
  results[i,4]<-as.numeric(summary(aov(a$Biomass~a$Treatment))[[1]]$"Pr(>F)"[1])
  b<-(TukeyHSD(aov(a$Biomass~a$Treatment)))
  results[i,5]<-as.numeric(b$`a$Treatment`[1,4])
  results[i,6]<-as.numeric(b$`a$Treatment`[2,4])
  results[i,7]<-as.numeric(b$`a$Treatment`[3,4])}

results
for (i in 1:length(colnames(results))) print(c(colnames(results[i]),class(results[,i])))


myplots <- list()  # new empty list
for (i in 1:length(subset_list))
  local({
    i <- i
    a<-get(subset_list[i]) 
    Species<-a$Sp.
    Biomass<-a$Biomass
    Treatment_type<-a$Treatment
    bp<-ggplot(data=a, aes(x=Treatment_type, y =Biomass,fill=factor(Treatment_type)))+
      geom_boxplot(outlier.shape = NA)+
      xlab("Treatment") + ylab("Biomass production")+
      ggtitle(Species) +
      theme(legend.position="none")
    print(i)
    print(bp)
    myplots_Cu[[i]] <<- bp  # add each plot into plot list
  })



