library(car)
library(FSA)
library(ggplot2)
library(gridExtra)
setwd("/Users/polina/Desktop/R")
fungi.raw<-read.csv (file="Biomass+Copper.csv",header=TRUE)
fungi.raw<-na.omit(fungi.raw)

# Check that each variable is in the rigth class
for (i in 1:length(colnames(fungi.raw))) print(c(colnames(fungi.raw[i]),class(fungi.raw[,i])))

#====Create subsets====
#This loop creates a subset of data for each species and saves it in a data frame with the name "biomass_name.of.sp."
#It also creates a list with the names of those subsets, named "subset_list"
subset_list<-vector()
for (i in 1:length(levels(fungi.raw$Sp.))) {
  assign(paste("Cu_uptake",levels(fungi.raw$Sp.)[i],sep="_"),
         subset(fungi.raw,fungi.raw$Sp.==levels(fungi.raw$Sp.)[i]))
  subset_list[i]<-paste("Cu_uptake",levels(fungi.raw$Sp.)[i],sep="_")
}
subset_list

#======Assumptions for ANOVA======#
assumptions_Cu<-matrix(NA,length(subset_list),4,
                    dimnames = list(levels(fungi.raw$Sp.)
                                    ,c("p_normality of residuals","normality of residuals"
                                       ,"p_homoscedasticity","homoscedasticity")))

for (i in 1:length(subset_list)) {
  a<-get(subset_list[i])
  assumptions_Cu[i,1]<- as.numeric(unlist(shapiro.test(residuals(aov(a$Copper~a$Treatment))))[2])
  assumptions_Cu[i,3]<- as.numeric(unlist(leveneTest(a$Copper~a$Treatment))[3])
  qqnorm(residuals(aov(a$Copper~a$Treatment)))
  qqline(residuals(aov(a$Copper~a$Treatment)))}

for (i in 1:length(subset_list)) {
  if (assumptions_Cu[i,1]<=0.05) {assumptions_Cu[i,2]<-"no"}
  else assumptions_Cu[i,2]<-"yes"}
for (i in 1:length(subset_list)) {
  if (assumptions_Cu[i,3]<=0.05) {assumptions_Cu[i,4]<-"no"}
  else assumptions_Cu[i,4]<-"yes"}
assumptions_Cu

#====Analysis====
#Creates a matrix where we will store our results, names rows with species
results_Cu<-matrix(NA,length(subset_list),7,dimnames = list(levels(fungi.raw$Sp.),
                                                         c("mean_C_Cu","mean_A_Cu","mean_G_Cu", "p_value ANOVA_Cu",
                                                           "C-A","G-A","G-C")))


for (i in 1:length(subset_list)) {
  a<-get(subset_list[i])
  results_Cu[i,1]<-as.numeric(mean(a$Copper[a$Treatment=="C"]))
  results_Cu[i,2]<-as.numeric(mean(a$Copper[a$Treatment=="A"]))
  results_Cu[i,3]<-as.numeric(mean(a$Copper[a$Treatment=="G"]))
  results_Cu[i,4]<-as.numeric(summary(aov(a$Copper~a$Treatment))[[1]]$"Pr(>F)"[1])
  b<-(TukeyHSD(aov(a$Copper~a$Treatment)))
  results_Cu[i,5]<-as.numeric(b$`a$Treatment`[1,4])
  results_Cu[i,6]<-as.numeric(b$`a$Treatment`[2,4])
  results_Cu[i,7]<-as.numeric(b$`a$Treatment`[3,4])}


#Check results
results_Cu
for (i in 1:length(colnames(results_Cu))) print(c(colnames(results_Cu[i]),class(results_Cu[,i])))

#====Graphs===#

myplots_Cu <- list()  # new empty list
for (i in 1:length(subset_list))
  local({
    i <- i
    a<-get(subset_list[i]) 
    Species<-a$Sp.
    Copper_content<-a$Copper
    Treatment_type<-a$Treatment
    bp<-ggplot(data=a, aes(x=Treatment_type, y =Copper_content,fill=factor(Treatment_type)))+
    geom_boxplot(outlier.shape = NA)+
    xlab("Treatment") + ylab("Copper content (mg/g DW)")+
    scale_y_continuous(limits = c(0, 40))+
    ggtitle(Species) +
    theme(legend.position="none")
    print(i)
    print(bp)
    myplots_Cu[[i]] <<- bp  # add each plot into plot list
  })



