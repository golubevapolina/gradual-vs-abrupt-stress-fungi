library(car)
library(FSA)
library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("/Users/polina/Desktop/Frontiers/R") 
means<-read.csv (file="Means_Frontiers.csv",header=TRUE)
means$Trt<-factor(means$Trt,levels = c ("Control", "Gradual", "Abrupt"))

############################################################################################################
# Creating subsets

B_Ctrl<-as.numeric(means$Mean_B[means$Trt=="Control"]) #subset Biomass for control group
G_ES_B<-as.numeric(means$ES_B[means$Trt=="Gradual"]) #ES gradual Biomass
A_ES_B<-as.numeric(means$ES_B[means$Trt=="Abrupt"]) #ES abrupt Biomass
G_ES_Cu<-as.numeric(means$ES_Cu[means$Trt=="Gradual"]) #ES gradual Copper
A_ES_Cu<-as.numeric(means$ES_Cu[means$Trt=="Abrupt"]) #ES abrupt Biomass

d<-cbind.data.frame(B_Ctrl, G_ES_B, A_ES_B, G_ES_Cu, A_ES_Cu) # Create data frame with these subsets
d
d[,6]<-(d[,3]-d[,2])#Add column of ES.difference (A-G)
colnames(d)<-c("Biomass production in Control group","ES Biomass under Gradual treatment", 
               "ES Biomass under Abrupt treatment", "ES Copper accumulation under Gradual treatment",
               "ES Copper accumulation under Abrupt treatment",
               "Difference in ES Biomass between the treatments")
############################################################################################################
#Correlation tests

B_G<-cor.test(G_ES_B,B_Ctrl, method = "pearson",use = "complete.obs")
B_G

B_A<-cor.test(A_ES_B,B_Ctrl, method = "pearson",use = "complete.obs")
B_A

Cu_G<-cor.test(G_ES_Cu,B_Ctrl, method = "pearson",use = "complete.obs")
Cu_G

Cu_A<-cor.test(A_ES_Cu,B_Ctrl, method = "pearson",use = "complete.obs")
Cu_A

res_ES_Diff<-cor.test(d[,1],d[,6] , method = "pearson",use = "complete.obs")
res_ES_Diff
############################################################################################################
#Plots for correlation Biomass for control group and ES

i<-ggplot(d, aes(x=B_Ctrl, y = G_ES_B))+
  geom_point()+
  geom_smooth(method=lm,alpha= 0.3)+  
  stat_cor(method="pearson")
i

ii<-ggplot(d, aes(x=B_Ctrl, y = A_ES_B))+
  geom_point()+
  geom_smooth(method=lm,alpha= 0.3)+  
  stat_cor(method="pearson")
ii

iii<-ggplot(d, aes(x=B_Ctrl, y = G_ES_Cu))+
  geom_point()+
  geom_smooth(method=lm,alpha= 0.3)+  
  stat_cor(method="pearson")
iii

iiii<-ggplot(d, aes(x=B_Ctrl, y = A_ES_Cu))+
  geom_point()+
  geom_smooth(method=lm,alpha= 0.3)+  
  stat_cor(method="pearson")
iiii

gridExtra::grid.arrange(i,ii,iii,iiii)

############################################################################################################
#Plot for Correlation between Biomass in control group and ES Difference

ES_diff<-ggplot(d, aes(x= `Biomass production in Control group`, 
                 y =`Difference in ES Biomass between the treatments`))+
  xlim(30,155)+
  geom_point(size=3)+
  geom_smooth(method=lm,alpha= 0.3, aes())+ 
  stat_cor(aes(),method="pearson", label.x = 110.5, label.y = 17.5, size=5)+
  theme(text = element_text(size=15))+
  theme(legend.position="none")+
  theme_bw()

ES_diff


############################################################################################################
#Plot for Correlation between Biomass in control group and ES Biomass
setwd("/Users/polina/Desktop/R.files/Liquid medium experiment/ES~GR") #laptop
GR_data<-read.csv(file="GR-ES_Frontiers.csv",header=TRUE)
colnames(GR_data)<-c("Isolate",
                     "Treatment",
                     "Biomass", 
                     "Biomass production in Control group",
                     "ES Biomass")


ES_B_GR<-ggplot(GR_data, aes(x=`Biomass production in Control group`,
                             y = `ES Biomass`, color= `Treatment`), show.legend = FALSE)+
  xlim(30,155)+
  geom_point(aes(color=`Treatment`), size=3, show.legend = FALSE) +
  #geom_text(aes(label=Species),hjust=0, vjust=0)+
  geom_smooth(method=lm,alpha= 0.3, aes(fill=`Treatment`),show.legend = FALSE)+  
  scale_fill_manual(values=c("#f8766d","#00bfc4"))+  
  scale_color_manual(values=c("#f8766d","#00bfc4"))+
  stat_cor(aes(color=`Treatment`),method="pearson", label.x = 110.5,  size=5, show.legend = FALSE)+
  #theme(text = element_text(size=15), legend.position="none")+
  theme_bw()

ES_B_GR

gridExtra::grid.arrange(ES_B_GR,ES_diff, nrow=2)

