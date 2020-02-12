#####################################################################
#
#  		    Correlation between Biomass and Copper accumulation
#
#####################################################################
setwd("/Users/polina/Desktop/Frontiers/R") 
means2<- read.csv (file = "Means_Frontiers.csv", header=TRUE)
means2<-na.omit(means2)

colnames(means2)<- c("Species", "Treatment", "Biomass_production","ES_Biomass", "Copper_uptake","ES_Copper_uptake")

###########---Cor.coeff. for Gradual and Abrupt---###########
G_ES_B<-as.numeric(means2$ES_Biomass[means2$Treatment=="Gradual"])
G_ES_Cu<-as.numeric(means2$ES_Copper_uptake[means2$Treatment=="Gradual"])
res_ES_G<-cor.test(G_ES_B, G_ES_Cu, method = "pearson",use = "complete.obs")
res_ES_G

A_ES_B<-as.numeric(means2$ES_Biomass[means2$Treatment=="Abrupt"])
A_ES_Cu<-as.numeric(means2$ES_Copper_uptake[means2$Treatment=="Abrupt"])
res_ES_A<-cor.test(A_ES_B, A_ES_Cu, method = "pearson",use = "complete.obs")
res_ES_A

cor.results<-matrix(nrow = 2, ncol = 4)
colnames(cor.results)<-c("cor.coef","p-value","conf.int 0.05", "conf.int 0.95")
row.names(cor.results)<-c("Gradual","Abrupt")
cor.results[1,1]<-as.numeric(res_ES_G$estimate)
cor.results[2,1]<-as.numeric(res_ES_A$estimate)
cor.results[1,2]<-as.numeric(res_ES_G$p.value)
cor.results[2,2]<-as.numeric(res_ES_A$p.value)
cor.results[1,3]<-(res_ES_G$conf.int[1])
cor.results[2,3]<-(res_ES_A$conf.int[1])
cor.results[1,4]<-(res_ES_G$conf.int[2])
cor.results[2,4]<-(res_ES_A$conf.int[2])
cor.results

###########---Plot---###########
library(ggplot2)
library(ggpubr)

io<-ggplot(means2, aes(x= ES_Biomass, y = ES_Copper_uptake, color= Treatment))+
  geom_point(aes(color=Treatment)) +geom_text(aes(label=Species),hjust=0, vjust=0)+
  geom_smooth(method=lm,alpha= 0.3, aes(fill=Treatment))+  
  scale_fill_manual(values=c("#f8766d","#00bfc4"))+  
  scale_color_manual(values=c("#f8766d","#00bfc4"))+
  theme(legend.position=c(0.9,0.15))+
  stat_cor(aes(color=Treatment),method="pearson", label.x = 0)
io

# Shapiro-Wilk normality test for ES_Copper
shapiro.test(G_ES_Cu) 
# Shapiro-Wilk normality test for ES_Biomass
shapiro.test(G_ES_B)
# Shapiro-Wilk normality test for ES_Copper
shapiro.test(A_ES_Cu) 
# Shapiro-Wilk normality test for ES_Biomass
shapiro.test(A_ES_B) 
