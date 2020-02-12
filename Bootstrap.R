#####################################################################
#  Nonparametric test for the difference in the effect sizes
#    between Gradual and Abrupt treatments
#    considering the variability in control as well    
#                               
#    2019.03.13 Masahiro Ryo

#####################################################################

library(ggplot2)
library(ggridges)

# devtools::install_github("thomasp85/patchwork") 
library(patchwork)

#devtools::install_github("lwjohnst86/ggepi")
library(ggepi)

###########---Functions ---###########
#----------------------------------------------
# Distributions of effect sizes
#----------------------------------------------
#
#  This function, ES_distribution, calculates ES values for all possible replicate combinations
#  e.g. if #s of replicate are 5 for control and 4 for treatment, 20 combinations are made.
#  thereby it can consider the variabilities in both control and treatment.
#  It generates The ES distribution for each of Gradual-Ctrl pair and Abrupt-Ctrl pair.
#

ES_bootstrap = function(response, data, n_iter = 10){
  
  #_______________________________________________________  
  # Part 1. Basic set-up
  
  # this stores results & information
  output = list()
  
  # row numbers and # of replicates for each treatment
  # not all of them are 5 replicates and some show NA
  id_ctrl = which(data["Treatment"]=="C" & !is.na(data[response]))
  id_trtG = which(data["Treatment"]=="G" & !is.na(data[response])) 
  id_trtA = which(data["Treatment"]=="A" & !is.na(data[response]))
  
  #_______________________________________________________  
  # Part 2. ES value calculation for G-C & A-C pairs, and the difference. Using bootstrap.
  
  bs = c()  
  
  for(k in 1:n_iter){
    
    # calculating the difference in the median values between treatment and control
    trt_G  = mean(sample(data[id_trtG,  response], length(id_trtG),  replace = T))
    trt_A  = mean(sample(data[id_trtA,  response], length(id_trtA),  replace = T))
    ctrl =   mean(sample(data[id_ctrl,  response], length(id_ctrl),  replace = T))
    
    ### ES is prepresented by the difference from CTRL
    effectG = (trt_G - ctrl)
    effectA = (trt_A - ctrl)
    deltaAG = effectA - effectG
    bs = rbind(bs, c(effectG, effectA, deltaAG))
  }
  colnames(bs) = c("effectG", "effectA", "effectA-G")
  
  #_______________________________________________________  
  # Part 3. Statiatical summary & test
  
  # NHST test  
  p = apply(bs, 2, function(x) min(length(which(x>0))/n_iter, 1 - length(which(x>0))/n_iter))  
  
  # stats summary
  output[["summary"]] = rbind(apply(bs, 2, function(x) quantile(x, c(.05, .50, .95))),
                              p)
  colnames(output[["summary"]]) = c("effectG", "effectA", "effectA-G")
  output[["ES_bootstrap"]] = bs
  
  return(output)
}

###########--- Main part ###########

# Firstly it applies the function above to all the species.
# Secondly, it visualizes the ES distributions using ggplot2 and ggridges 

#_______________________________________________________  
# Part 1. Basic set-up
# reading the dataset

setwd("/Users/polina/Desktop/Frontiers/R")
df =read.csv(file="Biomass+Copper.csv", sep = "," )

# extracting the names of species
sp_list = sapply(unlist(unique(df["Sp."])), as.character)

response = "Biomass"
ES = list()
df_plot = numeric()
stats_summary = matrix(NA, ncol=4,nrow=length(sp_list))
colnames(stats_summary) = c("p", "mean_difference_G_minus_A", "sd_difference", "n_perm")
rownames(stats_summary) = sp_list

#_______________________________________________________  
# Part 2. Apply the function for all species

n_iter = 1000

# apply the function for all species
for(i_sp in sp_list){
  ES[[i_sp]] = ES_bootstrap(data=df[which(df$Sp.==i_sp),], response=response, n_iter = n_iter)
}

ES[["global"]] = ES_bootstrap(data=df, response=response, n_iter = n_iter)

df.plot   = c()
df.plotbs = c() 
for(i_sp in c(sp_list,"global")){
  df.plot = rbind(df.plot,
                  data.frame(
                    species   = rep(i_sp, 3),
                    treat     = c("Gradual", "Abrupt", "Difference"), 
                    t(data.frame(ES[[i_sp]]["summary"]))
                  )
  )
  
  df.plotbs = rbind(df.plotbs,
                    data.frame(
                      species = rep(i_sp, n_iter*3),
                      treat   = rep(c("Gradual", "Abrupt", "Difference"), each=n_iter),
                      values  = unlist(ES[[i_sp]]["ES_bootstrap"])
                    )
  )
}

rownames(df.plot) = c()

write.csv(df.plot, file = "statistical summary.csv")

###########--- ggplot  for Biomass---###########

g1 = 
  ggplot(df.plotbs[df.plotbs$treat!="Difference", ], aes(y=species, x=values, fill = treat)) +
  xlab("Effect size [Treatment - Control]") +
  ylab("Isolates") +
  theme(legend.position = c(0.8,0.1)) +
  scale_fill_manual(values=c("plum3","turquoise3"))+   
  geom_density_ridges(rel_min_height = 0.001,
                      scale = 0.7, alpha = .6, color="#000000",
                      jittered_points = F) +
  geom_estci(data=df.plot[df.plot$treat=="Gradual", ], 
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), color = "turquoise3",
             size=0.7, ci.linesize = 0.9, alpha=0.6, position=position_nudge(y = -0.2)) +
  geom_estci(data=df.plot[df.plot$treat=="Abrupt", ], 
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), colour = "plum3",
             size=0.7, ci.linesize = 0.9, alpha=0.7, position=position_nudge(y = -0.1))+
  theme_bw()
g1

g2 = 
  ggplot(df.plotbs[df.plotbs$treat=="Difference", ], aes(y=species, x=values)) +
  xlab("Difference: ES[Abrupt] - ES[Gradual]") +
  ylab("Isolates") +
  theme(legend.position = "none") +
  geom_density_ridges(rel_min_height = 0.001,
                      scale = 0.7, alpha = .7, color="#000000", fill="slateblue3",
                      jittered_points = F) +
  geom_estci(data=df.plot[df.plot$treat=="Difference", ],
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), colour = "slateblue4", alpha=0.7,
             size=0.7, ci.linesize = 0.9, position=position_nudge(y = -0.2))+
  theme_bw()
g2

g1 | g2

#Repeat the same for Copper accumulation
response = "Copper"
ES = list()
df_plot = numeric()
stats_summary = matrix(NA, ncol=4,nrow=length(sp_list))
colnames(stats_summary) = c("p", "mean_difference_G_minus_A", "sd_difference", "n_perm")
rownames(stats_summary) = sp_list

n_iter = 1000

# apply the function for all species
for(i_sp in sp_list){
  ES[[i_sp]] = ES_bootstrap(data=df[which(df$Sp.==i_sp),], response=response, n_iter = n_iter)
}

ES[["global"]] = ES_bootstrap(data=df, response=response, n_iter = n_iter)

df.plot   = c()
df.plotbs = c() 
for(i_sp in c(sp_list,"global")){
  df.plot = rbind(df.plot,
                  data.frame(
                    species   = rep(i_sp, 3),
                    treat     = c("Gradual", "Abrupt", "Difference"), 
                    t(data.frame(ES[[i_sp]]["summary"]))
                  )
  )
  
  df.plotbs = rbind(df.plotbs,
                    data.frame(
                      species = rep(i_sp, n_iter*3),
                      treat   = rep(c("Gradual", "Abrupt", "Difference"), each=n_iter),
                      values  = unlist(ES[[i_sp]]["ES_bootstrap"])
                    )
  )
}

rownames(df.plot) = c()

write.csv(df.plot, file = "statistical summary.csv")

###########---ggplot for Copper  ---###########	                                      	          
g3 = 
  ggplot(df.plotbs[df.plotbs$treat!="Difference", ], aes(y=species, x=values, fill = treat)) +
  xlab("Effect size [Treatment - Control]") +
  ylab("Isolates") +
  theme(legend.position = c(0.8,0.1)) +
  scale_fill_manual(values=c("plum3","turquoise3"))+   
  geom_density_ridges(rel_min_height = 0.001,
                      scale = 0.7, alpha = .6, color="#000000",
                      jittered_points = F) +
  geom_estci(data=df.plot[df.plot$treat=="Gradual", ], 
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), color = "turquoise3",
             size=0.7, ci.linesize = 0.9, alpha=0.6, position=position_nudge(y = -0.2)) +
  geom_estci(data=df.plot[df.plot$treat=="Abrupt", ], 
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), colour = "plum3",
             size=0.7, ci.linesize = 0.9, alpha=0.7, position=position_nudge(y = -0.1))+
  theme_bw()
g3

g4 = 
  ggplot(df.plotbs[df.plotbs$treat=="Difference", ], aes(y=species, x=values)) +
  xlab("Difference: ES[Abrupt] - ES[Gradual]") +
  ylab("Isolates") +
  theme(legend.position = "none") +
  geom_density_ridges(rel_min_height = 0.001,
                      scale = 0.7, alpha = .7, color="#000000", fill="slateblue3",
                      jittered_points = F) +
  geom_estci(data=df.plot[df.plot$treat=="Difference", ],
             aes(x = X50., y = species, xmin=X5., xmax=X95., xintercept=0), colour = "slateblue4", alpha=0.7,
             size=0.7, ci.linesize = 0.9, position=position_nudge(y = -0.2))+
  theme_bw()
g4

g3 | g4
