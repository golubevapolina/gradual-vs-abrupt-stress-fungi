#####################################################################
#
#  				     Phylogenetic signal  
#
#####################################################################

##Read the new tree
library(ape)
setwd("/Users/polina/Desktop/Frontiers/R") 
tree<-readRDS("seqDSMZ.NJtree.root.RDS")
plot(tree)
tree$tip.label[29]<-"DSM100409_FOX_RLCS32"
tree$tip.label

ftree<-drop.tip(tree, c("DSM101518_DF13", "DSM100285_M_RLCS03","DSM100289_DF19_RLCS11", 
                        "DSM100331_C35_RLCS19" ,"DSM100284_C23_RLCS07" ,"DSM100288_DF17_RLCS29",
                        "DSM100292_DF37_RLCS25", "DSM100323_C31_RLCS28" ,"DSM100324_C28_RLCS17",
                        "DSM100327_A_RLCS21","DSM100328_DF35_RLCS31", "DSM100329_DF58_RLCS20",
                        "DSM100330_C21_RLCS26", "DSM100406_DF24_RLCS09", "DSM100410_J_RLCS24",
                        "DSM101518_DF13","DSM101519_C11_RLCS23"), root.edge = 0)
plot(ftree)

##Phylogenetic signal with package "phylosignal" 
library(phylosignal)
library(ape)
library(phylobase)

###sorting the trait data table to make it match the tree

##Read the trait table
library(readxl)
setwd("/Users/polina/Desktop/R.files/Liquid medium experiment/Phylogeny") #laptop
TRAITS<-read.csv("Phyl_ESmean.csv") #trait table
ftree # pruned tree (in this case including 32 strains)

TRAITS_A<-TRAITS[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33),]

row.names(TRAITS_A)<-c("DSM100400_C13_RLCS06","DSM100401_C29_RLCS22","DSM100326_C33_RLCS27","DSM100402_C34_RLCS15", "DSM100322_C40_RLCS04","DSM100403_C41_RLCS05","DSM100404_DF04_RLCS14",       "DSM100286_DF09_RLCS10","DSM100405_DF10_RLCS12","DSM100287_DF16_RLCS18","DSM100407_DF25_RLCS02","DSM100290_DF32_RLCS13", "DSM100291_DF36_RLCS30","DSM100408_DF42_RLCS16","DSM100293_DF56_RLCS01","DSM100325_DISC_RLCS08",  "DSM100409_FOX_RLCS32")

TRAITS_G<-TRAITS[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34),]
row.names(TRAITS_G)<-c("DSM100400_C13_RLCS06","DSM100401_C29_RLCS22","DSM100326_C33_RLCS27","DSM100402_C34_RLCS15", "DSM100322_C40_RLCS04","DSM100403_C41_RLCS05","DSM100404_DF04_RLCS14", "DSM100286_DF09_RLCS10","DSM100405_DF10_RLCS12","DSM100287_DF16_RLCS18","DSM100407_DF25_RLCS02","DSM100290_DF32_RLCS13", "DSM100291_DF36_RLCS30","DSM100408_DF42_RLCS16","DSM100293_DF56_RLCS01","DSM100325_DISC_RLCS08", "DSM100409_FOX_RLCS32")

idx <- sapply(ftree$tip.label,function(x){
  which(rownames(TRAITS_A)==x)
})
TRAITS_A <- TRAITS_A[idx,]

idx <- sapply(ftree$tip.label,function(x){
  which(rownames(TRAITS_G)==x)
})
TRAITS_G <- TRAITS_G[idx,]

###########---Biomass---###########
phyloTrait_B_A<-TRAITS_A[, c(4)] 
phyloTrait_B_G<-TRAITS_G[, c(4)] 
plot(ftree) #pruned tree, see above
treeP4D<-as(ftree, "phylo4")
phyTrait4d_B_A<-phylo4d(treeP4D, phyloTrait_B_A)
phyTrait4d_B_G<-phylo4d(treeP4D, phyloTrait_B_G)
barplot.phylo4d(phyTrait4d_B_A, center=FALSE, scale=FALSE, tree.type="phylogram") #plot phylogeny on trait values, for visualization
barplot.phylo4d(phyTrait4d_B_G, center=FALSE, scale=FALSE, tree.type="phylogram") #plot phylogeny on trait values, for visualization

Biomass_A<-phyloSignal(phyTrait4d_B_A) #compute phylogenetic signals and p-values (it can spit you out results for multiple traits at once)
Biomass_G<-phyloSignal(phyTrait4d_B_G)

###########---Cu accumulation---###########
phyloTrait_Cu_A<-TRAITS_A[, c(5)] 
phyloTrait_Cu_G<-TRAITS_G[, c(5)] 
plot(ftree) #pruned tree, see above
treeP4D<-as(ftree, "phylo4")
phyTrait4d_Cu_A<-phylo4d(treeP4D, phyloTrait_Cu_A)
phyTrait4d_Cu_G<-phylo4d(treeP4D, phyloTrait_Cu_G)

barplot.phylo4d(phyTrait4d_Cu_A, center=FALSE, scale=FALSE, tree.type="phylogram") #plot phylogeny on trait values, for visualization
barplot.phylo4d(phyTrait4d_Cu_G, center=FALSE, scale=FALSE, tree.type="phylogram") #plot phylogeny on trait values, for visualization

Cu_A<-phyloSignal(phyTrait4d_Cu_A) #compute phylogenetic signals and p-values (it can spit you out results for multiple traits at once)
Cu_G<-phyloSignal(phyTrait4d_Cu_G)

Phylogeny_P<-c(Biomass_A$pvalue[3], Biomass_G$pvalue[3], Cu_A$pvalue[3], Cu_G$pvalue[3])
P.bnf<-p.adjust(Phylogeny_P,method = "bonferroni")
P.bnf<-as.matrix(P.bnf)
colnames(P.bnf)<-c("P-value")
row.names(P.bnf)<-c("ES of Gradual treatment on Biomass","ES of Abrupt treatment on Biomass","ES of Gradual treatment on Cu uptake","ES of Abrupt treatment on Cu uptake")
P.bnf

Phylogeny_K<-c(Biomass_A$stat[3], Biomass_G$stat[3], Cu_A$stat[3], Cu_G$stat[3])

