#############Phylogenetic diversity##########
##
##
library(ape)
library(picante)
library(geiger)
library(ggplot2)
library(dplyr)

library(nlme)
library(ape)
library(picante)
library(geiger)
library(phytools)
library(mvtnorm)
library(brms)
library(phytools)

phy <- "(((t1:0.15,t2:0.15):0.4,t3:0.55):0.5,(t4:0.25,t5:0.25):0.8);"
phy <- read.tree(text=phy)
plot(phy, label.offset=0.05)
edgelabels(c(0.5,0.4,0.15,0.15,0.55,0.8,0.25,0.25),adj=c(0.3,-0.3),frame="none",bg="",cex=0.8)
axisPhylo() # put up a scale bar

data(phylocom)

phy$edge.length


pd_tree <- sum(phy$edge.length)
pd_tree

# GMPD 2.0
gmpd <- read.csv("data/GMPD_datafiles/GMPD_main.csv", as.is=T)
# Removing parasites not reported to species
Sys.setlocale('LC_ALL','C') 
gmpd <- gmpd[grep("sp[.]",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("ABOLISHED",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("no binomial name",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("not identified to genus",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("SPLITTED in ICTV",gmpd$ParasiteCorrectedName, invert=TRUE),]
gmpd <- gmpd[grep("Diphyllobothrium sp",gmpd$ParasiteCorrectedName, invert=TRUE),]

# Removing hosts with no binomial name reported
gmpd <- gmpd[grep("no binomial name",gmpd$HostCorrectedName, invert=TRUE),]

# Tree
fritz_tree <- read.nexus("data/Fritz_2009.tre")[[1]]

# Binomial cleaning & matching
gmpd$HostCorrectedName <- gsub(" ", "_", gmpd$HostCorrectedName)
species.to.exclude <- fritz_tree$tip.label[!(fritz_tree$tip.label %in% 
                                               gmpd$HostCorrectedName)]
# Subsetting tree
gmpd_tree <- drop.tip(fritz_tree,species.to.exclude)

# Community matrix
com <- table(gmpd$HostCorrectedName, gmpd$ParasiteCorrectedName)

##Do this with canids only
canid <- subset(gmpd, gmpd$HostFamily == 'Canidae') 
canid
##########

# Binomial cleaning & matching
canid$HostCorrectedName <- gsub(" ", "_", canid$HostCorrectedName)
species.to.exclude <- fritz_tree$tip.label[!(fritz_tree$tip.label %in% 
                                               canid$HostCorrectedName)]
# Subsetting tree
canid_tree <- drop.tip(fritz_tree,species.to.exclude)

# Community matrix
canidcom <- table(canid$HostCorrectedName, canid$ParasiteCorrectedName)
canidcom



#canidcom$presence <- ifelse(canidcom$Freq > 0, 1, 0)


canidcom[canidcom > 1] <- 1



# pd expects rows to be "sites", and columns to be "species"
# in our case, since we are interested in host PD, our parasites are "sites"
# and we need to transpose our community matrix.

com <- canidcom

com <- t(com)

# Now calculate PD
com.pd <- pd(com, canid_tree, include.root=TRUE)
head(com.pd)

# Compare PD and species richness
plot(com.pd$PD ~ com.pd$SR, xlab = "Species richness", ylab = "Faith's PD", pch=19)











###################################################
###now subset carnivora##
#



##Do this with canids only
canid2 <- subset(gmpd, gmpd$HostFamily == 'Canidae') 
canid2
str(canid2)

paracom <- subset(canid2, canid2$HostOrder == 'Carnivora' )
paracom
