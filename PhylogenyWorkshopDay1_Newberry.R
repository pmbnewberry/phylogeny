library(ape)
library(picante)
library(geiger)

ftree <- read.tree("sample_newick_tree.txt")

ftree2 <- read.nexus("sample_NEXUS_tree.txt ")

plot(ftree2)

str(ftree2)

plot(ftree)

str(ftree)
#drop.tip

plot(ftree)

ftree3 <- drop.tip(ftree, 4)
plot(ftree3)
#one way to drop a tip

plot(ftree)
ftree3 <- drop.tip(ftree, "A_boreas")
plot(ftree3)
#another way to drop tip

#drop the outgroup (there just to root the tree)
plot(ftree)
ftree3 <- drop.tip(ftree, c(21,22))
plot(ftree3)

ftree4 <- drop.tip(ftree2, c(13, 14, 15))
plot(ftree4)


ftree4 <- drop.tip(ftree2, c("R_aurora","R_cascadae", "R_luteiven"))
plot(ftree4)

plot(ftree2)

#extract clade
plot(ftree)
ftree3 <- extract.clade(ftree, 37)
plot(ftree3)

plot(ftree)
nodelabels() ############ tells you the node numbers

#pull out the P and H clade (position 29)
plot(ftree)
ftree3 <- extract.clade(ftree, 29)
plot(ftree3)




#lining up data and a tree
#four assumtions of many phylo methods:

#data and tree match perfectly (ie., all species in data occur 
#in tree, and the spelling of species labels in data matches 
#that of tip labels)

#no missing data (complete case)

#tree is fully dichotomous

#order of species in data matches that of the tip labels in 
#the tree


#loading sample tree
ctree <- read.tree("canid_tree.txt")
str(ctree)
plot(ctree)

#tree is not fully resolved, see "polytomies" when not sure of the relationship
ctree <- multi2di(ctree)  #how does this work?? adds no interior branches, still zero length branches 
#can also add a tiny bit of tolerance
str(ctree)

#loading sample data
cdata <- read.csv("Canid_traits.csv")
str(cdata)

###using treedata to get body mass data lined up with tree###

#extracting and cleaning body mass data mass
cmass <- cdata$AdultBodyMass_g

#assigning sprecies names to the mass data
names(cmass) <- cdata$Binomial


#eliminating missing data
hist(cmass)
range(cmass)
cmass <- cmass[cmass > 0]

#log transforming (optional)
hist(cmass)
hist(log10(cmass))
cmass <- log10(cmass)  ##Log transform the data after checking the graph
hist(cmass)
shapiro.test(cmass)  ##p value 0.8624, sig. similar to normal distribution

#using treedata, geuiger allows lining up the tree and data
library(geiger)

#treedata command (tree, data to line up with the branches, pullout non matching data)

CMout <- treedata(ctree, cmass)
str(CMout)
plot(CMout[[1]])  ##hard brackets indexes elements of a list
hist(CMout[[2]])

#final cleanup
CMtree <- CMout[[1]]
plot(CMtree)

#returns dataframe
CMass2 <- CMout[[2]]
head(CMass2) #returns a vector, we want a dataframe

#returns vector
CMass2 <- CMout[[2]][,1] 

head(CMass2)
str(CMtree)

#note that data is not in same order as tip labels
#to fix
CMass2 <- CMass2[CMtree$tip.label]

head(CMass2)
str(CMtree)


####################################
#Goal: assing eye opening data to phylogeny
cdata <- read.csv("Canid_traits.csv")
str(cdata)

#extracting and cleaning eye opening data
ceyes <- cdata$AgeatEyeOpening_d

#assigning species names to the eye opening data
names(ceyes) <- cdata$Binomial



#eliminating missing data
hist(ceyes)
range(ceyes)
ceyes <- ceyes[ceyes > 0]

#treedata command (tree, data to line up with the branches, pullout non matching data)

CEout <- treedata(ctree, ceyes)
str(CEout)
plot(CEout[[1]])  ##hard brackets indexes elements of a list
hist(CEout[[2]])

#final cleanup
CEtree <- CEout[[1]]
plot(CEtree)

#returns dataframe
CEyes2 <- CEout[[2]]
head(CEyes2) #returns a vector, we want a dataframe

#returns vector
CEyes2 <- CEout[[2]][,1] 

head(CEyes2)
str(CEtree)


#note that data is not in same order as tip labels
#to fix
CEyes2 <- CEyes2[CEtree$tip.label]

head(CEyes2)
str(CEtree)
##################################################

#testing for phylogenetic signal 
#if no phylogenetic signal, easier to just use normal statistics
#so test for it, to see whether to use these tools

#look back at Mass

library(picante)

phylosignal(CMass2, CMtree)

#returns:          K              PIC.variance.obs      PIC.variance.rnd.mean PIC.variance.P PIC.variance.Z
#              1 1.133742       0.02153078            0.06593864          0.001      -3.752943

#now run eye opeing data for phylo symbol

phylosignal(CEyes2, CEtree)

#returns   K                 PIC.variance.obs      PIC.variance.rnd.mean PIC.variance.P PIC.variance.Z
#       1 0.5829063         2.109657              2.492379          0.229     -0.7524992

#weaker signal than with mass
#neither show much signal 
######################################
#Tree plotting basics

#adding a scale bar

#edata2 is CEyes2
#etree in CEtree


plot(ftree)
add.scale.bar()

plot(ftree)
axisPhylo()  #plot across whole tree

################################# Making a plot ########################################################

#continuous character reconstruction

#creating a color pallete
mypallete <- c("red","orange","yellow","green","blue","purple")

#creating bins for species traits
divider <- (max(CMass2+0.0001)-min(CMass2))/6               #choose a binning scheme (divider returns 0.2511593)
index <- floor((CMass2-min(CMass2))/divider)                #assign index values to match the bins

#ACE
CMtree$edge.length <- CMtree$edge.length+0.001              #add some length to the zero length branches
recon <- ace(CMass2, CMtree, CI = F)
recon.values <- recon$ace
recon.index <- floor((recon.values-min(CMass2))/divider)   #converting interior values to binning scheme

#creating a plot
plot(CMtree, label.offset = 0.1, edge.width = 2)
nodelabels(pch = 21, cex = 2, bg = mypallete[recon.index+1])
tiplabels(pch = 21, cex = 2, bg = mypallete[index+1])

#add a legend
pts = c(0,1,2,3,4,5)   #temporary vector to feed into legend
legend("topright", pch = 21, pt.cex = 2, c("0.00-0.25","0.25-0.50","0.50-0.75","0.75-1.00","1.00-1.26","1.26-1.51"), title = "Log 10 Mass", pt.bg = mypallete[pts+1])

##################################

#Discrete reconstruction

#getting a continuous character

AC <- cdata$ActivityCycle
names(AC) <- cdata$Binomial
AC <- AC[AC > 0]         #pull out -999's

#aligning tree and data
out2 <- treedata(ctree, AC)

CAtree <- out2[[1]]   
CAtree <- multi2di(CAtree)  #resolve polytomies
CAtree$edge.length <- CAtree$edge.length+0.001 #add length to branches

CAC <- out2[[2]][,1]
CAC <- CAC[CAtree$tip.label]

#peforming reconstruction
anC <- ace(CAC, CAtree, type = "d")  #respective probabilities belonging to a particular state
states <- anC$lik.anc

#making the plot
plot(CAtree, edge.width = 2, label.offset = 1)
co <- c("white", "gray", "black")
tiplabels(pch = 22, bg = co[as.numeric(CAC)], cex = 2, adj = 1)
nodelabels(pie = states, piecol = c("white", "gray", "black"), cex = 0.5)
axisPhylo()

nodelabels()
tiplabels()

#adding a legend
pts = c(1,2,3)
legend("topright", pch = 22, pt.cex = 2, c("nocturnal","crepuscular","diurnal"), title = "Activity Cycle", pt.bg = co[pts])

###########conduct a reconstruction as above with any trait in canid database. Export the image. 
#
#COntinuous trait to reconstruct:AdultHeadBodyLen_mm

#AHBL <- cdata$AdultHeadBodyLen_mm
#names(AHBL) <- cdata$Binomial
#AHBL <- AHBL[AHBL > 0]

#aligning tree and data
#out3 <- treedata(ctree, AHBL)



#AHtree <- out3[[1]]   
#AHtree <- multi2di(AHtree)  #resolve polytomies
#AHtree$edge.length <- AHtree$edge.length+0.001 #add length to branches with 0.0 length

#AHC <- out3[[2]][,1]
#AHC <- AHC[AHtree$tip.label]

#peforming reconstruction
#anC <- ace(AHC, AHtree, type = "d")  #respective probabilities belonging to a particular state
#states <- anC$lik.anc

####################Attempt 2
#continuous character reconstruction

AHBL <- cdata$AdultHeadBodyLen_mm
#AHBL <- as.numeric(cdata$AdultHeadBodyLen_mm)
names(AHBL) <- cdata$Binomial
AHBL <- AHBL[AHBL > 0]

#aligning tree and data
out3 <- treedata(ctree, AHBL)

AHtree <- out3[[1]]   
AHtree <- multi2di(AHtree)  #resolve polytomies
AHtree$edge.length <- AHtree$edge.length+0.001 #add length to branches with 0.0 length

#AHC <- out3[[2]][,1]
#AHC <- AHC[AHtree$tip.label]

#peforming reconstruction
anC <- ace(AHC, AHtree, type = "d")  #respective probabilities belonging to a particular state
states <- anC$lik.anc

#creating a color pallete
mypallete <- c("red","orange","yellow","green","blue","purple")


#creating bins for species traits
              #choose a binning scheme (divider returns 0.2511593)
index <- floor((AHBL-min(AHBL))/divider)                #assign index values to match the bins

#ACE
AHtree$edge.length <- AHtree$edge.length+0.001              #add some length to the zero length branches
recon <- ace(anC, AHtree, CI = F)
recon.values <- recon$ace
recon.index <- floor((recon.values-min(anC))/divider)   #converting interior values to binning scheme

#creating a plot
plot(AHtree, label.offset = 0.1, edge.width = 2)
nodelabels(pch = 21, cex = 2, bg = mypallete[recon.index+1])
tiplabels(pch = 21, cex = 2, bg = mypallete[index+1])

#add a legend
pts = c(0,1,2,3,4,5)   #temporary vector to feed into legend
legend("topright", pch = 21, pt.cex = 2, c("0.00-145.2 mm","145.2-290.5 mm","290.5-435.7 mm","435.7-581.0 mm","581.0-726.3 mm","726.3-871.5 mm"), title = "Adult Head Body Length", pt.bg = mypallete[pts+1])

axisPhylo()


##############Add all Carnivora#########








