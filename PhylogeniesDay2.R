####################
##Phylogenies day 2
####################
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


interval.length <- 1
times <- seq(0, 100, interval.length);

rate <- 0.1
root <- 0

normal.dev <- c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))

#sd calculation comes from the sigma squared t term from the brownian model

# take the cumulative sum and add in the root state to get trait values through time

traits <- cumsum(normal.dev);


plot(times, traits, type = "l", xlab = "time", ylab = "trait value", ylim = c(-10,10))


normal.dev2 <- 0

lastvalue <- C(1:1000)

for(i in 1:1000){
  normal.dev <- c(root, rnorm(n=(length(times)-1), mean = 0, sd = sqrt(rate * interval.length)))

  traits <- cumsum(normal.dev);
  lines(times, traits, type = "l", xlab = "time", ylab = "trait value", ylim = c(-10,10))
  lastvalue[i] <- (traits[100])
}

lastvalue
hist(lastvalue)
####################
#Example tree


phy <- "(((t1:0.15,t2:0.15):0.4,t3:0.55):0.5,(t4:0.25,t5:0.25):0.8);"
phy <- read.tree(text=phy)
plot(phy, label.offset=0.05)
edgelabels(c(0.5,0.4,0.15,0.15,0.55,0.8,0.25,0.25),adj=c(0.3,-0.3),frame="none",bg="",cex=0.8)
axisPhylo() # put up a scale bar

vcv(phy) # this is ape's vcv function

##Bray Curtis treats spp as separate in microbiome work
##or do a PCA cophenetic taking shared phylogeny into account, divide this 
#value by 2 to gt distance between spp


require(phytools)
set.seed(100)
d <- fastBM(tree=phy, a=root, sig2=rate)
d # take a look - we have five trait values

# visualize trait evolution on the tree
phenogram(phy,d,spread.labels=TRUE)


fancyTree(phy, x=d,type = "phenogram95" )

# Let's simulate another trait and compare to the first
d2 <- fastBM(tree=phy, a=root, sig2=rate)
phenogram(phy,d2,spread.labels=TRUE)

##########################################
##Probability of Observing 2 traits: use dmvnorm from mvtnorm package
#fit a Brownian motion model to trait data###################

#d is the vector of traits to observe
d

#use 2 loops, one for values of the rate and one for values of the root


# we will also need the variance covariance matrix
v <- vcv(phy)
v

require(mvtnorm)

# we can use dmvnorm to compute the likelihood of getting our data with our actual values

#dmvnorm parts: the traits vector, the mean vector

#dmvnorm(x=d, mean = rep(root, length(phy$tip.label)), sigma = v * rate[i], log = T)


rate_range <-seq(0, 0.2, by = 0.01)

  
root_range <-seq(0.1, 1.05, by = 0.02) 



probability <- expand.grid(rate_range, root_range)
probability

#Using Data Frame

for(i in 1:length(probability$Var1)){
  print(dmvnorm(x=d, mean = rep(probability$Var2[i], length(phy$tip.label)), sigma = v * probability$Var1[i], log = T))
}



##Useing Nested Lops 
ll <- c()

for (i in c(1:length(root_range))){ 
  for(j in c(1:length(rate_range))) {
     ll <- c(ll,dmvnorm(x=d, mean = rep(root_range[i], length(phy$tip.label)), 
                   sigma = (v * rate_range[j]), log = T))
  }
}


ll

###Can use fitContinuous instead of loops to fit evolutionary models ##########
require(geiger)

bm <- fitContinuous(phy=phy, dat=d, model="BM")
bm
##################################################

##Phylogenetic Independent Contrasts (PICs)

library(phytools)
set.seed(999)

## simulate a coalescent shaped tree
tree<-rcoal(n=100)
plotTree(tree,ftype="off")

## simulate uncorrelated Brownian evolution
x<-fastBM(tree, a=0, sig2=1)
y<-fastBM(tree, a=0, sig2=1)
plot(x,y,pch=20)
fit<-lm(y~x)
abline(fit)

summary(fit)
##but we know that we shouldn't be expecting a strong correlation, bc they
#were generated independently
#so why is there a strong correlation? Type 1 error

#Pic is phylogenetic independent contrasts 

ix<-pic(x, tree, scaled=TRUE)
iy<-pic(y, tree, scaled=TRUE)
ix
plot(ix,iy,pch=20)
fit<-lm(iy ~ ix - 1) ## we have to fit the model without an intercept term (this treats the contrasts as vectors)
abline(fit)

#relationship disappears when  using this method.
#the change in traits is consistent with evolutionary Brownian motion,
#not a true significant correlation with the 2 traits
#so Independent contrats help remove some Type 1 error rate

#show the sum of squared contrasts divided by n gives the
#maximum liklihood estimate of sigma^2

sum(ix^2)/length(ix)
mlx <- fitContinuous(phy = tree, dat= x, model= "BM")
mlx

sum(iy^2)/length(iy)
mly <- fitContinuous(phy = tree, dat= y, model= "BM")
mly
#####it's close

###Phylogenetic Generalized Leasst Squares#####
##OLS ordinary least squares

##show that OLS regression for Brownian MOtin traits(BM) simulataed on tree results elevate
#liklihood of Type 1 error

#alpha is the likelihood of type 1 ,want it to be below 0.05
#
#yi = alpha = beta(Xi) + error(i)
#
## simulate uncorrelated Brownian evolution
## simulate a coalescent shaped tree


sumtype1 <- 0


for (i in 1:200){
  x<-fastBM(tree, a=0, sig2=1)
  y<-fastBM(tree, a=0, sig2=1)
  #plot(x,y,pch=20)
  fit<-lm(y~x)
  #abline(fit)
  pvalue <- summary(fit)$coefficients[2,4] 
  
  if(summary(fit)$coefficients[2,4] < 0.05){
    sumtype1 <- sumtype1+ 1
    }
  
}

sumtype1
##returns lots of type 1 errors...every time p < 0.05 (significant) you know it's 
##wrong bc there trait values were randomly generated

tree<-rcoal(n=100)

fit

#show that PIC regression results correct the type 1 error

######################################
##Challenge##
##fix the WBC vs Body Size phylogeny

setwd("C:/Users/workshop/Desktop/IDEAS_PCM_Workshop-master/Intro_to_phylogenies_IDEAS")

dat <- read.csv("data/Cooper_2012.csv", as.is=T)
fritz_tree <- read.nexus("data/Fritz_2009.tre")[[1]]

# Remove species for which we don't have complete data
dat <- na.omit(dat)

# Match data to tree names
dat$Species_W.R05 <- gsub(x = dat$Species_W.R05 , replacement = "_", pattern = " ")

species.to.exclude <- fritz_tree$tip.label[!(fritz_tree$tip.label %in% dat$Species_W.R05)]


tree <- drop.tip(fritz_tree,species.to.exclude)


tree

# Order tree to make it nicer when plotting
tree <- ladderize(tree, right = FALSE)

# Name the rows of dat with the species codes remove obsolete columns
rownames(dat) <- dat$Species_W.R05
dat <- subset(dat, select=-c(Species_W.R05,Species_W.R93))

# Check that the order of the tree and data match
name.check(tree, dat)

tree <- multi2di(tree)

dat <- dat[tree$tip.label, ]


# Great! Time for analysis!

vcv(phy)

# Convert the covariance matrix to a correlation matrix
corrmat <- cov2cor(vcv(phy))
# Print the matrix, rounding the numbers to three decimals
round(corrmat,3)

corrmat <- vcv(phy,corr=TRUE)
round(corrmat,3)

require(nlme)
wbc.gls <- gls(WBC ~ log(AdultBodyMass_g), data=dat)
summary(wbc.gls)




vcv(phy)

# Convert the covariance matrix to a correlation matrix
corrmat <- cov2cor(vcv(phy))
# Print the matrix, rounding the numbers to three decimals
round(corrmat,3)

corrmat <- vcv(phy,corr=TRUE)
round(corrmat,3)

require(nlme)
wbc.gls <- gls(WBC ~ log(AdultBodyMass_g), data=dat)
summary(wbc.gls)


plot(tree, show.tip.label = FALSE)



resids <- as.numeric(residuals(wbc.gls))/2.646484
resids

tiplabels(cex = abs(resids), pch = 21, bg = ifelse(resids > 0, 'red', 'blue'))

wbc.gls$residuals

#add a legend
pts = c(0,1,2,3,4,5)   #temporary vector to feed into legend


legend("topleft", pch = 21, pt.cex = 2, c("Red = Positive Residual", "BLUE = Negative Residual"), title = "Tree WBC Phylogeny") 

axisPhylo()

###########################################################
##DAy 3##

# Calculate the correlation matrix from the tree
mat <- vcv(tree, corr=TRUE)
mat

# Create the correlation structure for gls
corr.struct <- corSymm(mat[lower.tri(mat)],fixed=TRUE)
corr.struct
# Run the pgls
wbc.pgls1 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=corr.struct)

# Note, summary(wbc.pgls1) returns the entire correlation matrix. In this example it is quire large with 213 species, so we will just look at the coefficient estimates.
summary(wbc.pgls1)$tTable


# Get the correlation structure
bm.corr <- corBrownian(phy=tree)
bm.corr


# PGLS
wbc.pgls1 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=bm.corr)
summary(wbc.pgls1)

###Challengs:Try the PGLS function on the Cooper2012 data

# Calculate the correlation matrix from the tree

# Get the correlation structure
bm.corr <- corBrownian(phy=tree)
bm.corr
wbc.pgls1 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=bm.corr)
summary(wbc.pgls1)

##Force through the origin (subtract 1)##


wbcPic <- pic(dat$WBC, tree, scaled=TRUE)
massPic <- pic(log(dat$AdultBodyMass_g), tree, scaled=TRUE)

regress <- lm(wbcPic ~ massPic - 1)


summary(regress)



##Why would you need to relax correlation structure?
##Brownian motion may not apply as much for different reasons

##see the lambda transformations (1, borwnina tree, 05. more time for evolution,
##0.1 max time for independednt spp evolution Karstar phylogenetic tree, no correlation structure)
##



names_wbc <- dat$WBC
names(names_wbc) <- rownames(dat)
names_wbc


lm_wbc <- fitContinuous(phy = tree, dat = names_wbc, model = c("lambda"))
bm_wbc <- fitContinuous(phy = tree, dat = names_wbc, model = c("BM"))

lm_wbc
bm_wbc

# the lambda model has a lower AIC value than the brownian model, so
# lamba is more useful

#try rescaling using lambda = 0.5 (the intermediate value)
rescale_tree = rescale(tree, model = c("lambda"), 0.5)

plot(rescale_tree,show.tip.label = FALSE)

#speciation model, discount branch length. This is what the punctuated equilibrium model uses


# Get the correlation structure using corPagel, similar to corBrownian
pagel.corr <- corPagel(0.3, phy=tree, fixed=FALSE)

# PGLS with corPagel
wbc.pgls2 <- gls(WBC ~ log(AdultBodyMass_g), data = dat, correlation=pagel.corr)
summary(wbc.pgls2)

lm_wbc_o <- gls((WBC ~ log(AdultBodyMass_g)+ Order), data = dat, correlation=pagel.corr)
summary(lm_wbc_o)

##AIC is lower when accounting for Order as a co-predictor. 
#but it is only a little bit lower
#so adding Order did not correct for non-independence

pagel.corr2 <- corPagel(0.3, phy=tree, fixed=FALSE)

wbc.pgls3 <- gls(PSR ~ log(AdultBodyMass_g), data = dat, correlation=pagel.corr2)
summary(wbc.pgls3)

#correlation is much weaker between PSR and body mass as compared to WBC and body mass.
#this evident due to the higher standard error, lower correlation

###Phylogenetic Uncertainty
##Load the 100 phylogenetic trees representing uncertainty

tree100 <- read.nexus("Faurby_2015_100trees.nex")
tree100



for (i in c(1:100)) {
  lambda <- fitContinuous(tree100[i,0], dat = dat, model = c("lambda"))
  lambda <- summary(lambda)$Parameter estimate(s)[5,1]
  print("lambda")
} 

print()


