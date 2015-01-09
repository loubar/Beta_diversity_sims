# Louise Barwell & Nick Isaac
# 7 October 2011 - 24 July 2014

### Measuring pairwise beta-diversity with abundance data

## functions used in simulating assemblage pairs with different compositions and abundance structures

###############################################################

library(vegan)
library(vegetarian)
library(fossil)
library(untb)

##############################################################################################################################



# Simulate a Fisher log series rank abundance distribution.
# This will be the focal assemblage that is manipulated in the following simulations 
# The Fisher logseries describes a species abundance distribution.  Converting this to a rank abundance 
# distribution generates uncertainty as the exact abundances of species in
# in each abundance class is unknown.  There is considerable variation in the assemblages
# generated for the same number of species and individuals

# to ensure that the assemblage always has exactly N individuals and S species, repeat the fisher.ecosystem function until the
# assemblage returned falls within these limits N=10000+-1%, S=100:

require(untb) # contains function fisher.ecosystem for generating synthetic Fisherian assemblages 

FisherRAD<-function(N, S){
  repeat {
    z<-fisher.ecosystem(N=N, S=S, nmax=N) # from package 'untb'
    if (length(z)==S & sum(z)>9990 & sum(z)<10010) break # N is fixed between 9990 and 10010 to avoid confounding effect of sample size and effect of alpha diversity
  }
  return(as.numeric(z))
}
#  S is fixed exactly
#  fixing N exactly takes much longer than setting tolerance limits of 10000+-10 

#########################################################################################################################################
##################################################################################################################################

# Function to simulate undersampling of both assemblages
# Used to evaluate desirable property D1: unbiased by sample size 

# Comm1 is the focal assemblage
# N is the sample size 
# p is the proportion of species to turnover
sampling_effect<-function(Comm1, N, p){ #  a function to randomly subsample from BOTH assemblages
  turn<-random_composition_change(Comm1, p=p) # generate a second community with turnover p
  return(rbind(tabulate(sample(rep(1:length(turn[1,]), turn[1,]), N), nbin=length(turn[1,])), tabulate(sample(rep(1:length(turn[2,]), turn[2,]), N), nbin=length(turn[2,])))) # sample Comms 1 and 2
}

###################################################################################################################################

# A function to simulate unequal sample sizes
# Used to evaluate desirable property D2: unbiased by sample size

unequal_sampling<-function(Comm1, N, p){ # a function to randomly subsample from ONE assemblage
  turn<-random_composition_change(Comm1, p=p) # generate a second community with turnover p
  return(rbind(turn[1,], tabulate(sample(rep(1:length(turn[2,]), turn[2,]), N), nbin=length(turn[2,])))) # undersample Comm2 only
}
###################################################################################################################################
# A function to simulate assemblage pairs with a different level of alpha-diversity in each pair
# note that alpha diversity within pairs is equal
# Used to evaluate desirable property D3: unbiased by alpha-diversity 

#  N is the expected number of individuals  
#  S is the expected number of species (must be divisible by 10 so that the number of species to turover over is an integer).  
gen_alphas<-function(N, S){ # a function to generate a Fisherian ecosystem with N individuals and S species
  repeat {
    z<-fisher.ecosystem(N=N, S=S, nmax=N) # from package 'untb'
    if (length(z)==S & sum(z)>(N*0.99) & sum(z)<(N*1.01)) break # N is fixed between 9990 and 10010 to avoid confounding effect of sample size and effect of alpha diversity
  }
  #  S is fixed
  #  fixing N exactly takes much longer than setting tolerance limits of 10000+-10 
  return(as.numeric(z))  # an assemblage with N = 10000+-10 and S = S: together these values determine Fisher's alpha
  #return(hist(rep(1:length(x), x))) # can use this line to plot the species abundance distribution
}
###################################################################################################################################
# a function to simulate an environmental gradient generating directional species turnover across three assemblages 
# Used to evaluate desirable property D4: additivity under directional turnover

# Comm1 is the starting community
# gradient specifies the how many times more likely are species shared between Comm 1 and Comm2 to be turned over  
# this makes the turnover directional, rather than random (e.g. Comm1 and Comm3 will be more different than Comm1 and Comm2, as they lie further apart along the environmental gradient)
# p is the proportion turnover between each assemblage pair
# I assume a constant rate of proportion of turnover (p) in space or along a environmental gradient but....
# .... I could specify p1 and p2 in the function to allow the rate of turnover to vary.  
# together p and gradient determine the the spatial or environmental gradient
# 

add_comp_change <- function(Comm1, p, gradient){ # a function to generate 3 assemblages with directional turnover
  #takes a vector of community abundances
  #randomly changes the identity of a set proportion
  Comm<-Comm1
  #Community 2 is the same as comm1 but with the identities of certain species changed
  Comm2 <- Comm
  #Randomly select a proportion, p, of species (by rank).  
  to.change <- sample(1:length(Comm), p * length(Comm))
  #Reassign the identity of these species. (e.g. give new rank id not present in community 1, a “new” species that is absent in community 1) 
  #set the abundance of these selected species to zero
  Comm2[to.change] <- 0
  #append the abundances of these 'new' species
  Comm2 <- c(Comm2, Comm[to.change])
  a<-which(Comm1>0) # Comm1 spp.
  b<-which(Comm2>0) # Comm2 spp.
  x<-setdiff(b,a) # unique to Comm2 
  y<-setdiff(a,b) # unique to Comm1
  z<-intersect(b,a) # shared spp. between comm 1 and comm2.  these should be disproportionately likely to be lost in Comm 3 along an environmental gradient
  probs1<-c(rep(gradient, length(z)), rep(1, length(x))) # species still present from Community 1 are 10 times as likely to be lost than species only present in Community 2 - a gradient of turnover
  probs<-probs1/(sum(probs1))
  Comm3 <- Comm2
  to.change2<- sample(b, p*length(Comm), prob=probs)
  Comm3[to.change2] <- 0
  Comm3<-c(Comm3, Comm2[to.change2])
  #these 'new' species now have zero abundance in community 1
  Comm <- c(Comm, rep(0, length(to.change)+length(to.change2)))
  Comm2<-c(Comm2, rep(0, length(to.change2)))
  return(rbind(Comm, Comm2, Comm3))
  
} 
####################################################################################################################################
# A function to simulate turnover when samples are orthogonal in space
# We expect similarity (the complement of beta-diversity) to be probabilistic
# Used to evaluate desirable property D5:  multiplicative similarity under orthogonal sampling


mult_comp_change <- function(Comm1, p){ # p is the proportion of species to conserve between assemblage pairs 
  # we expect similarity, not beta diversity, to be multiplicative if turnover is random - the probability,p, of observing a species in the third quadrat is (prob in assemblage 1 * prob in assemblage 2)  
  Comm<-Comm1
  Comm2 <- Comm1
  to.conserve <- sample(1:length(Comm), p * length(Comm)) # randomly choose p species to conserve
  if(length(to.conserve)==0){
    to.turnover<-1:length(Comm)
  } else 
    to.turnover<-(1:length(Comm))[-to.conserve]
  Comm2[to.turnover] <- 0
  Comm2 <- c(Comm2, Comm[to.turnover])
  a<-which(Comm>0) # all Comm1 spp.
  b<-which(Comm2>0) # all Comm2 spp.
  x<-setdiff(b,a) # unique to Comm2 
  y<-setdiff(a,b) # unique to Comm1
  z<-intersect(b,a) # shared between assemblages 1 and 2 
  Comm3 <- Comm2
  to.conserve2<- sample(b, p*length(b))
  if(length(to.conserve2)==0){
    to.turnover2<-b
  } else 
    to.turnover2<-b[-to.conserve2]
  Comm3[to.turnover2] <- 0
  Comm3<-c(Comm3, Comm2[to.turnover2])
  Comm <- c(Comm, rep(0, length(Comm3)-length(Comm)))
  Comm2<-c(Comm2, rep(0, length(Comm3)-length(Comm2)))
  return(rbind(Comm, Comm2, Comm3))  
} 
################################################################################################################################

# functions to generate species turnover between assemblage pairs 
# Used to evaluate desirable property D6: strict monotonic increase in the value of beta with species turnover

random_composition_change <- function(Comm1, p){ 
	#takes a vector of community abundances
	#randomly changes the identity of a set proportion, p
	#Community 2 is the same as comm1 but with the identities of randomly chosen species changed
	Comm2 <- Comm1
  #Randomly select a proportion, p, of species (by rank).  
	to.change <- sample(1:length(Comm1), p * length(Comm1))		
	#Reassign the identity of these species. (e.g. give new rank id not present in community 1, a “new” species that is absent in community 1) 
	#set the abundance of these selected species to zero
	Comm2[to.change] <- 0
	#append the abundances of these 'new' species
	Comm2 <- c(Comm2, Comm1[to.change])
  #these 'new' species have zero abundance in community 1
	Comm1 <- c(Comm1, rep(0, length(to.change)))
	return(rbind(Comm1, Comm2))
	}

# as random_composition_change function above, but simulating both turnover and species richness differences simultaneously 
# species richness differences are also called nestedness (all species present in the species poor assemblage are a subset of those in the species rich assemblage)
# personality property P1: broad-sense or narrow sense?  
# personality property P3: the relative sensitivity to species turnover and species richness differences

# functions to generate assemblages with both nestedness and turnover
nested<-function(S, p, Comm1){ # Comm2 is nested in Comm1
  Comm1a<-Comm1
  choose.nested<-sort(sample(1:length(Comm1a), S)) # randomly select which species are present in the nested assemblage
  Comm2a<-rep(0, length(Comm1a))
  Comm2a[choose.nested]<-Comm1a[choose.nested]
  # scale the abundances back up to 10000 (same as reference assemblage, Comm1): this avoids confounding effect of nestedness with sample size
  comm2b<-round((Comm2a/sum(Comm2a))*sum(Comm1))
  x<-sum(Comm1)-sum(comm2b)
  comm2b[which(comm2b>0)[1]]<-comm2b[which(comm2b>0)[1]]+x # add any missing individuals to the most abundance species 
  # now introduce some turnover, too
  to.turnover<-sample(which(comm2b>0), size=round(p*length(which(comm2b>0)), 0))
  comm2b<-c(comm2b, comm2b[to.turnover])
  comm2b[to.turnover]<-0
  Com1<-c(Comm1a, rep(0, length(comm2b)-length(Comm1a)))
  return(rbind(Com1, comm2b, deparse.level=0)) # Comm2 is nested within Comm1
}

################################################################################################################################

############################################################################################################################## 

# a function to simulate decoupling of species ranks between assemblage pairs 
# used to evaluate 
# desirable property D7: strict monotonic increase with decoupling of species ranks
# desirable property D9: value of beta under extreme decouplign of species ranks < beta when species turnover is complete
# personality property P4: Relative sensitivity to decoupling of species ranks and species turnover components of beta   

# a function to generate assemblage pairs with a given partial correlation between species ranks
par_cor<-function(Comm1, r, p){ # r is the partial correlation between species rank abundance in assemblage pairs (e.g. the degree to which the ranks in assemblage 2 are determined by the ranks in assemblage 1, versus how much they are determined at random)
  rran<-sample(1:length(Comm1))
  oran<-as.numeric(1:length(Comm1))
  new<-order((r*oran)+((1-(abs(r)))*rran)) # the rank order is partially determined by the rank in the first assemblage and partially at random.
  Com2<-Comm1[new]
  to.change <- sample(1:length(Comm1), p * length(Comm1))  
  Com2[to.change] <- 0
  #append the abundances of these 'new' species
  Comm2 <- c(Com2, Comm1[to.change])
  Comm1 <- c(Comm1, rep(0, length(to.change)))
  return(rbind(Comm1, Comm2))
}
###################################################################################################################################

###################################################################################################################################

# a function to simulate differences in evenness among assemblage pairs at different levels of abundance 
# used to evaluate desirable properties: 
# D8 : strict monotonic increase with differences in evenness among assemblage pairs
# D10: beta under extreme differences in evenness < beta when species turnover is complete
# P5: relative sensitivity to evenness differences and species turnover components of beta

evenness<-function(Comm1, even, turnover){ # Comm1 is the Fisher log series assemblage, even is the chosen starting assemblage (from perfectly even to extremely uneven), turnover is the propertion of species to turnover
  ref<-rep(100, 100)
  redist<-c(sum(Comm1)-sum(rep(1, length(Comm1)-1)), rep(rep(1, length(Comm1)-1)))
  if (all(even==ref) | all(even==c(sum(Comm1)-sum(rep(1, length(Comm1)-1)), rep(rep(1, length(Comm1)-1))))){ # perfectly even assemblage and perectly uneven assemblage
    comm2<-even
  } else { # else manipulate the probabilities of being assigned to each species
    redist2<-rep(1, length(Comm1))
    comm2<-sort(tabulate(sample(1:100, redist[1]-1, replace=T, prob=even), nbin=100)+redist2, decreasing=T)
  }
  to.change <- sample(1:length(comm2), turnover * length(comm2))
  comm2a <- c(comm2, comm2[to.change])
  comm2a[to.change] <- 0
  refa<-c(ref, rep(0, length(comm2a)-length(ref)))
  return(rbind(refa, comm2a, deparse.level=0))
}


########################################################################################################################################

# function to generate species replication(Legendre and de Caceres 2013 - P7)
 

Srep_inv<-function(Comm, sploss, turnover){
  reps <- vector("list", 10)
  z <- nestedness(Comms=Comm, sploss, turnover, rev=FALSE)
  reps[[1]]<-z
  x<-z
  for (i in 1:9){
    x<-cbind(x,z)
    reps[[i+1]]<-x
  }
  return(reps)
}


################################################################################################################################


# a function to generate a pair of assemblages with different alpha diversities and with different levels of turnover
# Alpha 1 is the starting assemblage and Alpah 2 is the second assemblage 
# The function is given pairs of assemblages with increasing differences in  alpha diversity
# the assemblages of different alpha -diversity are created using the function gen_alphas above
# used to evaluate the personality trait P2: sensitivity to differences in alpha-diversity
# this is the abundance-based equivalent of sesnitivity to differences in species richness
# and is intended to distinguish between 'narrow-sense' and 'broad-sense' abundance metrics (sensu Koleff et al. 2003)

Diff_alpha<-function(Alpha1, Alpha2, turnover){
  Comm1<-Alpha1
  Comm2<-c(Alpha2, rep(0, length(Comm1)-length(Alpha2)))
  to_turn<-sample(1:length(Comm2[Comm2>0]), turnover*length(Alpha2))
  Comm2[to_turn]<-0
  new_sp<-Alpha2[to_turn]
  Comm2<-c(Comm2, new_sp)
  Comm1<-c(Comm1, rep(0, length(Comm2)-length(Comm1)))
  return(rbind(Comm1, Comm2))
}
###############################################################################################################################

# function to generate nestedness (species loss and different levels of turnover)
# We use this to test whether metrics are broad sense (sensitive to species turnover and differences in alpha-diversity) or narrow sense (only sensitive to species turnover)
 
nestedness<-function(Comms, sploss, turnover, rev){
  Comm2<-Comms
  to_lose<-sample(1:length(Comms), sploss)
  Comm2[to_lose]<-0 # Comm2 is a nested subset of the species in Comm1
  # but this time we do not scale the abundances back up.  Both the species and their abundances are nested
  replace<-sample(which(Comms>0 & Comm2>0), turnover*length(which(Comm1>0 & Comm2>0))) # Of the shared species, a proportion, turnover, are replaced by species of a different identity in Comm2 
  Comm2[replace]<-0
  Com2<-c(Comm2, Comm1[replace])
  Com1<-c(Comms, rep(0, length(replace)))
  if(rev==TRUE){
    return(rbind(Com2, Com1))
  }
  else 
    return(rbind(Com1, Com2))
}

########################################################################################################################################
# a function to generate turnover, while assuming a positive occupancy-abundance relationship (ONR)
# a positive ONR is a near-ubiquitous macroecological pattern
# it reflects that species with low abundance are also typically range-limited
# I simulate this pattern here by allowing rare species to be more likely to be turned over
# the probability of being turned over is inversely proportional to abundance of a species
# as opposed to random turnover for all species, as specified in all previus simulations
# I am interested in the difference this makes to the value of beta-diversity for the different metrics
# The prediction is that the value of beta-diversity will be much lower for the same level of turnover under a positive ONR,
# because mainly species will low abundances will have been turned (therefore carrying less weight with abundance-based metrics) 
# Metrics influences predominantly by turnover in dominant species are expected to show significant decreases in the value of beta for a given
# level of turnover.
# This might be problematic for detecting beta-diversity patterns under a positive ONR.

# used to evaluate personality property P6: effect of a positive ONR on the value of beta-diversity

# Comm1 is the focal assemblage
# p is the proportion of turnover

rare_composition_change <- function(Comm1, p){ 
  #takes a vector of community abundances
  #randomly changes the identity of a set proportion
  Comm<-Comm1
  #Community 2 is the same as comm1 but with the identities of p species changed
  Comm2 <- Comm
  probs<- (1/(Comm)^0.65) # prob of a species being turned over is proportional to the inverse of its relative abundance.  Equates to Singletons being around 116x more likely to be turned over than the most common species
  # this translates to an occupancy-abundance relationship with a slope of ~1.54 (on log-log axes with occupancy on the x axis and abundance on the y)
  # this is consistent with empirical interspecific occupancy-abundance relationships in the literature.
  #select a proportion, p, of species (by rank).  
  to.change <- sample(1:length(Comm), p * length(Comm), prob=probs)   
  #Reassign the identity of these species. (e.g. give new rank id not present in community 1, a “new” species that is absent in community 1) 
  #set the abundance of these selected species to zero
  Comm2[to.change] <- 0
  #append the abundances of these 'new' species
  Comm2 <- c(Comm2, Comm[to.change])
  #these 'new' species have zero abundance in community 1
  Comm <- c(Comm, rep(0, length(to.change)))
  return(rbind(Comm, Comm2))
}
####################################################################################################################################


##################################################################################################################################### 

# a function to turnover a single species in a Fisher log series species abundance distribution
# and examine the way that metrics use abundance information (how much is turnover weighted by the abundance?)

# weighting of turnover by abundance 

# Comm1 is  the focal assemblage
# sp.rank is the rank of the single species to be turned over

abd_sens<-function(Comm1, sp.rank){ # a function to generate an assemblage pair with a single species turned over
  Comm2<-Comm1
  Comm2[sp.rank]<-0
  Comm2<-c(Comm2, Comm1[sp.rank])
  Com1<-c(Comm1, 0)
  return(rbind(Com1, Comm2))
}



##################################################################################################################################### 

# a function to calculate the value of the beta-diversity metrics not implemented in vegdist
# see the functions to calculate each of the metrics 
new_metrics<-function(Comms){ 
	Hor<-Horn2(Comms) # Horn
	Renk<-Renkonen(Comms)
	AvEuc<-AvEuclidean(Comms)
	ClassicJac<-as.vector(1-(betadiver(Comms, method = 10)))
	ClassicSor<-as.vector(1-(betadiver(Comms, method = 11)))
	sim<-as.vector(betadiver(Comms, method="sim"))
	ChaoSor<-chaojs(Comms, method="sorensen")
	ChaoJac<-chaojs(Comms, method="jaccard")
	NESS50<-NESSm50(Comms)
	JostShan<-JostShannon(Comms)
	JostSimp<-JostSimpson(Comms)
	LandeShan<-LandeShannon(Comms)
	LandeSimp<-LandeSimpson(Comms)
	Part <- part(Comms)
  BaselgaBC_turn<-Part$BaselgaBC_turn 
  BaselgaBC_nest<-Part$BaselgaBC_nest
  BaselgaR_turn<-Part$BaselgaR_turn
  BaselgaR_nest<-Part$BaselgaR_nest
  PodaniBC_turn<-Part$PodaniBC_turn
  PodaniBC_nest<-Part$PodaniBC_nest
  PodaniR_turn<-Part$PodaniR_turn
  PodaniR_nest<-Part$PodaniR_nest
	return(c(Hor, Renk, AvEuc, ClassicJac, ClassicSor, sim, ChaoSor, ChaoJac, NESS50, JostShan, JostSimp, LandeShan, LandeSimp, BaselgaBC_turn, BaselgaBC_nest, BaselgaR_turn, BaselgaR_nest, PodaniBC_turn, PodaniBC_nest, PodaniR_turn, PodaniR_nest))
} # returns a vector with the value of beta for each of the new indices	

##############################################################################################################################

# a function to calculate the value of the all 25 beta-diversity metrics for a focal assemblage pair 
beta_metrics_all <- function(comms){ 
	#takes a pair of community abundance vectors and calculates multiple beta diversity indices
	method <- list('jaccard','canberra','bray','gower','kulczynski','morisita','horn','euclidean','manhattan','altGower','binomial', 'cao')
	y<-sapply(method, function(x) vegdist(comms, method=x), USE.NAMES=T)
	names(y) <- c('Ruzicka','Canberra','Bray-Curtis','Gower','Kulczynski','Morisita','Morisita-Horn','Euclidean','Manhattan','alt. Gower','Binomial', 'CYd')
	z<-new_metrics(comms)
	newmeth<-list("Horn", "Renkonen", "Av. Euclidean", "Classic Jaccard", "Classic Sorensen", "sim", "Chao Sorensen", "Chao Jaccard",  "NESS", "Jost Shannon", "Jost Simpson", "Lande Shannon", "Lande Simpson", "Baselga B-C turn", "Baselga B-C nest", "Baselga Ruzicka turn", "Baselga Ruzicka nest", "Podani B-C turn", "Podani B-C nest", "Podani Ruzicka turn", "Podani Ruzicka nest")
	names(z)<-newmeth
	return(c(y,z))
	} # returns a vector containing the value of beta for all 33 metrics in the analysis
###################################################################################################################################

# functions to calculate the metrics not implemented in vegdist in vegan package

# Renkonen

Renkonen<-function(Comms){
  a<-as.data.frame(rbind(Comms[1,]/(sum(Comms[1,])),Comms[2,]/(sum(Comms[2,]))))
  b<-sapply(1:length(a), function(i) min(a[,i]))
  return(1-sum(b))
}

# Horn 
Horn2<-function(Comms){
  k<-log(colSums(Comms), 10)
  k[k==-Inf]<-0
  a<-sum((colSums(Comms))*k)
  b<-sum(na.omit(Comms[1,]*log(Comms[1,],10)))
  d<-sum(na.omit(Comms[2,]*log(Comms[2,],10)))
  e<-(sum(Comms[1,])+sum(Comms[2,]))*(log(sum(Comms[1,])+sum(Comms[2,]),10))
  f<-sum(Comms[1,])*(log(sum(Comms[1,]),10))
  g<-sum(Comms[2,])*(log(sum(Comms[2,]),10))
  return(1-((a-b-d)/(e-f-g)))
}


# Average Euclidean Distance
AvEuclidean<-function(Comms){
  return(sqrt(sum(((Comms[1,]-Comms[2,])^2)/(sum(Comms[1,]+Comms[2,]>0)))))
}


# NESS 
# (normalised expected number of shared species) - when exponent m = 1 this is the same as the Morisita index.  Increasing m increases the weight of rare species

NESSm50<-function(Comms){
  p1<-Comms[1,]/(sum(Comms[1,]))
  p2<-Comms[2,]/(sum(Comms[2,]))
  mu1<-1-((1-p1)^50)
  mu2<-1-((1-p2)^50)
  NESS<-(2*sum(mu1*mu2))/(sum((mu1)^2)+sum((mu2)^2))
  return(1-NESS)
}


# true beta diversity, sensu Jost et al. 2007, Hill 1973

JostShannon<-function(Comms){ # 
  return(d(Comms, lev = "beta", q = 1, wts=c(sum(Comms[1,]), sum(Comms[2,])))-1) # q = 1 Shannon Entropy with weights to account for sample area or sampling effort 
} 
# the theoretical maximum is the number of communities, in this case 2.  The value will be 2 when they are completely distinct and  1 when they are identical  	
# the minus 1 puts comparisons between 2 communities on a scale of 0-1 rather than 1-2

JostSimpson<-function(Comms){
  return(d(Comms, lev = "beta", q = 2)-1) # q = 2 Simpson index 
}	

# additive partitioning sensu Lande 1996, using the formulae in Magurran and McGill (2011) pp 69-70


LandeShannon<-function(Comms){
  return(adipart(y=Comms, index="shannon"))
}


LandeSimpson<-function(Comms){
  return(adipart(y=Comms, index="simpson"))
}


# chao jaccard and chao Sorensen
# These metrics are designed to use abundance data to estimate the number of unseen shared species

chaojs<-function(Comms, method){
  shared<-as.matrix(Comms[,Comms[1,]>0 & Comms[2,]>0])
  if(ncol(as.matrix(shared))==0){
    return(1)
  } 
  if (ncol(as.matrix(shared))>0){ 
    In<-shared
    In[In!=1]<-0
    # U terms
    mterm<-(sum(Comms[2,])-1)/sum(Comms[2,])
    f_2<-ncol(as.matrix(shared[,shared[2,]==2 & shared[1,]>=1]))
    if(f_2==0){
      f_2<-1
    }
    ftermU<-ncol(as.matrix(shared[,shared[2,]==1 & shared[1,]>=1]))/(2*f_2)
    U<-sum(shared[1,]/sum(Comms[1,]))+(mterm*ftermU*sum(shared[1,]/sum(shared[1,])*In[1,]))
    if (U>1){
      U<-1
    }
    #V terms
    nterm<-(sum(Comms[1,])-1)/sum(Comms[1,])
    f2_<-ncol(as.matrix(shared[,shared[1,]==2 & shared[2,]>=1]))
    if(f2_==0){
      f2_<-1
    }
    ftermV<-ncol(as.matrix(shared[,shared[1,]==1 & shared[2,]>=1]))/(2*f2_)
    V<-sum(shared[2,]/sum(Comms[2,]))+(nterm*ftermV*sum((shared[2,]/sum(shared[2,]))*In[2,]))
    
    if(V>1){
      V<-1
    }
    
    if (method=="jaccard"){
      return(1-(U*V/(U+V-U*V))) # jaccard
    }
    if(method=="sorensen"){
      return(1-(2*U*V/(U+V)))
    }
  }
}

# function to return abundance-based metrics that are designed to partition the Bray-Curtis and Ruzicka metrics into 
# a) balanced variation in abundance (variation due to replacement of individuals by individuals of another species:the abundance equivalent of turnover)
# b) abundance gradients (variation due to the loss of individuals from one site: the abundance equivalent of nestedness)
# BC is an abundance-based extension of Classic Sorensen
# Ruzicka is an abundance-based extension of Classic Jaccard
part<-function (x) { # function to return metrics designed to partition the Ruzicka and Bray-Curtis metrics into turnover and nestedness components 
  A<-sum(pmin(x[1,],x[2,]))
  B<-sum(x[1,])-sum(pmin(x[1,],x[2,]))
  C<-sum(x[2,])-sum(pmin(x[1,],x[2,]))
  
  BaselgaBC_turn <- min(B, C)/(A+min(B,C)) #(Baselga 2013)
  BaselgaBC_nest <- (abs(B-C)/(2*A+(B+C)))*(A/(A+min(B,C))) #(Baselga 2013)
  BaselgaR_turn <- (2*min(B, C))/(A+2*min(B,C)) # (Legendre 2014)
  BaselgaR_nest <- (abs(B-C)/(A+B+C))*(A/(A+2*min(B,C))) # (Legendre 2014)
  
  PodaniBC_turn <- (2*min(B,C))/(2*A+(B+C)) # (Legendre 2014)
  PodaniBC_nest <- abs(B-C)/(2*A+B+C) # (Legendre 2014)
  PodaniR_turn <- 2*min(B,C)/(A+B+C) # (Podani et al. 2013)
  PodaniR_nest <- abs(B-C)/(A+B+C) # (Podani et al. 2013)
  
  
  results<-list(BaselgaBC_turn=BaselgaBC_turn, 
                BaselgaBC_nest=BaselgaBC_nest,
                BaselgaR_turn=BaselgaR_turn,
                BaselgaR_nest=BaselgaR_nest,
                PodaniBC_turn=PodaniBC_turn,
                PodaniBC_nest=PodaniBC_nest,
                PodaniR_turn=PodaniR_turn,
                PodaniR_nest=PodaniR_nest)
  return(results)
}
 


##########################################################################################################


