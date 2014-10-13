# Measuring taxonomic beta-diversity with abundance data

# Louise Barwell and Nick Isaac 
# 24/10/2011 to present
# Simulations to explore and score the properties and personalities of 25 beta-diversity metrics:
# Warning: on my computer, the entire script takes around 5 days to run 

######## Eleven desirable properties for all abundance-based metrics 

# D1)  unbiased by sample size
# D2)  unbiased by unequal sample sizes
# D3)  unbiased by alpha diversity
# D4)  additive under directional turnover
# D5)  multiplicative (probabilistic) under orthogonal sampling
# D6)  monotonic increase with species turnover
# D7)  monotonic increase with decoupling of species ranks 
# D8)  monotonic increase with differences in Shannon's Evenness
# D9)  value of beta under extreme decoupling of species ranks < beta when species turnover is complete 
# D10) value of beta when extreme differences in evenness < beta when species turnover is complete
# D11) Defined minima and maxima

####### eight 'personality' traits that capture some of the variation between metrics ########

# P1) Proportionality to species turnover
# P2) Sensitivity to differences in alpha diversity
# P3) Relative sensitivity to nestedness and turnover components of beta
# P4) Relative sensitivity to decoupling of species ranks and species turnover components of beta 
# P5) Relative sensitivity to evenness differences and species turnover components of beta
# P6) Decreased sensitivity to turnover in rare species under a positive ONR
# !P7)Sensitivity to species turnover
# !P8)Sensitivity to species loss (nestedness)
# !P9) How do metrics use abundance information?  How is the value of beta for turnover in a single species affected by the abundance of that species?

rm(list=ls())
library(ggplot2) 
library(reshape2) 
library(RCurl)
library(spatstat)

setwd("C:/Users/loubar/Dropbox/Manuscripts/Beta diversity/Beta_sims")

# create folders to organise the results by simulation
for (i in 1:11){
  dir.create(path=paste("D", i, sep=""))
}

for (i in c(1:6, 9, 10)){
  dir.create(path=paste("P", i, sep=""))
}



# source the functions in the file beta_diversity_funcs.r on github (https://github.com/loubar/Beta_diversity_sims/blob/master/Beta_diversity_funcs.r)
# these are functions to calculate the beta-diversity metrics and to simulate assemblage pairs under different scenarios
u <- "https://raw.githubusercontent.com/loubar/Beta_diversity_sims/master/Beta_diversity_funcs.r"
script <- getURL(u, ssl.verifypeer = FALSE)
eval(parse(text = script))


# give each metric an index
# 1:8 have no upper limit and will need to be plotted separately
metrics <- c(14,16,15,23, 20,17,18,1, 3, 4,6,5, 19,21,2,11,10, 9,12,13,22,24,25,7,8)
names(metrics)<- c('Jaccard abd.','Canberra','Bray-Curtis','Gower', 'Kulczynski','Morisita','Morisita-Horn','Euclidean', 'Manhattan', 'alt. Gower','Binomial', 'CYd','Horn','Renkonen','Av. Euclidean','Classic Jaccard', 'Classic Sorensen', 'sim','Chao Sorensen','Chao Jaccard', 'NESS','Jost Shannon','Jost Simpson','Lande Shannon', 'Lande Simpson')  



# generate 100 Fisher log series rank abundance distributions with ~10000 individuals and 100 species
# Use the average of 100 replicates.
# Comm1 is the starting assemblage in all of the following simulations
assemblages<-replicate(100, FisherRAD(N=10000, S=100))
Comm1<-round(rowMeans(assemblages), 0) #  Fisherian starting assemblage with ~ 10000 individuals and 100 species
#load('Comm1.rData')

######################################################################################################################################
# Simulations to test for each of the properties described above.
# Each metric is given a score to quantify how well it satisfies each property.

##################### D1) unbiased by sample size ######################################################################################

turnover<-(0:5*2)/10
N<-c(10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, sum(Comm1))
by.N<-list()
diffall<-list()
for (j in 1:length(turnover)){
  sim_refs<-replicate(10000, beta_metrics_all(sampling_effect(Comm1, N=sum(Comm1), p=turnover[j])))
  z<-list()
  diff<-list()
  for(k in 1:length(N)){
    sims<-replicate(10000, beta_metrics_all(sampling_effect(Comm1, N=N[k], p=turnover[j])))
    Metric<-rownames(sims)
    ID<-metrics
    sp_turnover<-turnover[j]
    N_size<-N[k]
    median<-sapply(1:nrow(sims), function(i) median(sims[i,], na.rm=T))
    uq<-as.numeric(sapply(1:nrow(sims), function(i) quantile(sims[i,], 0.75, na.rm=T)))
    lq<-as.numeric(sapply(1:nrow(sims), function(i) quantile(sims[i,], 0.25, na.rm=T)))
    min<-as.numeric(sapply(1:nrow(sims), function(i) min(cbind(sims[i,], sim_refs[i,]),na.rm=T)))
    max<-as.numeric(sapply(1:nrow(sims), function(i) max(cbind(sims[i,], sim_refs[i,]),na.rm=T)))
    z[[k]]<-data.frame(Metric, ID, sp_turnover, N_size, median, lq, uq, min, max, row.names=NULL)
    diff[[k]]<-t(sims-sim_refs)
  }
  a<-do.call(rbind,z)
  b<-do.call(rbind, diff)
  by.N[[j]]<-a
  diffall[[j]]<-b
}
D1_sim<-do.call(rbind, by.N)
D1_diffs<-do.call(rbind, diffall)

D1_diffs<-D1_diffs[,names(sort(metrics))]
ranges<-sapply(1:25, function(i) max(D1_sim[D1_sim$ID==i,"max"])-min(D1_sim[D1_sim$ID==i,"min"])) 

# score the performances using 10000reps*6turnover*7gradients*25 metrics = 10 500 000 simulations

# 1) the difference between median beta in fully censused assemblages and beta when assemblages are undersampled: score is the RMSE across all unique combinations of sample size and turnover
y<-acast(D1_sim, sp_turnover ~ N_size ~ ID, value.var='median')
D1_score<-sapply(1:25, function(i) sqrt(mean(((y[,16,i]-y[,,i])/(ranges[i]))^2)))
names(D1_score)<-names(sort(metrics))


D1_diffs<-sapply(1:ncol(D1_diffs), function(i) D1_diffs[,i]/ranges[i]) # first standardise the errors by the range of values for this simulation
colnames(D1_diffs)<-names(sort(metrics))



# add the scores into the data
for (i in 1:length(D1_score)){
  D1_sim[D1_sim$Metric==names(D1_score)[i],"D1_score"]<-paste(names(D1_score)[i], "\nD1 RMSE=", round(D1_score[i],4))
}
D1_sim$D1_score<-factor(D1_sim$D1_score, levels=sapply(1:length(sort(D1_score)), function(i) paste(names(sort(D1_score))[i], "\nD1 RMSE=", round(sort(D1_score)[i],4))))


# plot the response to undersampling both assemblages
a<-ggplot(D1_sim[D1_sim$ID>8,],aes(x=N_size,y=median, ymin=lq, ymax=uq, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D1_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.1, 1.1), breaks=((0:5)*2/10))+scale_x_continuous(limits=c(0, 10005))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))+theme(strip.text=element_text(size=8))
f<-d+ylab (expression(beta))+xlab(expression(italic(N)))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D1_sim$N), colour="black", linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#j
ggsave("D1/D1_st.png", height=8, width=6, dpi=500, pointsize=10)


a<-ggplot(D1_sim[D1_sim$ID<=8,],aes(x=N_size,y=median, ymin=lq, ymax=uq, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D1_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 10005))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))+theme(axis.title.y=element_text(angle=90))
f<-d+ylab (expression(beta))+xlab(expression(italic(N)))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D1_sim$N), colour="black", linetype=2)+theme(strip.text=element_text(size=8))
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#j
ggsave("D1/D1_unst.png", height=8/2, width=6, dpi=500, pointsize=10)


save(D1_diffs, file="D1/D1_diffs.rData") # error in all simulations
save(D1_sim, file="D1/D1_sim.rData") #  median, uq and lq of betas at each unqiue combination of sample size and turnover: this is for the plots
save(D1_score, file="D1/D1_score.rData") # RMSE

###########################################################################################################################################################################################################################

##################### D2) unbiased by unequal sample size ##################################################################################################################################################################

turnover<-(0:5*2)/10
N<-c(10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, sum(Comm1))
by.N<-list()
diffall<-list()
for (j in 1:length(turnover)){
  sim_refs<-replicate(10000, beta_metrics_all(unequal_sampling(Comm1, N=N[length(N)], p=turnover[j])))
  z<-list()
  diff<-list()
  for(k in 1:length(N)){
    sims<-replicate(10000, beta_metrics_all(unequal_sampling(Comm1, N=N[k], p=turnover[j])))
    Metric<-rownames(sims)
    ID<-metrics
    sp_turnover<-turnover[j]
    N_diff<-sum(Comm1)-N[k]
    median<-sapply(1:nrow(sims), function(i) median(sims[i,], na.rm=T))
    uq<-as.numeric(sapply(1:nrow(sims), function(i) quantile(sims[i,], 0.75, na.rm=T)))
    lq<-as.numeric(sapply(1:nrow(sims), function(i) quantile(sims[i,], 0.25, na.rm=T)))
    min<-as.numeric(sapply(1:nrow(sims), function(i) min(sims[i,],na.rm=T)))
    max<-as.numeric(sapply(1:nrow(sims), function(i) max(sims[i,],na.rm=T)))
    z[[k]]<-data.frame(Metric, ID, sp_turnover, N_diff, median, lq, uq, min, max, row.names=NULL)
    diff[[k]]<-t(sims-sim_refs)
  }
  a<-do.call(rbind,z)
  b<-do.call(rbind, diff)
  by.N[[j]]<-a
  diffall[[j]]<-b
}
D2_sim<-do.call(rbind, by.N)
D2_diffs<-do.call(rbind, diffall)

# D2_diffs are the differences between beta for fully censused assemblages and beta when one assemblage is undersampled 
# across all simulations 10000reps*16sample size differences*6 levels of turnover*25 metrics = 24 000 000 rows
D2_diffs<-D2_diffs[,names(sort(metrics))]
ranges<-sapply(1:25, function(i) max(D2_sim[D2_sim$ID==i,"max"])-min(D2_sim[D2_sim$ID==i,"min"])) 


# score the performance
# 1) differences in median beta: score as the RMSE
y<-acast(D2_sim, sp_turnover ~ N_diff ~ ID, value.var='median')
D2_score<-sapply(1:25, function(i) sqrt(mean(((y[,1,i]-y[,,i])/(ranges[i]))^2)))
names(D2_score)<-names(sort(metrics))


D2_diffs<-sapply(1:ncol(D2_diffs), function(i) D2_diffs[,i]/ranges[i]) # first standardise the errors by the range of values for this simulation
colnames(D2_diffs)<-names(sort(metrics))

# add the scores into the data
for (i in 1:length(D2_score)){
  D2_sim[D2_sim$Metric==names(D2_score)[i],"D2_score"]<-paste(names(D2_score)[i], "\nD2 RMSE =", round(D2_score[i],4))
}
D2_sim$D2_score<-factor(D2_sim$D2_score, levels=sapply(1:length(sort(D2_score)), function(i) paste(names(sort(D2_score))[i], "\nD2 RMSE =", round(sort(D2_score)[i],4))))


# plot the response to undersampling one assemblage

options(scipen=20)
a<-ggplot(D2_sim[D2_sim$ID>8,],aes(x=N_diff,y=median, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover), ymin=lq,ymax=uq)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D2_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.1), breaks=((0:5)*2/10))+scale_x_continuous(limits=c(0, 10000))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab(expression(beta))+xlab(expression(paste(Delta, italic(N), sep="")))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, colour="black", linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#j
ggsave("D2/D2_st.png", height=8, width=6, dpi=500, pointsize=10)

a<-ggplot(D2_sim[D2_sim$ID<=8,],aes(x=N_diff,y=median, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover), ymin=lq,ymax=uq)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D2_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 10000))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression(beta))+xlab(expression(paste(Delta, italic(N), sep="")))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, colour="black", linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#j
ggsave("D2/D2_unst.png", height=8/2, width=6, dpi=500, pointsize=10)

save(D2_diffs, file="D2/D2_diffs.rData") # differences in all simulations
save(D2_sim, file="D2/D2_sim.rData") # summary : median, uq and lq of betas at each unqiue combination of sample size and turnover
save(D2_score, file="D2/D2_score.rData") # difference in medians

#############################################################################################################################################################################################################

##################### D3) unbiased by alpha diversity ##################################################################################################################################################################

#  generate 10 assemblages with different alpha diversities 
#  N is fixed between 9990 and 10010 to avoid confounding effect of sample size and effect of alpha diversity
#  fixing N exactly takes much longer than setting tolerance limits of 10000+-10 
#  S is fixed as an integer divisible by 10 so that the number of species to turnover is always an integer


# The Fisher logseries describes a species abundance distribution.  Converting this to a rank abundance 
# distribution generates uncertainty as to the exact abundances of species in
# in each abundance class.
# Therefore, we use the average of 100 assemblages with S species and N individuals (using the function
# fisher.ecosystem) to
# obtain Fisherian assemblages with different alpha diversities.  

S<-c(300, 250, 200, 150, 100, 80, 60, 40, 20, 10)
turnover<-c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

# warning: the following loop takes in excess of 10 hours on my computer
alphas<-list() # a list of 10 aseemblages with different alpha-diversities are generated.
for (i in 1:length(S)){
  Spp<-as.data.frame(matrix(nrow=S[i], ncol=100))
  for (k in 1:100){
    repeat{
      z<-try(gen_alphas(N=10000, S=S[i])) # the try function catches the error in rep(1:sum(jj), rep(j, jj)): invalid "times" argument.
      # weirdly, this error only appears when I put the fisher.ecosystem function into a repeat loop (in function gen_alphas). 
      if(class(z)!="try-error") break
    }
    Spp[,k]<-z  
  }
  alphas[[i]]<-round(rowMeans(Spp), 0)
}  


by.alpha<-list() # each element in the list is for a different alpha.  Each is a data frame containing the names of the metrics in column 1, the value of alpha in col 2 and the value of beta for multiple levels of turnover
diffall<-list()
for (j in 1:length(turnover)){
  # first simulate under highest levels of alpha diversity
  sim_refs<-replicate(10000, beta_metrics_all(random_composition_change(Comm1=alphas[[1]], p=turnover[[j]])))
  z<-list()
  diff<-list()
  for(k in 1:length(alphas)){
    sims<-replicate(10000, beta_metrics_all(random_composition_change(Comm1=alphas[[k]], p=turnover[[j]])))
    Metric<-rownames(sims)
    ID<-metrics
    sp_turnover<-turnover[j]
    Alpha<-fishers.alpha(N=sum(alphas[[k]]), S=length(alphas[[k]]))
    median<-sapply(1:nrow(sims), function(s) median(sims[s,]))
    uq <- sapply(1:nrow(sims), function(s) quantile(sims[s,], probs=0.75)) # 25th percentile
    lq <- sapply(1:nrow(sims), function(s) quantile(sims[s,], probs=0.25)) # 75th percentile
    min <- sapply(1:nrow(sims), function(s) min(sims[s,]))
    max <- sapply(1:nrow(sims), function(s) max(sims[s,]))
    z[[k]]<-data.frame(Metric, ID, Alpha, sp_turnover, median, uq, lq, min, max, row.names=NULL)
    diff[[k]]<-t(sims-sim_refs)
  }
  a<-do.call(rbind,z)
  b<-do.call(rbind, diff)
  by.alpha[[j]]<-a
  diffall[[j]]<-b
}

D3_sim<-do.call(rbind, by.alpha)
D3_diffs<-do.call(rbind, diffall)

D3_diffs<-D3_diffs[,names(sort(metrics))]
ranges<-sapply(1:25, function(i) max(D3_sim[D3_sim$ID==i,"max"])-min(D3_sim[D3_sim$ID==i,"min"])) 

# score the performances
# 1) differences in median beta: take the mean across all unique levels of turnover and alpha diversity
x<-acast(D3_sim, sp_turnover ~ Alpha ~ ID, value.var='median')
D3_score<-sapply(1:25, function(i) sqrt(mean(((x[,10,i]-x[,,i])/ranges[i])^2)))
names(D3_score)<-names(sort(metrics))

# D3_diffs are the differences between beta when assemblages have highest alpha diversity and beta when alpha diversity is low 
# across all simulations 10000reps*10 levels of alpha-diversity*6 levels of turnover*25 metrics = 15 000 000 rows
D3_diffs<-sapply(1:ncol(D3_diffs), function(i) D3_diffs[,i]/ranges[i]) # first standardise the errors by the range of values for this simulation
colnames(D3_diffs)<-names(sort(metrics))


# add the scores into the data
for (i in 1:length(D3_score)){
  D3_sim[D3_sim$Metric==names(D3_score)[i],"D3_score"]<-paste(names(D3_score)[i], "\nD3 RMSE =", round(D3_score[i],4))
}
D3_sim$D3_score<-factor(D3_sim$D3_score, levels=sapply(1:length(sort(D3_score)), function(i) paste(names(sort(D3_score))[i], "\nD3 RMSE =", round(sort(D3_score)[i],4))))


##### Plot the response to equal alpha-diversity in assemblage pairs 


a<-ggplot(D3_sim[D3_sim$ID>8,],aes(x=Alpha,y=median, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover), ymin=lq,ymax=uq)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D3_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.01), breaks=((0:6)*2/10))+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression(beta))+xlab (expression(alpha[Fisher]))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D3_sim$Alpha), linetype=2, colour="black")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#j
ggsave("D3/D3_st.png", height=8, width=6, dpi=500, pointsize=10)

a<-ggplot(D3_sim[D3_sim$ID<=8,],aes(x=Alpha,y=median, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover), ymin=lq,ymax=uq)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D3_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression(beta))+xlab (expression(alpha[Fisher]))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")+theme(axis.title.y=element_text(angle=0))
h<-g+geom_vline(xintercept=max(D3_sim$Alpha), linetype=2, colour="black")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("D3/D3_unst.png", height=8/2, width=6, dpi=500, pointsize=10)

save(alphas, file='D3/alphas.rData')
save(D3_sim, file='D3/D3_sim.rData')
save(D3_diffs, file='D3/D3_diffs.rData')
save(D3_score, file='D3/D3_score.rData')

#################################################################################################################################################################################################################################################################


########################### D4) additive under directional turnover ########################################################

# generate a hypothetical environmental gradient across three assemblages
# species shared between assemblages 1 and 2 are more likely to be absent in assemblage three
# the probability of being turned over between assemblages 2 and 3 is much higher for species that were present in assemblage 1

grad<-c(1, 5, 10, 50, 100, 500, 1000) # how many times more likely are assemblage 1 species to be turned over between assemblages 2 and 3
turnover<-0:5/10
addtu<-list()
diffall<-list()
for (k in 1:length(grad)){ 
  add<-list()
  diff<-list()
  for (i in 1:length(turnover)){
    x1<-replicate(10000, add_comp_change(Comm1, p=turnover[i], gradient=grad[k]), simplify=F) 
    # gradient = Species shared between Comm1 and Comm2 are gradient times more likely to be turned over in Comm3
    # gradient specifies the steepness of the spatial gradient - the severity of spatial turnover
    diffadd<-t(sapply(1:length(x1), function(j) (beta_metrics_all(x1[[j]][1:2,])+beta_metrics_all(x1[[j]][2:3,]))-beta_metrics_all(x1[[j]][c(1,3),])))
    mediandiff<-sapply(1:ncol(diffadd), function(j) median(diffadd[,j], na.rm=T))
    uqdiff<-sapply(1:ncol(diffadd), function(j) quantile(diffadd[,j], 0.75, na.rm=T))
    lqdiff<-sapply(1:ncol(diffadd), function(j) quantile(diffadd[,j], 0.25, na.rm=T))
    y<-cbind(sapply(1:length(x1), function(j) beta_metrics_all(x1[[j]][c(1,3),])), sapply(1:length(x1), function(j) beta_metrics_all(x1[[j]][c(1,2),])+beta_metrics_all(x1[[j]][c(2,3),])))
    max<-sapply(1:nrow(y), function(j) max(y[j,]))
    min<-sapply(1:nrow(y), function(j) min(y[j,]))
    Metric<-names(metrics)
    ID<-metrics
    gr<-grad[k]
    sp_turnover<-turnover[i]
    add[[i]]<-data.frame(Metric, ID, gr, sp_turnover, mediandiff, uqdiff, lqdiff, min, max, row.names=NULL)
    diff[[i]]<-diffadd
  }
  addtu[[k]]<-do.call(rbind, add)
  diffall[[k]]<-do.call(rbind, diff)
}
D4_sim<-do.call(rbind, addtu)
D4_diff<-do.call(rbind, diffall)


D4_diff<-D4_diff[,names(sort(metrics))]
ranges<-sapply(1:25, function(i) max(D4_sim[D4_sim$ID==i,"max"])-min(D4_sim[D4_sim$ID==i,"min"])) 

D4_diff<-sapply(1:ncol(D4_diff), function(i) D4_diff[,i]/ranges[i]) # first standardise the errors by the range of values for this simulation
colnames(D4_diff)<-names(sort(metrics))

D4_RMSE<-sapply(1:ncol(D4_diff), function(i) sqrt(mean(D4_diff[,i]^2)))
names(D4_RMSE)<-names(sort(metrics))
# then take the median of the errors
D4_bias<-colMeans(D4_diff, na.rm=T)
names(D4_bias)<-names(sort(metrics))

# add the scores into the data, both RMSE and bias
for (i in 1:length(D4_bias)){
  D4_sim[D4_sim$Metric==names(D4_bias)[i],"D4_scores"]<-paste(names(D4_bias)[i], "\nD4 bias =", round(D4_bias[i],4), "\nD4 RMSE =", round(D4_RMSE[i],4))
}

#plot the difference between added beta (B_12+B23)  and observed beta (B_13)
a<-ggplot(D4_sim[D4_sim$ID>8,],aes(x=gr,y=mediandiff, ymin=lqdiff, ymax=uqdiff, group=factor(sp_turnover), fill=factor(sp_turnover), colour=factor(sp_turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D4_scores, nrow=6, ncol=3)+scale_x_log10(breaks=c(1,10,100,1000))+geom_hline(yintercept=0, linetype=2)+scale_y_continuous(limits=c(-0.4, 1), breaks=-2:5*2/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 8, colour = 'black'))+theme(axis.title.y =element_text(angle=0,size = 8, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression((beta["1,2"]+beta["2,3"])-beta["1,3"]))+xlab (expression(italic(g)))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
#h
ggsave("D4/D4_st.png", height=8, width=6, dpi=500, pointsize=10)

a<-ggplot(D4_sim[D4_sim$ID<=8,],aes(x=gr,y=mediandiff, ymin=lqdiff, ymax=uqdiff, group=factor(sp_turnover), fill=factor(sp_turnover), colour=factor(sp_turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D4_scores, scales="free_y", nrow=3, ncol=3)+scale_x_log10(breaks=c(1,10,100,1000))+geom_hline(yintercept=0, linetype=2)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 8, colour = 'black'))+theme(axis.title.y =element_text(angle=0,size = 8, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression((beta["1,2"]+beta["2,3"])-beta["1,3"]))+xlab (expression(italic(g)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")+labs(fill=expression(italic(t)), colour=expression(italic(t)))
h<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))+theme(strip.text=element_text(size=8))
#h
ggsave("D4/D4_unst.png", height=8/2, width=6, dpi=500, pointsize=10)

save(D4_sim, file='D4/D4_sim.rData')
save(D4_diff, file='D4/D4_diff.rData')
save(D4_RMSE, file='D4/D4_RMSE.rData')
save(D4_bias, file='D4/D4_bias.rData')

###########################################################################################################################################################################################


######################## D5) multiplicative (probabilistic) under orthogonal sampling ######################################################################################################

turnover<-(0:5)*2/10

mult<-list()
diffall<-list()
for (i in 1:length(turnover)){
  x<-replicate(10000, mult_comp_change(Comm1, p=turnover[i]), simplify=F)
  diffmult<-t(sapply(1:length(x), function(j) ((1-beta_metrics_all(x[[j]][1:2,]))*(1-beta_metrics_all(x[[j]][2:3,])))-(1-beta_metrics_all(x[[j]][c(1,3),]))))
  mediandiff<-sapply(1:ncol(diffmult), function(j) median(diffmult[,j], na.rm=T))
  uqdiff<-sapply(1:ncol(diffmult), function(j) quantile(diffmult[,j], 0.75, na.rm=T))
  lqdiff<-sapply(1:ncol(diffmult), function(j) quantile(diffmult[,j], 0.25, na.rm=T))
  y<-cbind(sapply(1:length(x), function(j) (1-beta_metrics_all(x[[j]][c(1,3),]))), sapply(1:length(x), function(j) (1-beta_metrics_all(x[[j]][c(1,2),]))*(1-beta_metrics_all(x[[j]][c(2,3),]))))
  max<-sapply(1:nrow(y), function(j) max(y[j,]))
  min<-sapply(1:nrow(y), function(j) min(y[j,]))
  Metric<-names(metrics)
  ID<-metrics
  sp_turnover<-turnover[i]
  mult[[i]]<-data.frame(Metric, ID, sp_turnover, mediandiff, uqdiff, lqdiff, min, max, row.names=NULL)
  diffall[[i]]<-diffmult
}
D5_sim<-do.call(rbind, mult)
D5_diff<-do.call(rbind, diffall)
# It is not possible to test whether similarity is multiplicative for unstandardised indices: they do not have a similraity complement
# so the unstandardised metrics are removed from the data set (this counts as a fail for this property)
D5_sim<-D5_sim[D5_sim$ID>8,]

D5_diff<-D5_diff[,names(sort(metrics))[9:25]]
ranges<-sapply(9:25, function(i) max(D5_sim[D5_sim$ID==i,"max"])-min(D5_sim[D5_sim$ID==i,"min"])) 


# score the performance  
D5_diff<-sapply(1:ncol(D5_diff), function(i) D5_diff[,i]/ranges[i]) # first standardise the errors by the range of values for this simulation
colnames(D5_diff)<-names(sort(metrics))[9:25]
# then take the median of the errors

D5_RMSE<-c(rep(NA, 8), sapply(1:ncol(D5_diff), function(i) sqrt(mean(D5_diff[,i]^2))))
names(D5_RMSE)<- names(sort(metrics))
D5_bias<-c(rep(NA, 8), colMeans(D5_diff, na.rm=T))
names(D5_bias)<-names(sort(metrics))

# add the scores into the data
for (i in 1:length(D5_bias)){
  D5_sim[D5_sim$Metric==names(D5_bias)[i],"D5_scores"]<-paste(names(D5_bias)[i], "\nD5 bias =", round(D5_bias[i],4), "\nD5 RMSE =", round(D5_RMSE[i],4))
}




# plot the difference between the value of beta expected under multiplicative behaviour and the observed beta
a<-ggplot(D5_sim,aes(x=sp_turnover,y=mediandiff, ymin=lqdiff, ymax=uqdiff))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D5_scores,scales="fixed", nrow=6, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))+geom_hline(linetype=2, yintercept=0)
d<-c+theme(axis.text.x=element_text(angle=0, hjust=1, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression((1-beta["1,2"])(1-beta["2,3"])-(1-beta["1,3"])))+xlab (expression(italic(t)))
#f
ggsave("D5/D5_st.png", height=8, width=6, dpi=500, pointsize=10)

save(D5_sim, file='D5/D5_sim.rData')
save(D5_diff, file='D5/D5_diff.rData')
save(D5_bias, file='D5/D5_bias.rData')
save(D5_RMSE, file='D5/D5_RMSE.rData')
##########################################################################################################################################################################################


################# # D6) monotonic increase with species turnover ###################################################################################################################

S<-(10:1)*10
p<-0:10/10

n.and.t<-list() 
for(j in 1:length(p)){
  nestedness<-list()
  for (i in 1:length(S)){
    sims<-replicate(10000, nested(S=S[i], p=p[j], Comm1), simplify=F)
    all_s<-sapply(1:length(sims), function(s) beta_metrics_all(sims[[s]]))
    Beta<-apply(all_s, 1, FUN=median)
    upper<-sapply(1:nrow(all_s), function(k) quantile(all_s[k,], probs=0.75))
    lower<-sapply(1:nrow(all_s), function(k) quantile(all_s[k,], probs=0.25))
    min<-sapply(1:nrow(all_s), function(k) min(all_s[k,]))
    max<-sapply(1:nrow(all_s), function(k) max(all_s[k,]))
    Metric<-names(Beta)
    ID<-metrics
    sp_turnover=p[j]
    sp_loss=S[1]-S[i]
    nestedness[[i]]<-data.frame(Metric, ID, sp_turnover, sp_loss, Beta, upper, lower, min, max, row.names=NULL)   
  }
  nestedness<-do.call(rbind, nestedness)
  n.and.t[[j]]<-nestedness
}
D6_sim<-do.call(rbind, n.and.t)

# Score the performance
x<-acast(D6_sim, sp_turnover ~ sp_loss ~ ID, value.var="Beta")
# desirable property D6 monotonic increase with species turnover (strictly increasing: if turnover1,2 > turnover3,4, then beta1,2 > beta3,4) 
# score: TRUE / FALSE
D6_score<-sapply(1:25, function(i) all(sapply(1:dim(x)[2], function(j) rank(x[,j,i])==1:length(p))))
D6_score<-as.character(D6_score)
names(D6_score)<-names(sort(metrics))

# add the scores in to the data.frame
for (i in 1:length(D6_score)){
  D6_sim[D6_sim$Metric==names(D6_score)[i],"D6_score"]<-paste("D6 =", D6_score[i]) 
}

save(D6_sim, file="D6/D6_sim.rData")
save(D6_score, file="D6/D6_score.rData")
# plot this later including the the scores for P2 relative sensitivity to nestedness and species turnover 

#############################################################################################################################################################################################################################################################


########### D7)  monotonic increase with decoupling of species ranks ###################################################################################################################################################### 


turnover<-(0:5*2)/10
r<-c(-10:10/10)

partial.r<-list()


for (k in 1:length(turnover)){
  all.r<-list()
  for(i in 1:length(r)){
    x<-replicate(10000, beta_metrics_all(par_cor(Comm1, r=r[i], p=turnover[k]))) # 1000 random samples from a population with rank correlation coefficient r (with Comm1)
    Metric<-rownames(x)
    ID<-metrics
    p.turnover<-rep(turnover[k], length(Metric))
    partial.cor<-rep(r[i], length(Metric))
    Beta<-sapply(1:nrow(x), function(j) median(x[j,]))
    upper<-sapply(1:nrow(x), function(m) quantile(x[m,], probs=0.75))
    lower<-sapply(1:nrow(x), function(n) quantile(x[n,], probs=0.25))
    all.r[[i]]<-data.frame(Metric, ID, p.turnover, partial.cor, Beta, lower, upper, row.names=NULL)
  }
  allr2<-do.call(rbind, all.r)
  partial.r[[k]]<-allr2
}
D7_sim<-do.call(rbind, partial.r)


# Score the performances (TRUE/FALSE)
x<-acast(D7_sim, p.turnover ~ partial.cor ~ ID, value.var="Beta")
D7_score<-as.character(sapply(1:25, function(i) all(sapply(1:4, function(j) order(x[j,,i])==length(r):1))))
names(D7_score)<-names(sort(metrics))

for (i in 1:length(D7_score)){
  D7_sim[D7_sim$Metric==names(D7_score)[i],"D7_score"]<-paste("D7 =", D7_score[i]) 
}

save(D7_score, file="D7/D7_score.rData")
# plot later with D9 and P3

##########################################################################################################################################

########################## D8) monotonic increase with differences in Shannon's Evenness ###############################################

Comm1[1]<-Comm1[1]+1 # make Comm1 a round 10000


turnover<-(0:5*2)/10

even<-list() #  create a list of starting assemblages with different levels of evenness
pr<-c(1+-4:5*2/10, 4, 6, 8) # the power to raise the abundances in the starting assemblage to

even[[1]]<-rep(100, 100) # perfectly even
even[[15]]<-c(sum(Comm1)-sum(rep(1, length(Comm1)-1)), rep(rep(1, length(Comm1)-1)))# extremely uneven (all except the dominant species have just one individual)
for(i in 1:length(pr)){
  even[[i+1]]<-Comm1^pr[i] # and some evennesses in between.
} 

ed<-rep(NA, length(even))
for (i in 1:length(even)){ # obtain the mean evenness difference across 10000 simulations at each level of evenness
  comms<-replicate(10000, evenness(Comm1, turnover=0, even=even[[i]]), simplify=F)
  ed[i]<-mean(sapply(1:length(comms), function(j) abs(diversity(comms[[j]][1,], index="shannon")/log(length(Comm1))-diversity(comms[[j]][2,], index="shannon")/log(length(Comm1)))))
}


Evenness<-list()

for (j in 1:length(turnover)){
  z<-list()
  for(k in 1:length(even)){
    sim<-replicate(10000, evenness(Comm1, even=even[[k]], turnover=turnover[j]))
    bma<-sapply(1:dim(sim)[3], function(i) beta_metrics_all(sim[,,i]))
    Beta<-apply(bma, 1, median)
    upper<-sapply(1:nrow(bma), function(i) quantile(bma[i,], probs=0.75, na.rm=T))
    lower<-sapply(1:nrow(bma), function(i) quantile(bma[i,], probs=0.25, na.rm=T))
    sp_turnover<-rep(turnover[j], length(metrics))
    Even_diff<-ed[k]
    Metric<-rownames(bma)
    ID<-metrics
    z[[k]]<-data.frame(Metric, ID, Even_diff, sp_turnover, Beta, upper, lower, row.names=NULL)
  }
  a<-do.call(rbind,z)
  Evenness[[j]]<-a
}

D8_sim<-do.call(rbind, Evenness) 
# score the performances


# score: TRUE / FALSE
x<-acast(D8_sim, sp_turnover ~ Even_diff ~ ID, value.var="Beta")
D8_score<-as.character(sapply(1:25, function(i) all(sapply(1:5, function(j) (rank(x[j,,i])==1:length(even))))))
names(D8_score)<-names(sort(metrics))

for (i in 1:length(D8_score)){
  D8_sim[D8_sim$Metric==names(D8_score)[i],"D8_score"]<-paste("D8 =", D8_score[i]) 
}

save(D8_score, file="D8/D8_score.rData")
#plot later with D10 and P4

########################################################################################################################################################################################


##################################### D9) value of beta under extreme decoupling of species ranks < beta when species turnover is complete #########################################

# use the simulations in D7


x<-acast(D7_sim, p.turnover ~ partial.cor ~ ID, value.var="Beta")
D9_score<-as.character(sapply(1:25, function(i) x[1,1,i]<x[6,length(r),i]))
names(D9_score)<-names(sort(metrics))

# add scores into D7_sim data frame
for (i in 1:length(D9_score)){
  D7_sim[D7_sim$Metric==names(D9_score)[i],"D9_score"]<-paste("D9 =", D9_score[i]) 
}

D7_sim$title<-paste(D7_sim$Metric, "\n", D7_sim$D7_score, "\n", D7_sim$D9_score, sep="")

a<-ggplot(D7_sim[D7_sim$ID>8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)+scale_y_continuous(limits=c(-0.04, 1.03), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(r)))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=1, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("D7/D7_st.png", height=8, width=6, dpi=500, pointsize=10)


a<-ggplot(D7_sim[D7_sim$ID<=8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(r)))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=1, linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("D7/D7_unst.png", height=8/2, width=6, dpi=500, pointsize=10)

save(D9_score, file="D9/D9_score.rData")

######################################################################################################################################################################################


########################### D10) value of beta under extreme differences in evenness < beta when species turnover is complete ######################################################


# use the evenness difference simulations in D8

# score: TRUE / FALSE
x<-acast(D8_sim, sp_turnover ~ Even_diff ~ ID, value.var="Beta")
D10_score<-as.character(sapply(1:25, function(i) x[1,length(even),i]<x[6,1,i]))                                                                         
names(D10_score)<-names(sort(metrics))


# add the scores into the D8_sim data frame

for (i in 1:length(D10_score)){
  D8_sim[D8_sim$Metric==names(D10_score)[i],"D10_score"]<-paste("D10 =", D10_score[i]) 
}

save(D10_score, file="D10/D10_score.rData")

#######################################################################################################################################################################################


############################# D11) Defined minima and maxima ############################################################################################################ 


# see Table 1
D11_score<-as.character(as.logical(c(rep(0, 8), rep(1, 17))))
names(D11_score)<-names(sort(metrics))
save(D11_score, file="D11/D11_score.rdata")
##############################################################################################################################################################



##############################################################################################################################################################

############################################## P1)   Sensitivity to differences in alpha diversity ################################


turnover<-(0:5*2)/10
diff_alpha<-list() 
for (j in 1:length(turnover)){
  z<-list()
  
  for(i in 1:length(alphas)){
    p.turnover<-rep(turnover[j], length(metrics)) # build a data frame for each unique combination of alpha difference and turnover  
    Alpha_D<-rep(fishers.alpha(N=sum(alphas[[1]]), S=length(alphas[[1]]))-fishers.alpha(N=sum(alphas[[i]]), S=length(alphas[[i]])), length(metrics))
    sim.1000<-replicate(10000, beta_metrics_all(Diff_alpha(alphas[[1]], alphas[[i]], turnover=turnover[j])))
    Beta<-sapply(1:nrow(sim.1000), function(s) median(sim.1000[s,]))
    upper<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.75)) # 25th percentile
    lower<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.25)) # 75th percentile
    max<-sapply(1:nrow(sim.1000), function(s) max(sim.1000[s,]))
    min<-sapply(1:nrow(sim.1000), function(s) min(sim.1000[s,]))
    Metric<-rownames(sim.1000)
    ID<-metrics
    z[[i]]<-data.frame(Metric, ID, Alpha_D, p.turnover, Beta, upper, lower, max, min, row.names=NULL)
  } # stick these data frames  for each level of alpha diversity difference together
  a<-do.call(rbind,z)
  diff_alpha[[j]]<-a
}

P1_sim<-do.call(rbind, diff_alpha)

# score the perfomance
ranges<-sapply(1:25, function(i) max(P1_sim[P1_sim$ID==i,"max"])-min(P1_sim[P1_sim$ID==i,"min"])) 


# score: mean difference between median beta and a reference level of beta diversity under equal alpha diversity (see Table 3)
y<-acast(P1_sim, p.turnover ~ Alpha_D ~ ID, value.var='Beta')
P1_score<-sapply(1:25, function(i) sqrt(mean(((y[,1,i]-y[,,i])/(ranges[i]))^2)))
names(P1_score)<-names(sort(metrics))


# add the scores into the data.frame P1_sim
for (i in 1:length(P1_score)){
  P1_sim[P1_sim$Metric==names(P1_score)[i],"P1_score"]<-paste(names(P1_score)[i], "\nP1 RMSE =", round(P1_score[i],4))
}
P1_sim$P1_score<-factor(P1_sim$P1_score, levels=sapply(1:length(sort(P1_score)), function(i) paste(names(sort(P1_score))[i], "\nP1 RMSE =", round(sort(P1_score)[i],4))))

# plot the responses to differences in alpha-diversity

a<-ggplot(P1_sim[P1_sim$ID>8,],aes(x=Alpha_D,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~P1_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.01), breaks=((0:6)*2/10))+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste(Delta, alpha[Fisher], sep="")))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=min(P1_sim$Alpha_D), linetype=2, colour="black")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("P1/P1_st.png", height=8, width=6,dpi=500, pointsize=10)

a<-ggplot(P1_sim[P1_sim$ID<=8,],aes(x=Alpha_D,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~P1_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste(Delta, alpha[Fisher], sep="")))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=min(P1_sim$Alpha_D), linetype=2, colour="black")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
ggsave("P1/P1_unst.png", height=8/2, width=6,dpi=500, pointsize=10)

save(P1_sim, file="P1/P1_sim.rData")
save(P1_score, file="P1/P1_score.rData")

############################################################################################################################################################################################################################################################



################################ P2) Relative sensitivity to nestedness and turnover components of beta #############################################################


# use the simulations in D6_sim

# score  the performmaces
# personality iv) Relative sensitivity to species loss and species turnover
# score: median beta for extreme species loss (s=90) / median beta for complete turnover (t=1)  
x<-acast(D6_sim, sp_turnover ~ sp_loss ~ ID, value.var="Beta")
P2_score<-sapply(1:25, function(i) x[1,10,i]/x[6,1,i])
names(P2_score)<-names(sort(metrics))


# add the scores to the D6_sim data frame
for (i in 1:length(P2_score)){
  D6_sim[D6_sim$Metric==names(P2_score)[i],"P2_score"]<-paste("P2 =", round(P2_score[i],4))
}


D6_sim$title<-paste(D6_sim$Metric, "\n", D6_sim$D6_score, "\n", D6_sim$P2_score, sep="")


a<-ggplot(D6_sim[D6_sim$ID>8 & D6_sim$sp_loss %in% c(0, 10, 30, 50, 70, 90),],aes(x=sp_turnover,y=Beta, ymin=lower, ymax=upper, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=(0:5*2/10))+scale_y_continuous(limits=c(-0.04, 1.101), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(t)))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=0, linetype=2)
j<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
ggsave("P2/P2_st.png", height=8, width=6,dpi=500, pointsize=10)

#unstandardised metrics
#x11()
a<-ggplot(D6_sim[D6_sim$ID<=8 & D6_sim$sp_loss %in% c(0, 10, 30, 50, 70, 90),],aes(x=sp_turnover,y=Beta, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss), ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title,scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=(0:5*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(t)))+labs(colour=expression(italic(s)), fill=expression(italic(s)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=0, linetype=2)
j<-g+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("P2/P2_unst.png", height=8/2, width=6,dpi=500, pointsize=10)

save(P2_score, file="P2/P2_score.rData")
save(D6_sim, file="P2/D6_sim.rData")
###############################################################################################################################################################################


################################ P3) Relative sensitivity to decoupling of species ranks and species turnover components of beta ######################################



# use the simulations in D7_sim

# score the performances
# score: value of beta for extreme decoupling of ranks / value of beta for complete species turnover.  
x<-acast(D7_sim, p.turnover ~ partial.cor ~ ID, value.var="Beta")
P3_score<-sapply(1:25, function(i) x[1,1,i]/x[6,length(r),i])
names(P3_score)<-names(sort(metrics))


# add the scores into the data.frame D7_sim
for (i in 1:length(P3_score)){
  D7_sim[D7_sim$Metric==names(P3_score)[i],"P3_score"]<-paste("P3 =", round(P3_score[i],4))
}
D7_sim$title<-paste(D7_sim$Metric, "\n", D7_sim$D7_score, "\n", D7_sim$D9_score, "\n", D7_sim$P3_score, sep="")
save(D7_sim, file="P3/D7_sim.rData")

a<-ggplot(D7_sim[D7_sim$ID>8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)+scale_y_continuous(limits=c(-0.04, 1.03), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(r)))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=1, linetype=2)
j<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("P3/P3_st.png", height=8, width=6,dpi=500, pointsize=10)


a<-ggplot(D7_sim[D7_sim$ID<=8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(r)))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=1, linetype=2)
j<-g+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("P3/P3_unst.png", height=8/2, width=6,dpi=500, pointsize=10)

save(P3_score, file="P3/P3_score.rData")
save(D7_sim, file="P3/D7_sim.rData")

################################ P4) Relative sensitivity to evenness differences and species turnover components of beta ######################################

# use the simulations in D8_sim

# score the performances
# score: value of beta for extreme decoupling of ranks / value of beta for complete species turnover
x<-acast(D8_sim, sp_turnover ~ Even_diff ~ ID, value.var="Beta")
P4_score<-sapply(1:25, function(i) x[1,length(even),i]/x[6,1,i])
names(P4_score)<-names(sort(metrics))


for (i in 1:length(P4_score)){
  D8_sim[D8_sim$Metric==names(P4_score)[i],"P4_score"]<-paste("P4 =", round(P4_score[i],4))
}
D8_sim$title<-paste(D8_sim$Metric, "\n", D8_sim$D8_score, "\n", D8_sim$D10_score, "\n", D8_sim$P4_score, sep="")



a<-ggplot(D8_sim[D8_sim$ID>8,],aes(x=Even_diff,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0, 1), breaks=0:5*2/10)+scale_y_continuous(limits=c(0, 1), breaks=(0:5*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste(Delta, E[Shannon], sep="")))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("P4/P4_st.png", height=8, width=6,dpi=500, pointsize=10)

a<-ggplot(D8_sim[D8_sim$ID<=8,],aes(x=Even_diff,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 1), breaks=(0:5*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste(Delta, E[Shannon], sep="")))+labs(fill=expression(italic(t)), colour=expression(italic(t)))
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("P4/P4_unst.png", height=8/2, width=6,dpi=500, pointsize=10)

save(P4_score, file="P4/P4_score.rData")
save(D8_sim, file="P4/D8_sim.rData")

####################################################################################################################################################################################################

################################ P5: relative sensitivity to turnover in rare versus common species ##########################
# a simulation to help visualise how  metrics use abundance information.

# turn over each species in turn.  How does the abundance of a single species turned over influence the value of beta?
Z<-sapply(1:length(Comm1), function(i) beta_metrics_all(abd_sens(Comm1, sp.rank=i)))
colnames(Z)<-Comm1
Z<-melt(Z)
colnames(Z)<-c("Metric", "Abundance", "Beta")
Z$Beta_abundance<-Z$Beta/Z$Abundance
Z$ID<-rep(metrics, 100)
Z$Abundance<-as.numeric(Z$Abundance)
Z$Metric<-reorder(Z$Metric,Z$ID)
Z$Rel_abd<-Z$Abundance/sum(Comm1)
Z$Beta_rel_abd<-round(Z$Beta/Z$Rel_abd, 3)

Z$Metric<-reorder(Z$Metric,Z$ID)

P5_sim<-Z
#P5_sim$Beta[P5_sim$Beta==0]<-0.000000001

P5_score<-rep(NA, length(metrics))
for (i in 1:25){
  P5_score[i]<-mean(P5_sim[P5_sim$ID==i & P5_sim$Abundance==min(P5_sim$Abundance),"Beta"])/mean(P5_sim[P5_sim$ID==i & P5_sim$Abundance==max(P5_sim$Abundance),"Beta"])
}
names(P5_score)<-names(sort(metrics))
save(P5_score, file="P5/P5_score.rData")

# scoring as the max/ min would be preferable (e.g. the score represents how many more times graeter is the value of beta for a common versus a rare species.  However, Morisita and some other care so little about rare species that the value is zero and we can't divide by it)

P5_sim$score<-NA
for (i in 1:length(metrics)){
  P5_sim[P5_sim$Metric==names(metrics)[i],"score"]<-round(P5_score[names(metrics)[i]], 4)
}

P5_sim$title<-paste(P5_sim$Metric, "\n", "P5 =", P5_sim$score)


# plot the response to a single species turned over 
a<-ggplot(P5_sim[P5_sim$ID>8,],aes(x=Rel_abd,y=Beta)) 
b<-a+geom_line(size=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))+xlim(0, 0.155)+ylim(0, 0.4)
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab(expression(italic(n)))
h<-f+theme(legend.position="none")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("P5/P5_st.png", height=8, width=6,dpi=500, pointsize=10)


a<-ggplot(P5_sim[P5_sim$ID<=8,],aes(x=Rel_abd,y=Beta)) 
b<-a+geom_line(size=0.5)+facet_wrap(~Metric, scales="free", nrow=3, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab(expression(italic(n)))
h<-f+theme(legend.position="none")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#j
ggsave("P5/P5_unst.png", height=8/2, width=6,dpi=500, pointsize=10)

save(P5_sim, file="P5/P5_sim.rData")

###############################################################################################################################################################################################################

############################# !P9) Decreased sensitivity to turnover in rare species under a positive ONR #######################################################################################################################

turnover<-0:10/10

sep.p<-list()

for (i in 1:length(turnover)){
  x<-replicate(10000, beta_metrics_all(rare_composition_change(Comm1, p=turnover[i])))
  Metric<-rownames(x)
  Beta<-sapply(1:nrow(x), function(j) median(x[j,]))
  ID<-metrics
  p.turnover<-rep(turnover[i], nrow(x))
  upper<-sapply(1:nrow(x), function(j) quantile(x[j,], probs=0.75))
  lower<-sapply(1:nrow(x), function(j) quantile(x[j,], probs=0.25))
  max<-sapply(1:nrow(x), function(j) max(x[j,]))
  min<-sapply(1:nrow(x), function(j) min(x[j,]))
  sep.p[[i]]<-data.frame(Metric, ID, p.turnover, Beta, upper, lower, max, min)
}

P9_sim<-do.call(rbind, sep.p)
rownames(P9_sim)<-NULL


# add the reference value of beta in the absence of a positive ONR (plus the upper and lower quartiles)
for (i in 1:25){
  P9_sim[P9_sim$ID==i,"lower_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"lower"]
  P9_sim[P9_sim$ID==i,"upper_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"upper"]
  P9_sim[P9_sim$ID==i,"Beta_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"Beta"]
  P9_sim[P9_sim$ID==i,"max_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"max"]
  P9_sim[P9_sim$ID==i,"min_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"min"]
}

P9_score<-sapply(1:25, function(i) with(P9_sim[P9_sim$ID==i,], sqrt(mean(((Beta_Test1-Beta)/(max(c(max_Test1, max))-min(c(min_Test1, min))))^2))))
names(P9_score)<-names(sort(metrics))


# add the scores into the data.frame P9_sim
for (i in 1:length(P9_score)){
  P9_sim[P9_sim$Metric==names(P9_score)[i],"P9_score"]<-paste(names(P9_score)[i], "\nP9 =", round(P9_score[i],4))
}
P9_sim$P9_score<-factor(P9_sim$P9_score, levels=sapply(1:length(sort(P9_score)), function(i) paste(names(sort(P9_score))[i], "\nP9 =", round(sort(P9_score)[i],4))))

# plot the decreased sensitivty to turnover in rare species under a positive ONR

# standardised
a<-ggplot(data=P9_sim[P9_sim$ID>8,],aes(x=p.turnover,y=Beta, ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.3)+facet_wrap(~P9_score,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5)*2)/10)+scale_y_continuous(limits=c(-0.04,1.101), breaks=((0:5)*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(italic(t)))
g<-f+geom_line(aes(p.turnover, Beta_Test1), linetype=2)+geom_ribbon(aes(ymin=lower_Test1, ymax=upper_Test1), colour=NA, alpha=0.3)
h<-g+theme(legend.position="none")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
#j
ggsave("P9/posONR_st.png", height=8, width=6, units="in") 

# unstandardised
a<-ggplot(data=P9_sim[P9_sim$ID<=8,],aes(x=p.turnover,y=Beta, ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.3)+facet_wrap(~P9_score,scales="free_y", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5)*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab (expression(italic(t)))
g<-f+geom_line(aes(p.turnover, Beta_Test1), linetype=2)+geom_ribbon(aes(ymin=lower_Test1, ymax=upper_Test1), colour=NA, alpha=0.3)
h<-g+theme(legend.position="none")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
ggsave("P9/posONR_unst.png", height=8/2, width=6, units="in") 

save (P9_sim, file="P9/P9_sim.rData")
save(P9_score, file="P9/P9_score.rData")

#############################################################################################################

############### !P10 Scale-dependence of beta-diversity metrics ####################################

# To test for scale-dependence, we simulate a patchy distribution using the Thomas point process, an inhomogeneous 
# Poisson point process within a square representing our study region.  This generates turnover implicitly.  
# We then sample from this study region using pairs of equal sized quadrats of width 'scales' with 100 repliactes for each.  The position of these quadrats is chosen at random from within the study region. 
# Calculate beta-diversity for each pair of quadrats and investigate how beta changes with scale.

# This test is closely related to the issue of relative sensitivity to rare versus common species.  We predict that
# metrics that are largely insensitive to turnover in rare species are likely to be highly independent of spatial scale. 
# because common species are widespread and are likely to be present in almost all quadrats, even those at the fine grain that, by definition have small sample sizes.
# By contrast rare species are very much more likely to be turned over within a pair of fine grain samples, while they are more likely to be present at a coarse grains.
# Therefore a metric that is sensitive to turnover on rare species is predicted to be scale-dependent, while a metric that cares only about common species will find very little difference between samples of diffrent spatial grain
# However, our results show that while Morisita (only cares about common sp turnover) is, on average, independent of scale (as predicted), it also shows
# a large amont of variation at fine grains.
# I think this is probably because the common species will be in almost all fine grain samples.  However, over multiple simulations, we eventually sample a quadrat pair where the common species is absent in one quadrat and present in the other
# So beta-diversity is occasionally very high, but on average very low, as only the common specioes matter to Morisita.





scales <- c(1,2,3,4,5,6,7,8,9,10, 20, 50)
# scale is the length in units of a square quadrat in a single simulation

# the intensity of the Thomas process is kappa*mu.  This is an inhomogeneous Poisson point process
# kappa is the intensity of a Poisson distribution describing the distribution of cluster centres 
# sigma is the diameter of each cluster

# To fix(ish) the number of individuals mu must be dependent on kappa and vice versa
# Here, kappa can be any where between 1 and the number of individuals of a given species
# and mu is the number of individuals/ kappa
agg<-function(Comm){
  ka <- runif(1, min=1, max=Comm)
  a<-rThomas(kappa=ka, sigma=runif(1, min=0.0001, max=0.9), mu=Comm/ka)
}

scale_agg<-function(Comm1, scale){
  y<-lapply(1:length(Comm1), function(i) agg(Comm1[i]))
  qcounts<-lapply(1:length(y), function(i) as.matrix(quadratcount(y[[i]], nx=scale, ny=scale)))
  q1<-sapply(1:length(qcounts), function(i) qcounts[[i]][sample(1:nrow(qcounts[[i]]), 1), sample(1:ncol(qcounts[[i]]), 1)])
  q2<-sapply(1:length(qcounts), function(i) qcounts[[i]][sample(1:nrow(qcounts[[i]]), 1), sample(1:ncol(qcounts[[i]]), 1)])
  z<-rbind(q1, q2)
  colnames(z)<-NULL
  return(z)
}

# now simulate this a few times and record how much each metric changes with spatial scale
scale_effect<-list()
for (i in 1:length(scales)){
  scaledata<-replicate(100, beta_metrics_all(scale_agg(Comm1, scale=scales[i])))
  Beta<-sapply(1:nrow(scaledata), function(j) median(scaledata[j,], na.rm=T))
  uq<-sapply(1:nrow(scaledata), function(j) quantile(scaledata[j,], prob=0.75,na.rm=T))
  lq<-sapply(1:nrow(scaledata), function(j) quantile(scaledata[j,], prob=0.25, na.rm=T))
  scale_effect[[i]]<-data.frame(Metric=rownames(scaledata), ID=metrics, Scale=scales[i], Beta=Beta, uq=uq, lq=lq)
}
scaling<-do.call(rbind, scale_effect)

# plot the effect of spatial scale on each of the metrics.  All should tend towards zero at the coarsest scales as this constitutes th entire 
a<-ggplot(scaling[scaling$ID>8,],aes(x=1/Scale,y=Beta, ymin=lq, ymax=uq))# group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~Metric, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.1, 1.1), breaks=((0:5)*2/10))+scale_x_continuous(limits=c(0, 1))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))+theme(strip.text=element_text(size=8))
f<-d+ylab (expression(beta))+xlab("Sample grain (proportion of region in each quadrat)")+theme(legend.position="none")
g<-f+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#f
ggsave("P10/scale_st.png", height=8, width=6, units="in")

a<-ggplot(scaling[scaling$ID<=8,],aes(x=1/Scale,y=Beta, ymin=lq, ymax=uq))# group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~Metric, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 1))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))+theme(strip.text=element_text(size=8))
f<-d+ylab (expression(beta))+xlab("Sample grain (proportion of region in each quadrat)")+theme(legend.position="none")
g<-f+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
#f
ggsave("P10/scale_unst.png", height=8/2, width=6, units="in")

save(scaling, file="P10/scaling.rData")

##################################################################################################################################################################################################################################################
load("D1/D1_score.rData")
load("D2/D2_score.rData")
load("D3/D3_score.rData")
load("D4/D4_RMSE.rData")
load("D5/D5_RMSE.rData")
load("D4/D4_bias.rData")
load("D5/D5_bias.rData")
load("D6/D6_score.rData")
load("D7/D7_score.rData")
load("D8/D8_score.rData")
load("D9/D9_score.rData")
load("D10/D10_score.rData")
load("D11/D11_score.rData")
load("P1/P1_score.rData")
load("P2/P2_score.rData")
load("P3/P3_score.rData")
load("P4/P4_score.rData")
load("P5/P5_score.rData")
##### Results figures and tables

# Table 3
D<-data.frame(D1=D1_score, D2=D2_score, D3=D3_score, D4=D4_RMSE, D5=D5_RMSE, D6=D6_score, D7=D7_score, D8=D8_score, D9=D9_score, D10=D10_score, D11=D11_score, check.rows=TRUE, stringsAsFactors=F)

# Table 4
P<-data.frame(P1=P1_score, P2=P2_score, P3=P3_score, P4=P4_score, P5=P5_score)
P<-round(P,4)

# Fig. 1
DP<-cbind(D[,1:4], P)
pca1<-prcomp(DP, scale=T)
summary(pca1)
plot(pca1)

# pca axes 1 and 2

tiff("pca.tif", height=6, width=6, res=500, units='in', pointsize=10)
  par(xpd=NA)
  biplot(pca1, las=1, tck=0.02, cex=0.8,  col=c("black", "dark grey"), asp=1.8, cex.axis=0.9)
  #mtext(outer=T, "a)", line=-6, side=3, adj=0.1)
dev.off()

# pca axis 2 and 3
par(xpd=NA)
biplot(pca1, choices=2:3, las=1, tck=0.02, cex=0.8,  col=c("black", "dark grey"), asp=1.8, cex.axis=0.9)
mtext(outer=T, "b)", line=-6, side=3, adj=0.1)

# pca axes 1 and 3
par(xpd=NA)
biplot(pca1, choices=c(1,3), las=1, tck=0.02, cex=0.8,  col=c("black", "dark grey"), asp=1.8, cex.axis=0.9)
mtext(outer=T, "b)", line=-6, side=3, adj=0.1)



############################################################################################################################################################################################################################################################
# Are any of the metrics outperformed by at least one other metric on all desirable properties?
# These are Pareto-dominated

# convert all to numeric variables so they can be ranked
# TRUE gets 1, FALSE gets 2
# NAs on test 5 get 2 (FAIL)
pareto<-D
#pareto$D5bias[is.na(pareto$D5bias)]<-'FALSE'
pareto$D5[is.na(pareto$D5)]<-'FALSE'
pareto[pareto=='FALSE']<-2
pareto[pareto=='TRUE']<-1

for (i in 5:11){
  pareto[,i]<-as.numeric(pareto[,i])
}
pareto<-round(pareto, 4)

for (i in 1:ncol(pareto)){
  pareto[,i] <- rank(pareto[,i])
}

# Create a matrix with a row and a column for each metric 
# the numbers are the number of tests in which a metric is beaten by each other metric
# If the value is eleven, that metric is beaten on all eleven desirable properties by another metric: pareto-dominated

Pdom<-matrix(nrow=25, ncol=25)
rownames(Pdom)<-rownames(D)
colnames(Pdom)<-rownames(D)

for (i in 1:nrow(pareto)){
  domin<-sapply(1:nrow(pareto), function(j) length(which(pareto[j,]<=pareto[i,])))
  Pdom[i,]<-domin
  Pdom[i,i]<-NA
}
write.csv(Pdom, file="Pdom.csv") # for the supplementary: so we can see which metrics were outperformed on all desirable properties


pd<-sapply(1:nrow(Pdom), function(i) length(which(Pdom[i,]==11))) # sixteen metrics are pareto dominated
names(pd)<-rownames(Pdom) 
# 1 Euclidean
# 2 Av. Euclidean
# 3 Manhattan
# 4 alt.Gower
# 5 CYd
# 6 Binomial
# 7 Lande Shannon
# 8 Lande Simpson
# 9 Classic Sorensen
# 10 Classic Jaccard
# 11 Bray-Curtis
# 12 Morisita-Horn
# 13 Kulczynski
# 14 Renkonen
# 15 Gower
# 16 Jost Simpson


Pdom2<-Pdom[which(pd==0),which(pd==0)] # exclude the pareto dominated metrics
# We are left with 9 metrics that are not outperformed by another metric on all desirable properties

# 1 sim
# 2 Chao Sorensen
# 3 Chao Jaccard
# 4 Jaccard  abd
# 5 Canberra
# 6 Morisita
# 7 Horn
# 8 NESS
# 9 Jost Shannon

# how many of these metrics are beaten on all but one test?
npd<-sapply(1:nrow(Pdom2), function(i) length(which(Pdom2[i,]==10)))
names(npd)<-rownames(Pdom2)
# a further five metrics are nearly pareto dominated: beaten by another metric on all but one test
# 1 Chao Sorensen
# 2 Chao Jaccard 
# 3 Jaccard abd
# 4 Canberra
# 5 Jost Shannon


# how much better are they on this one test?  Is it so much better that it can compensate for being worse on all other tests?
pareto[names(npd),]


# Chao Sorensen only better than Morisita and Horn on D3 (independence of alpha-diversity) and D1 (independence of sample size) respectively: 
pareto['Morisita','D3']-pareto['Chao Sorensen','D3'] # Chao Sorensen is 3.5 ranks better than Morisita on D3
(D['Morisita','D3']/D['Chao Sorensen','D3']) # Chao Sorensen score is 2.6 times better than Morisita on D3
pareto['Horn','D1']-pareto['Chao Sorensen','D1']  # Chao Sorensen is 5.5 ranks better than Horn on D1
(D['Horn','D1']/D['Chao Sorensen','D1'])  # Chao Sorensen score is 2.8 times better than Horn on D1

# Chao Jaccard better than sim, Horn, NESS, Jost Shannon and on D1 only (independence of sample size)
(D['sim','D1']/D['Chao Jaccard','D1']) # 2.3 times better than sim
(D['Horn','D1']/D['Chao Jaccard','D1']) # 2.1 times better than Horn
(D['NESS','D1']/D['Chao Jaccard','D1'])# 2.2 times better than NESS
(D['Jost Shannon','D1']/D['Chao Jaccard','D1'])# 2.0 times better than Jost Shannon
# Chao Jaccard putperformed by Chao Soresen on every test except D5 - multiplicativity
(D['Chao Sorensen','D5']/D['Chao Jaccard','D5']) # 1.2 times better than Chao Sorensen : negligable so discard



# Jaccard abd better than NESS only on D5 (multiplicativity), but the difference is negligable: discard
(D['NESS','D5']/D['Jaccard abd.','D5']) # 1.02 times better


# Canberra only better than sim on D8: because it is an abundance based metric (keep if you need an abundance metric)

D['sim','D8']
D['Canberra','D8'] 

# Jost-Shannon better than Horn on D1 only (independence of sample size)
D['Horn','D1']/D['Jost Shannon','D4'] # 2.66 times better


pareto[c("sim", "Chao Sorensen", "Morisita", "Horn", "NESS"),]




Pdom3<-Pdom2[which(npd==0),which(npd==0)]

# Now we are left with just five metrics

# 1 sim
# 2 Chao Soresen
# 3 Morisita
# 4 Horn
# 5 NESS

D[,1:5]<-round(D[,1:5], 4)
save(D, file="D.rData")
write.csv(D, file="D.csv")
save(P, file="P.rData")
write.csv(P, file="P.csv")
