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

# ! denotes simulations that are not scored as part of the analaysis.

library(ggplot2) 
library(reshape2) 
library(RCurl)

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


######################################################################################################################################
# Simulations to test for each of the properties described above.
# Each metric is given a score to quantify how well it satisfies each property.

##################### D1) unbiased by sample size ######################################################################################

turnover<-(0:5*2)/10
N<-c(10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, sum(Comm1))

by.N<-list()

for (j in 1:length(turnover)){
  z<-list()
  
  for(k in 1:length(N)){
    sims<-replicate(10000, sampling_effect(Comm1, N=N[k], p=turnover[j]), simplify=F)
    Betas<-sapply(1:length(sims), function(s) beta_metrics_all(sims[[s]]))
    Beta<-sapply(1:nrow(Betas), function(i) median(Betas[i,], na.rm=T))
    upper<-sapply(1:nrow(Betas), function(i) quantile(Betas[i,], probs=0.75, na.rm=T))
    lower<-sapply(1:nrow(Betas), function(i) quantile(Betas[i,], probs=0.25, na.rm=T))
    Metric<-rownames(Betas)
    ID<-metrics
    sp_turnover=turnover[j]
    p_sample<-N[k]/sum(Comm1)
    
    z[[k]]<-data.frame(Metric, ID, sp_turnover, p_sample, Beta, lower, upper, row.names=NULL)
  }
  a<-do.call(rbind,z)
  by.N[[j]]<-a
}

D1_sim<-do.call(rbind, by.N)



# score the performance  
x<-acast(D1_sim, sp_turnover ~ p_sample ~ ID, value.var="Beta")
D1_score<-sapply(1:25, function(i) mean(abs(x[,,i]-x[,16,i])/(max(x[,,i])-min(x[,,i]))))
names(D1_score)<-names(sort(metrics))

# add the scores into the data
for (i in 1:length(D1_score)){
  D1_sim[D1_sim$Metric==names(D1_score)[i],"D1_score"]<-paste(names(D1_score)[i], "\nD1 =", round(D1_score[i],4))
}
D1_sim$D1_score<-factor(D1_sim$D1_score, levels=sapply(1:length(sort(D1_score)), function(i) paste(names(sort(D1_score))[i], "\nD1 =", round(sort(D1_score)[i],4))))

# plot the response to undersampling both assemblages
a<-ggplot(D1_sim[D1_sim$ID>8,],aes(x=p_sample,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D1_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.1, 1.1), breaks=((0:5)*2/10))+scale_x_reverse(limits=c(1, 0))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
f<-d+ylab (expression(beta))+xlab("Proportion of individuals sampled")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D1_sim$p_sample), colour="black", linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
j

a<-ggplot(D1_sim[D1_sim$ID<=8,],aes(x=p_sample,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D1_score, scales="free_y", nrow=3, ncol=3)+scale_x_reverse(limits=c(1, 0))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))+theme(axis.title.y=element_text(angle=90))
f<-d+ylab (expression(beta))+xlab("Proportion of individuals sampled")+labs(fill="Proportion of\nspecies turned over", colour="Proportion of\nspecies turned over")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D1_sim$p_sample), colour="black", linetype=2)+theme(strip.text=element_text(size=8))
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
j
###########################################################################################################################################################################################################################

##################### D2) unbiased by unequal sample size ##################################################################################################################################################################

turnover<-(0:5*2)/10
N<-c(10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, sum(Comm1))
by.N<-list()

for (j in 1:length(turnover)){
  z<-list()
  
  for(k in 1:length(N)){
    p.turnover<-turnover[j]
    N.sample<-N[k]
    sim.1000<-replicate(10000, beta_metrics_all(unequal_sampling(Comm1, N=N[k], p=turnover[j])))
    Beta<-sapply(1:nrow(sim.1000), function(s) median(sim.1000[s,], na.rm=T))
    upper<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.75, na.rm=T)) # 25th percentile
    lower<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.25, na.rm=T)) # 75th percentile
    Metric<-rownames(sim.1000)
    ID<-metrics
    Difference<-sum(Comm1)-N[k]
    z[[k]]<-data.frame(Metric, ID, N.sample, Difference, p.turnover, Beta, upper, lower, row.names=NULL)
  }
  a<-do.call(rbind,z)
  by.N[[j]]<-a
}
D2_sims<-do.call(rbind, by.N)

# score the performance 
x<-acast(D2_sims, p.turnover ~ N.sample ~ ID, value.var="Beta")
D2_score<-sapply(1:25, function(i) mean(abs(x[,,i]-x[,16,i])/(max(x[,,i])-min(x[,,i]))))
names(D2_score)<-names(sort(metrics))

# add the scores into the data
for (i in 1:length(D2_score)){
  D2_sims[D2_sims$Metric==names(D2_score)[i],"D2_score"]<-paste(names(D2_score)[i], "\nD2 =", round(D2_score[i],4))
}
D2_sims$D2_score<-factor(D2_sims$D2_score, levels=sapply(1:length(sort(D2_score)), function(i) paste(names(sort(D2_score))[i], "\nD2 =", round(sort(D2_score)[i],4))))


#  plot the response to undersampling one assemblage

options(scipen=20)
a<-ggplot(D2_sims[D2_sims$ID>8,],aes(x=Difference,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D2_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.1), breaks=((0:5)*2/10))+scale_x_continuous(limits=c(0, 10000))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab("Sample size difference")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, colour="black", linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
j


a<-ggplot(D2_sims[D2_sims$ID<=8,],aes(x=Difference,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D2_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 10000))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))+theme(strip.text=element_text(size=8))
f<-e+ylab (expression(beta))+xlab ("Sample size difference")+labs(fill="Species turnover", colour="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, colour="black", linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
j

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
for (j in 1:length(turnover)){
  z<-list()
  
  for(k in 1:length(S)){
    p.turnover<-rep(turnover[j], length(metrics))
    Alpha<-rep(fishers.alpha(N=10000, S=S[k]), length(metrics))
    sim.1000<-replicate(10000, beta_metrics_all(random_composition_change(Comm1=alphas[[k]], p=turnover[[j]])))
    Beta<-sapply(1:nrow(sim.1000), function(s) median(sim.1000[s,]))
    upper<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.75)) # 25th percentile
    lower<-sapply(1:nrow(sim.1000), function(s) quantile(sim.1000[s,], probs=0.25)) # 75th percentile
    Metric<-rownames(sim.1000)
    ID<-metrics
    z[[k]]<-data.frame(Metric, ID, Alpha, p.turnover, Beta, upper, lower, row.names=NULL)
  }
  a<-do.call(rbind,z)
  by.alpha[[j]]<-a
}

D3_sim<-do.call(rbind, by.alpha)

# score the performances

x<-acast(D3_sim, p.turnover ~ Alpha ~ ID, value.var='Beta')
D3_score<-sapply(1:25, function(i) mean(abs(x[,10,i]-x[,,i])/(max(x[,,i])-min(x[,,i]))))
names(D3_score)<-names(sort(metrics))


# add the scores in to the D3_sim data.frame
for (i in 1:length(D3_score)){
  D3_sim[D3_sim$Metric==names(D3_score)[i],"D3_score"]<-paste(names(D3_score)[i], "\nD3 =", round(D3_score[i],4))
}
D3_sim$D3_score<-factor(D3_sim$D3_score, levels=sapply(1:length(sort(D3_score)), function(i) paste(names(sort(D3_score))[i], "\nD3 =", round(sort(D3_score)[i],4))))


##### Plot the response to equal alpha-diversity in assemblage pairs 


a<-ggplot(D3_sim[D3_sim$ID>8,],aes(x=Alpha,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D3_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.01), breaks=((0:6)*2/10))+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(alpha[Fisher]))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=max(D3_sim$Alpha), linetype=2, colour="black")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
j

a<-ggplot(D3_sim[D3_sim$ID<=8,],aes(x=Alpha,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~D3_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(alpha[Fisher]))+labs(fill="Species turnover", colour="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")+theme(axis.title.y=element_text(angle=0))
h<-g+geom_vline(xintercept=max(D3_sim$Alpha), linetype=2, colour="black")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j
#################################################################################################################################################################################################################################################################


########################### D4) additive under directional turnover ########################################################

# generate a hypothetical environmental gradient across three assemblages
# species shared between assemblages 1 and 2 are more likely to be absent in assemblage three
# the probability of being turned over between assemblages 2 and 3 is much higher for species that were present in assemblage 1

grad<-c(1, 5, 10, 50, 100, 500, 1000) # how many times more likely are assemblage 1 species to be turned over between assemblages 2 and 3
turnover<-0:5/10
addtu<-list()
for (k in 1:length(grad)){ 
  add<-list()
  for (i in 1:length(turnover)){
    x<-replicate(10000, add_comp_change(Comm1, p=turnover[i], gradient=grad[k]), simplify=FALSE) 
    # gradient = Species shared between Comm1 and Comm2 are gradient times more likely to be turned over in Comm3
    # gradient specifies the steepness of the spatial gradient - the severity of spatial turnover
    B1_2<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(1,2),])))
    B2_3<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(2,3),])))
    B1_3<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(1,3),])))
    Metric<-B1_2$Var1
    ID<-metrics
    P.turnover<-turnover[i]
    gr<-grad[k]
    Beta1_2<-B1_2$value
    Beta2_3<-B2_3$value
    Beta1_3<-B1_3$value
    Beta.add<-Beta1_2+Beta2_3
    add[[i]]<-data.frame(Metric, ID, gr, P.turnover, Beta1_2, Beta2_3, Beta1_3, Beta.add, row.names=NULL)
  }
  addtu[[k]]<-do.call(rbind, add)
}
D4_sim<-do.call(rbind, addtu)

for (i in 1:25){
  D4_sim[D4_sim$ID==i, "range"]<-max(D4_sim[D4_sim$ID==i,"Beta.add"], na.rm=T)-min(D4_sim[D4_sim$ID==i,"Beta.add"],na.rm=T)
}

D4_sim$diff<-(D4_sim$Beta.add-D4_sim$Beta1_3)/D4_sim$range # difference standardised by the range
qn<-function(x){
  quantile(x, probs=c(0.25, 0.5, 0.75))
}
summary_add<-aggregate(diff~Metric+ID+gr+P.turnover, data=D4_sim, FUN=qn)
D4_summary<-data.frame(summary_add[,1:4], LQ=summary_add$diff[,'25%'], median=summary_add$diff[,'50%'], UQ=summary_add$diff[,'75%'] )

D4_score<-sapply(1:25, function(i) mean(abs(D4_sim[D4_sim$ID==i & D4_sim$gr>1, 'diff'])))
names(D4_score)<-names(sort(metrics))

# add the scores in to the data.frame D4_summary
for (i in 1:length(D4_score)){
  D4_summary[D4_summary$Metric==names(D4_score)[i],"D4_score"]<-paste(names(D4_score)[i], "\nD4 =", round(D4_score[i], 4)) 
}
D4_summary$D4_score<-factor(D4_summary$D4_score, levels=sapply(1:length(sort(D4_score)), function(i) paste(names(sort(D4_score))[i], "\nD4 =", round(sort(D4_score)[i], 4))))

#plot the difference between added beta (B_12+B23)  and observed beta (B_13)
a<-ggplot(D4_summary[D4_summary$ID>8,],aes(x=gr,y=median, ymin=LQ, ymax=UQ, group=factor(P.turnover), fill=factor(P.turnover), colour=factor(P.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D4_score, nrow=6, ncol=3)+scale_x_log10(breaks=c(1,10,100,1000))+geom_hline(yintercept=0, linetype=2)+scale_y_continuous(limits=c(-0.4, 1), breaks=-2:5*2/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 8, colour = 'black'))+theme(axis.title.y =element_text(angle=0,size = 8, colour = 'black'))
f<-e+ylab (expression((beta["1,2"]+beta["2,3"])-beta["1,3"]))+xlab ("gradient")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))
h

a<-ggplot(D4_summary[D4_summary$ID<=8,],aes(x=gr,y=median, ymin=LQ, ymax=UQ, group=factor(P.turnover), fill=factor(P.turnover), colour=factor(P.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D4_score, nrow=3, ncol=3)+scale_x_log10(breaks=c(1,10,100,1000))+geom_hline(yintercept=0, linetype=2)+scale_y_continuous(limits=c(-0.4, 1),breaks=-2:5*2/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 8, colour = 'black'))+theme(axis.title.y =element_text(angle=0,size = 8, colour = 'black'))
f<-e+ylab (expression((beta["1,2"]+beta["2,3"])-beta["1,3"]))+xlab ("gradient")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(axis.title.y=element_text(angle=0))+theme(strip.text=element_text(size=8))
h
###########################################################################################################################################################################################


######################## D5) multiplicative (probabilistic) under orthogonal sampling ######################################################################################################

turnover<-(0:5)*2/10

mult<-list()
for (i in 1:length(turnover)){
  x<-replicate(10000, mult_comp_change(Comm1=Comm1, p=turnover[i]), simplify=F) 
  B1_2<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(1:2),])))
  B2_3<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(2,3),])))
  B1_3<-melt(sapply(1: length(x), function(j) beta_metrics_all(x[[j]][c(1,3),])))
  Metric<-B1_2$Var1
  ID<-metrics
  P.turnover<-1-turnover[i]
  Sim1_2<-1-B1_2$value # the similarity complement of beta is expected to be multiplicative
  Sim2_3<-1-B2_3$value
  Sim1_3<-1-B1_3$value
  Sim.mult<-(Sim1_2*Sim2_3)
  y<-data.frame(Metric, ID, P.turnover, Sim1_2, Sim2_3, Sim1_3, Sim.mult, row.names=NULL)
  mult[[i]]<-y
}

D5_sim<-do.call(rbind, mult)
# It is not possible to test whether similarity is multiplicative for unstandardised indices: they do not have a similraity complement
# so the unstandardised metrics are removed from the data set (this counts as a fail for this property)
D5_sim<-D5_sim[D5_sim$ID>8,]

for (i in 9:25){
  D5_sim[D5_sim$ID==i, "range"]<-max(D5_sim[D5_sim$ID==i,"Sim.mult"])-min(D5_sim[D5_sim$ID==i,"Sim.mult"])
}


D5_sim$diff<-(D5_sim$Sim.mult-D5_sim$Sim1_3)/D5_sim$range # difference standardised by the range

qn<-function(x){
  quantile(x, probs=c(0.25, 0.5, 0.75))
}
summary_mult<-aggregate(diff~Metric+ID+P.turnover, data=D5_sim, FUN=qn)
D5_summary<-data.frame(summary_mult[,1:3], LQ=summary_mult$diff[,'25%'], median=summary_mult$diff[,'50%'], UQ=summary_mult$diff[,'75%'] )


D5_score<-sapply(9:25, function(i) mean(abs(D5_sim[D5_sim$ID==i, 'diff'])))
D5_score<-c(rep(NA, 8), D5_score)
names(D5_score)<-names(sort(metrics))

# add the scores in to the data.frame D4_summary
for (i in 9:25){
  D5_summary[D5_summary$Metric==names(D5_score)[i],"D5_score"]<-paste(names(D5_score)[i], "\nD5 =", round(D5_score[i], 4)) 
}
D5_summary$D5_score<-factor(D5_summary$D5_score, levels=sapply(1:length(sort(D5_score)), function(i) paste(names(sort(D5_score))[i], "\nD5 =", round(sort(D5_score)[i], 4))))

# plot the difference between the value of beta expected under multiplicative behaviour and the observed beta
a<-ggplot(D5_summary,aes(x=P.turnover,y=median, ymin=LQ, ymax=UQ))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~D5_score,scales="fixed", nrow=6, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))+geom_hline(linetype=2, yintercept=0)
d<-c+theme(axis.text.x=element_text(angle=0, hjust=1, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression((1-beta["1,2"])(1-beta["2,3"])-(1-beta["1,3"])))+xlab ("Proportion turnover")
f

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
    Metric<-names(Beta)
    ID<-metrics
    sp_turnover=p[j]
    sp_loss=S[1]-S[i]
    nestedness[[i]]<-data.frame(Metric, ID, sp_turnover, sp_loss, Beta, upper, lower, row.names=NULL)   
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

# add the scores in to the data.frame D4_summary
for (i in 1:length(D6_score)){
  D6_sim[D6_sim$Metric==names(D6_score)[i],"D6_score"]<-paste("D6 =", D6_score[i]) 
}


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
f<-e+ylab (expression(beta))+xlab ("partial correlation")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=1, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j


a<-ggplot(D7_sim[D7_sim$ID<=8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("partial correlation")+labs(fill="Species turnover", colour="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=1, linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j


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



#######################################################################################################################################################################################


############################# D11) Defined minima and maxima ############################################################################################################ 


# see Table 1
D11_score<-as.character(as.logical(c(rep(0, 8), rep(1, 17))))
names(D11_score)<-names(sort(metrics))
##############################################################################################################################################################


################################ P1) Proportionality to species turnover ########################################################################## 


# use simulations in D6_sim (species turnover and species loss)
x<-acast(D6_sim, sp_turnover~sp_loss~ID, value.var="Beta")

P1_score <- sapply(1:25, function(i) mean(sapply(1:dim(x)[2], function(j) summary(lm(x[,j,i]~poly(as.numeric(rownames(x)), 2, raw=T)))$coefficients[3,4])))
names(P1_score)<-names(sort(metrics))


# add the scores to D6_sim data frame
for (i in 1:length(P1_score)){
  D6_sim[D6_sim$Metric==names(P1_score)[i],"P1_score"]<-paste("P1 =", round(P1_score[i], 4)) 
}



##############################################################################################################################################################

############################################## P2)   Sensitivity to differences in alpha diversity ################################



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
    Metric<-rownames(sim.1000)
    ID<-metrics
    z[[i]]<-data.frame(Metric, ID, Alpha_D, p.turnover, Beta, upper, lower, row.names=NULL)
  } # stick these data frames  for each level of alpha diversity difference together
  a<-do.call(rbind,z)
  diff_alpha[[j]]<-a
}

P2_sim<-do.call(rbind, diff_alpha)

# score the perfomance


# score: mean difference between median beta and a reference level of beta diversity under equal alpha diversity (see Table 3)
y<-acast(P2_sim, p.turnover ~ Alpha_D ~ ID, value.var='Beta')
P2_score<-sapply(1:25, function(i) mean(abs(y[,1,i]-y[,,i])/(max(y[,,i])-min(y[,,i]))))
names(P2_score)<-names(sort(metrics))


# add the scores into the data.frame P2_sim
for (i in 1:length(P2_score)){
  P2_sim[P2_sim$Metric==names(P2_score)[i],"P2_score"]<-paste(names(P2_score)[i], "\nP2 =", round(P2_score[i],4))
}
P2_sim$P2_score<-factor(P2_sim$P2_score, levels=sapply(1:length(sort(P2_score)), function(i) paste(names(sort(P2_score))[i], "\nP2 =", round(sort(P2_score)[i],4))))

# plot the responses to differences in alpha-diversity

a<-ggplot(P2_sim[P2_sim$ID>8,],aes(x=Alpha_D,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~P2_score, scales="fixed", nrow=6, ncol=3)+scale_y_continuous(limits=c(-0.0001, 1.01), breaks=((0:6)*2/10))+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste("Difference in ", alpha[Fisher])))+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=min(P2_sim$Alpha_D), linetype=2, colour="black")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

a<-ggplot(P2_sim[P2_sim$ID<=8,],aes(x=Alpha_D,y=Beta, group=factor(p.turnover), colour=factor(p.turnover), fill=factor(p.turnover), ymin=lower,ymax=upper)) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~P2_score, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 60), breaks=(0:6)*10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab (expression(paste("Difference in ", alpha[Fisher])))+labs(fill="Species turnover", colour="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=min(P2_sim$Alpha_D), linetype=2, colour="black")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j
############################################################################################################################################################################################################################################################



################################ P3) Relative sensitivity to nestedness and turnover components of beta #############################################################


# use the simulations in D6_sim

# score  the performmaces
# personality iv) Relative sensitivity to species loss and species turnover
# score: median beta for extreme species loss (s=90) / median beta for complete turnover (t=1)  
x<-acast(D6_sim, sp_turnover ~ sp_loss ~ ID, value.var="Beta")
P3_score<-sapply(1:25, function(i) x[1,10,i]/x[6,1,i])
names(P3_score)<-names(sort(metrics))

# add the scores to the D6_sim data frame
for (i in 1:length(P3_score)){
  D6_sim[D6_sim$Metric==names(P3_score)[i],"P3_score"]<-paste("P3 =", round(P3_score[i],4))
}


D6_sim$title<-paste(D6_sim$Metric, "\n", D6_sim$D6_score, "\n", D6_sim$P1_score, "\n", D6_sim$P3_score, sep="")

a<-ggplot(D6_sim[D6_sim$ID>8 & D6_sim$sp_loss %in% c(0, 10, 30, 50, 70, 90),],aes(x=sp_turnover,y=Beta, ymin=lower, ymax=upper, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=(0:5*2/10))+scale_y_continuous(limits=c(-0.04, 1.101), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("Proportion species turnover")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=0, linetype=2)
j<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

#unstandardised metrics
#x11()
a<-ggplot(D6_sim[D6_sim$ID<=8 & D6_sim$sp_loss %in% c(0, 10, 30, 50, 70, 90),],aes(x=sp_turnover,y=Beta, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss), ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title,scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=(0:5*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ('Proportion species turnover')+labs(colour="Species loss", fill="Species loss")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=0, linetype=2)
j<-g+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j



###############################################################################################################################################################################


################################ P4) Relative sensitivity to decoupling of species ranks and species turnover components of beta ######################################



# use the simulations in D7_sim

# score the performances
# score: value of beta for extreme decoupling of ranks / value of beta for complete species turnover.  
x<-acast(D7_sim, p.turnover ~ partial.cor ~ ID, value.var="Beta")
P4_score<-sapply(1:25, function(i) x[1,1,i]/x[6,length(r),i])
names(P4_score)<-names(sort(metrics))

# add the scores into the data.frame D7_sim
for (i in 1:length(P4_score)){
  D7_sim[D7_sim$Metric==names(P4_score)[i],"P4_score"]<-paste("P4 =", round(P4_score[i],4))
}
D7_sim$title<-paste(D7_sim$Metric, "\n", D7_sim$D7_score, "\n", D7_sim$D9_score, "\n", D7_sim$P4_score, sep="")


a<-ggplot(D7_sim[D7_sim$ID>8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)+scale_y_continuous(limits=c(-0.04, 1.03), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("partial correlation")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=1, linetype=2)
j<-g+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

a<-ggplot(D7_sim[D7_sim$ID<=8,],aes(x=partial.cor,y=Beta,ymin=lower,ymax=upper, group=factor(p.turnover), fill=factor(p.turnover), colour=factor(p.turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA, alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(-1,1), breaks=(-2:2*5)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("partial correlation")+labs(fill="Species turnover", colour="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
#h<-g+geom_vline(xintercept=1, linetype=2)
j<-g+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j

################################ P5) Relative sensitivity to evenness differences and species turnover components of beta ######################################

# use the simulations in D8_sim

# score the performances
# score: value of beta for extreme decoupling of ranks / value of beta for complete species turnover
x<-acast(D8_sim, sp_turnover ~ Even_diff ~ ID, value.var="Beta")
P5_score<-sapply(1:25, function(i) x[1,length(even),i]/x[6,1,i])
names(P5_score)<-names(sort(metrics))

for (i in 1:length(P5_score)){
  D8_sim[D8_sim$Metric==names(P5_score)[i],"P5_score"]<-paste("P5 =", round(P5_score[i],4))
}
D8_sim$title<-paste(D8_sim$Metric, "\n", D8_sim$D8_score, "\n", D8_sim$D10_score, "\n", D8_sim$P5_score, sep="")

a<-ggplot(D8_sim[D8_sim$ID>8,],aes(x=Even_diff,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title, scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0, 1), breaks=0:5*2/10)+scale_y_continuous(limits=c(0, 1), breaks=(0:5*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("Difference in Shannon's Evenness")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

a<-ggplot(D8_sim[D8_sim$ID<=8,],aes(x=Even_diff,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), colour=factor(sp_turnover), fill=factor(sp_turnover))) 
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~title, scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0, 1), breaks=(0:5*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("Difference in Shannon's Evenness")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j

####################################################################################################################################################################################################



############################# P6) Decreased sensitivity to turnover in rare species under a positive ONR #######################################################################################################################

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
  sep.p[[i]]<-data.frame(Metric, ID, p.turnover, Beta, upper, lower)
}

P6_sim<-do.call(rbind, sep.p)
rownames(P6_sim)<-NULL


# add the reference value of beta in the absence of a positive ONR (plus the upper and lower quartiles)
for (i in 1:25){
  P6_sim[P6_sim$ID==i,"lower_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"lower"]
  P6_sim[P6_sim$ID==i,"upper_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"upper"]
  P6_sim[P6_sim$ID==i,"Beta_Test1"]<-D6_sim[D6_sim$ID==i & D6_sim$sp_loss==0,"Beta"]
}

P6_score<-sapply(1:25, function(i) with(P6_sim[P6_sim$ID==i,], mean(abs((Beta_Test1-Beta)/(max(Beta_Test1)-min(Beta_Test1))))))
names(P6_score)<-names(sort(metrics))

# add the scores into the data.frame P6_sim
for (i in 1:length(P6_score)){
  P6_sim[P6_sim$Metric==names(P6_score)[i],"P6_score"]<-paste(names(P6_score)[i], "\nP6 =", round(P6_score[i],4))
}
P6_sim$P6_score<-factor(P6_sim$P6_score, levels=sapply(1:length(sort(P6_score)), function(i) paste(names(sort(P6_score))[i], "\nP6 =", round(sort(P6_score)[i],4))))

# plot the decreased sensitivty to turnover in rare species under a positive ONR

# standardised
a<-ggplot(data=P6_sim[P6_sim$ID>8,],aes(x=p.turnover,y=Beta, ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.3)+facet_wrap(~P6_score,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5)*2)/10)+scale_y_continuous(limits=c(-0.04,1.101), breaks=((0:5)*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ('Proportion species turnover')
g<-f+geom_line(aes(p.turnover, Beta_Test1), linetype=2)+geom_ribbon(aes(ymin=lower_Test1, ymax=upper_Test1), colour=NA, alpha=0.3)
h<-g+theme(legend.position="none")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

# unstandardised
a<-ggplot(data=P6_sim[P6_sim$ID<=8,],aes(x=p.turnover,y=Beta, ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.3)+facet_wrap(~P6_score,scales="free_y", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5)*2)/10)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=6))+theme(axis.text.y=element_text(size=6))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab ('Proportion species turnover')
g<-f+geom_line(aes(p.turnover, Beta_Test1), linetype=2)+geom_ribbon(aes(ymin=lower_Test1, ymax=upper_Test1), colour=NA, alpha=0.3)
h<-g+theme(legend.position="none")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j
###############################################################################################################################################################################################################

############################################## !P7) Sensitivity to species loss (nestedness) ################################

# not included in analysis: This is the presence-absence equivalent of sensitivity to differences in alpha-diversity
# We are interested predominantly in the sensitivity to abundance differences
# use simulations in D6
# 
# score: mean difference between median beta and a reference value of beta under no species loss 
x<-acast(D6_sim, sp_turnover ~ sp_loss ~ ID, value.var="Beta")
P7_score<-sapply(1:25, function(i) mean(abs(x[,,i]-x[,1,i])/(max(x[,,i])-min(x[,,i]))))
names(P7_score)<-names(sort(metrics))

# add the scores into the data.frame D6_sim
for (i in 1:length(P7_score)){
  D6_sim[D6_sim$Metric==names(P7_score)[i],"P7_score"]<-paste(names(P7_score)[i], "\nP7 score =", round(P7_score[i],4))
}
D6_sim$P7_score<-factor(D6_sim$P7_score, levels=sapply(1:length(sort(P7_score)), function(i) paste(names(sort(P7_score))[i], "\nP7 score =", round(sort(P7_score)[i],4))))

# plot the responses to species nestedness

a<-ggplot(D6_sim[D6_sim$ID>8 & D6_sim$sp_turnover %in% c(0:5*2/10),],aes(x=sp_loss,y=Beta, ymin=lower, ymax=upper, group=factor(sp_turnover), fill=factor(sp_turnover), colour=factor(sp_turnover)))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~P7_score,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,100), breaks=(0:5*20))+scale_y_continuous(limits=c(-0.04, 1.101), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ("Species loss")+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j
#unstandardised metrics
#x11()
a<-ggplot(D6_sim[D6_sim$ID<=8 & D6_sim$sp_turnover %in% c(0:5*2/10),],aes(x=sp_loss,y=Beta, group=factor(sp_turnover), fill=factor(sp_turnover), colour=factor(sp_turnover), ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~P7_score,scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0,100), breaks=(0:5*20))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ('Species loss')+labs(colour="Species turnover", fill="Species turnover")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j
############################# !P8)  Sensitivity to species turnover ################################################################################################################################

# use the simulations in D6_sim

# score the performances

# score: mean difference between median beta and a reference value of beta under no turnover 
x<-acast(D6_sim, sp_turnover ~ sp_loss ~ ID, value.var="Beta")
P8_score<-sapply(1:25, function(i) mean(abs(x[,,i]-x[1,,i])/(max(x[,,i])-min(x[,,i]))))
names(P8_score)<-names(sort(metrics))

# add the scores into the data frame
for (i in 1:length(P8_score)){
  D6_sim[D6_sim$Metric==names(P8_score)[i],"P8_score"]<-paste(names(P8_score)[i], "\nP8 score =", round(P8_score[i],4))
}
D6_sim$P8_score<-factor(D6_sim$P8_score, levels=sapply(1:length(sort(P8_score)), function(i) paste(names(sort(P8_score))[i], "\nP8 score =", round(sort(P8_score)[i],4))))


nest<-c(0, 10, 30, 50, 70, 90)
SpTurn<-D6_sim[D6_sim$sp_loss %in% nest,]

a<-ggplot(SpTurn[SpTurn$ID>8,],aes(x=sp_turnover,y=Beta, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss), ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~P8_score,scales="fixed", nrow=6, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5*2)/10))+scale_y_continuous(limits=c(-0.04, 1.101), breaks=((0:5)*2/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0, size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ('Species turnover')+theme(legend.position="none")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j
# unstandardised
a<-ggplot(SpTurn[SpTurn$ID<=8,],aes(x=sp_turnover,y=Beta, group=factor(sp_loss), fill=factor(sp_loss), colour=factor(sp_loss), ymin=lower, ymax=upper))
b<-a+geom_line(size=0.5)+geom_ribbon(colour=NA,alpha=0.5)+facet_wrap(~P8_score,scales="free_y", nrow=3, ncol=3)+scale_x_continuous(limits=c(0,1), breaks=((0:5*2)/10))
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 10, colour = 'black'))+theme(axis.title.y = element_text(angle=0, size = 12, colour = 'black'))
f<-e+ylab (expression(beta))+xlab ('Species turnover')+labs(colour="Species loss", fill="Species loss")
g<-f+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
h<-g+geom_vline(xintercept=0, linetype=2)
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j


## !P9: an extra simulation to help visualise how  metrics use abundance information.

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

P9_sim<-Z


# plot the response to a single species turned over 
a<-ggplot(P9_sim[P9_sim$ID>8,],aes(x=Rel_abd,y=Beta)) 
b<-a+geom_line(size=0.5)+facet_wrap(~Metric, scales="fixed", nrow=6, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))+xlim(0, 0.155)+ylim(0, 0.4)
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab("Relative abundance of single species turned over")
h<-f+theme(legend.position="none")
j<-h+ggtitle("a)")+theme(plot.title=element_text(hjust=0))
j

a<-ggplot(P9_sim[P9_sim$ID<=8,],aes(x=Rel_abd,y=Beta)) 
b<-a+geom_line(size=0.5)+facet_wrap(~Metric, scales="free", nrow=3, ncol=3)
c<-b+ theme_bw()+theme(panel.grid.major = element_line(colour ="white"), panel.grid.minor = element_line(colour ="white"))
d<-c+theme(axis.text.x=element_text(angle=0, size=8))+theme(axis.text.y=element_text(size=8))
e<-d+theme(axis.title.x = element_text(size = 12, colour = 'black'))+theme(axis.title.y = element_text(angle=0,size = 16, colour = 'black'))
f<-e+ylab (expression(italic(beta)))+xlab("Relative abundance of single species turned over")
h<-f+theme(legend.position="none")
j<-h+ggtitle("b)")+theme(plot.title=element_text(hjust=0))+theme(strip.text=element_text(size=8))
j

##################################################################################################################################################################################################################################################

##### Results figures and tables

# Table 3
D<-data.frame(D1=D1_score, D2=D2_score, D3=D3_score, D4=D4_score, D5=D5_score, D6=D6_score, D7=D7_score, D8=D8_score, D9=D9_score, D10=D10_score, D11=D11_score, check.rows=TRUE, stringsAsFactors=F)
D[,1:5]<-round(D[,1:5], 4)

# Table 4
P<-data.frame(P1=P1_score, P2=P2_score, P3=P3_score, P4=P4_score, P5=P5_score, P6=P6_score)#, P7=P7_score, P8=P8_score)
P<-round(P,4)

# Fig. 1
DP<-cbind(D[,1:4], P)
pca1<-prcomp(DP, scale=T)
summary(pca1)
plot(pca1)

# pca axes 1 and 2
par(xpd=NA)
biplot(pca1, las=1, tck=0.02, cex=0.8,  col=c("black", "dark grey"), asp=1.8, cex.axis=0.9)
mtext(outer=T, "a)", line=-6, side=3, adj=0.1)

# pca axis 2 and 3
par(xpd=NA)
biplot(pca1, choices=2:3, las=1, tck=0.02, cex=0.8,  col=c("black", "dark grey"), asp=1.8, cex.axis=0.9)
mtext(outer=T, "b)", line=-6, side=3, adj=0.1)




############################################################################################################################################################################################################################################################
# Are any of the metrics outperformed by at least one other metric on all desirable properties?
# These are Pareto-dominated

# convert all to numeric variables so they can be ranked
# TRUE gets 1, FALSE gets 2
# NAs on test 5 get 2 (FAIL)
pareto<-D
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



pd<-sapply(1:nrow(Pdom), function(i) length(which(Pdom[i,]==11))) # fifteen metrics are pareto dominated
names(pd)<-rownames(Pdom)

Pdom2<-Pdom[which(pd==0),which(pd==0)] # exclude the pareto dominated metrics
# We are left with ten metrics that are not outperformed by another metric on all desirable properties


# how many of these metrics are beaten on all but one test?
npd<-sapply(1:nrow(Pdom2), function(i) length(which(Pdom2[i,]==10)))
names(npd)<-rownames(Pdom2)

# how much better are they on this one test?  Is it so much better that it can compensate for being worse on the other tests?

# Chao Sorensen only better than sim and Horn on D1: 
pareto['sim','D1']-pareto['Chao Sorensen','D1'] # how many ranks better is Chao Sorensen than sim?
(D['sim','D1']/D['Chao Sorensen','D1']) # how many times better is Chao Sorensen score than sim?
pareto['Horn','D1']-pareto['Chao Sorensen','D1'] 
(D['Horn','D1']/D['Chao Sorensen','D1'])

# Chao Jaccard only better than sim on D1 : 
pareto['sim','D1']-pareto['Chao Jaccard','D1'] 
(D['sim','D1']/D['Chao Jaccard','D1'])

# Canberra only better than sim on D8: because it is an abundance based metric

D['sim','D8']
D['Canberra','D8'] 

# morisita-Horn only better than morisita on D4
pareto['Morisita','D4']-pareto['Morisita-Horn','D4'] 
D['Morisita','D4']/D['Morisita-Horn','D4']

# Jaccard abd only better than NESS on D5 : 
pareto['NESS','D5']-pareto['Jaccard abd.','D5'] 
D['NESS','D5']/D['Jaccard abd.','D5']

# Jost Shannon beats Horn only on D1: 
pareto['Horn','D1']-pareto['Jost Shannon','D1'] 
(D['Horn','D1']/D['Jost Shannon','D1'])





Pdom3<-Pdom2[which(npd==0),which(npd==0)]

# Now we are left with just four metrics: Morisita, sim, Horn and NESS 


