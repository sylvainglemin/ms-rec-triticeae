# Script to analyse SLiM simulation results
# Estimation of the DFE
# Sylvain Glémin (CNRS Rennes, France)
# March 2024

library(expint)

source("/Users/sglemin/Documents/Boulot/Recherche/Thematiques/DFE/FastDFE/fast_dfe.R")

# Function to generate SFS
generate_sfs <- function(i) {
  nS <- dim(df_sfs[df_sfs$selection_coeff==0 & df_sfs$count_in_sample==i,])[1]
  nNS <- dim(df_sfs[df_sfs$selection_coeff<0 & df_sfs$count_in_sample==i,])[1]
  # To avoid division by 0
  if(nS*nNS==0) {
    nS <- nS + 1
    nNS <- nNS + 1
  }
  return(c(nS,nNS))
}


# Reading simulation files

tabself <- c("0.0","0.5","0.9","0.99","1.0")

result <- c()

for(i in c(1:length(tabself))) {
  
  self <- tabself[i]
  data_rec <- read.table(paste("data/simulations/sigma",self,"_rep1_gene_map.txt",sep=""),header = T,sep = ",")
  data_param <- read.table(paste("data/simulations/sigma",self,"_rep1_params.txt",sep=""),header = F,sep = ":")
  # simulated parameters
  n_sample <- as.numeric(data_param[3,2])
  s <- as.numeric(self)
  Fis <- s/(2-s)
  z <- strsplit(data_param[9,2]," ")[[1]]
  p_neutral <- as.numeric(z[3])/(as.numeric(z[3])+as.numeric(z[4]))
  # Simulated theta
  Npop <- as.numeric(data_param[1,2])
  Ne <- Npop/(1+Fis)
  mu <- as.numeric(data_param[4,2])
  h <- 0.25
  smean <- -as.numeric(data_param[7,2])
  shape <- as.numeric(data_param[8,2])
  Smean <- 4*Ne*(h + Fis - h*Fis)*smean
  theta <- 4*Ne*mu*p_neutral
  n <- as.numeric(data_param[3,2])
  
  for(j in c(1:10)){
    
    data_sim <- read.table(paste("data/simulations/sigma",self,"_rep",j,"_data_10N.txt",sep=""),header = T,sep = ",")
    # Preparing data files
    data_sim$start <- 1000*floor(data_sim$position/10^3)
    data_sim$end <- 1000*ceiling(data_sim$position/10^3)-1
    mydata <- merge(data_sim,data_rec,by=c("start","end"))
    mydata$pi <- 2*mydata$count_in_sample*(n_sample-mydata$count_in_sample)/(n_sample*(n_sample-1))
    
    # Grouping SNPs by genes
    rec <- aggregate(mydata[mydata$selection_coeff==0,]$rate,by=list(mydata[mydata$selection_coeff==0,]$start),FUN = function(x) mean(x,na.rm=T))
    piS <- aggregate(mydata[mydata$selection_coeff==0,]$pi,by=list(mydata[mydata$selection_coeff==0,]$start),FUN = function(x) sum(x,na.rm=T)/1000)
    piN <- aggregate(mydata[mydata$selection_coeff<0,]$pi,by=list(mydata[mydata$selection_coeff<0,]$start),FUN = function(x) sum(x,na.rm=T)/1000)
    nbSNPs <- aggregate(mydata[mydata$selection_coeff==0,]$count_in_sample,by=list(mydata[mydata$selection_coeff==0,]$start),FUN = function(x) length(x))
    data_per_gene <- merge(data_rec,piS,by.x="start",by.y="Group.1",all = T)
    names(data_per_gene)[4] <- "piS"
    data_per_gene <- merge(data_per_gene,piN,by.x="start",by.y="Group.1",all = T)
    names(data_per_gene)[5] <- "piN"
    data_per_gene <- merge(data_per_gene,nbSNPs,by.x="start",by.y="Group.1",all = T)
    names(data_per_gene)[6] <- "nbSNPs"
    data_per_gene$piS <- ifelse(is.na(data_per_gene$piS),0,data_per_gene$piS)
    data_per_gene$piN <- ifelse(is.na(data_per_gene$piN),0,data_per_gene$piN)
    data_per_gene$nbSNPs <- ifelse(is.na(data_per_gene$nbSNPs),0,data_per_gene$nbSNPs)
    
    # Estimating the DFE
    
    # Whole dataset
    df_sfs <- mydata
    sfs <- sapply(c(1:(n-1)),generate_sfs)
    estim_tot <- least_square_tot(qns = 1/2,syn = sfs[1,],nonsyn = sfs[2,])
    pi_tot_obs <-  mean(data_per_gene$piS)
    Ne_obs <- pi_tot_obs/(4*mu)
    Smean_obs <- 4*Ne_obs*(h + Fis - h*Fis)*smean
    Sbound <- 10 # Proportion: 0 < S < 10
    f0_tot_obs <- 1 - gammainc(estim_tot$shape,Sbound*estim_tot$shape/estim_tot$Sdel)/gamma(estim_tot$shape)
    f0_tot_pred <- 1 - gammainc(shape,Sbound*shape/Smean)/gamma(shape)
    f0_tot_adj <- 1 - gammainc(shape,Sbound*shape/Smean_obs)/gamma(shape)
    omega_tot <- estim_tot$p*estim_tot$Sadv
    
    # High recombination
    df_sfs <- mydata[mydata$rate > quantile(mydata$rate,0.5),]
    sfs <- sapply(c(1:(n-1)),generate_sfs)
    estim_high <- least_square_tot(qns = 1/2,syn = sfs[1,],nonsyn = sfs[2,])
    pi_obs_high <- mean(data_per_gene[data_per_gene$rate > quantile(data_per_gene$rate,0.5),]$piS)
    Ne_obs <- pi_tot_obs/(4*mu)
    Smean_obs <- 4*Ne_obs*(h + Fis - h*Fis)*smean
    Sbound <- 10 # Proportion: 0 < S < 10
    f0_tot_obs <- 1 - gammainc(estim_tot$shape,Sbound*estim_tot$shape/estim_tot$Sdel)/gamma(estim_tot$shape)
    f0_tot_pred <- 1 - gammainc(shape,Sbound*shape/Smean)/gamma(shape)
    f0_tot_adj <- 1 - gammainc(shape,Sbound*shape/Smean_obs)/gamma(shape)
    omega_tot <- estim_tot$p*estim_tot$Sadv
    estim_high <- least_square_del(qns = 2,syn = sfs[1,],nonsyn = sfs[2,])
    
    
    
    
    
    1 - gammainc(estim$shape,estim$shape/estim$Smean)/gamma(estim$shape)
    1 - gammainc(0.5,0.5/Se)/gamma(0.5)
    
    
  }
}




























# Test of the estimation of the DFE ####

source("/Users/sglemin/Documents/Boulot/Recherche/Thematiques/DFE/FastDFE/fast_dfe.R")


df_sfs <- mydata[mydata$rate > quantile(mydata$rate,0.25),]
pi <- mean(data_per_gene[data_per_gene$rate > quantile(data_per_gene$rate,0.25),]$piS)
ratio <- (0.04/3)/pi
Ne <- pi /(4*10^(-6))
h <- 0.5
s <- 0.01
F <- 0
Se <- 4*Ne*(h+F-h*F)*s
Se

generate_sfs <- function(i) {
  nS <- dim(df_sfs[df_sfs$selection_coeff==0 & df_sfs$count_in_sample==i,])[1]
  nNS <- dim(df_sfs[df_sfs$selection_coeff<0 & df_sfs$count_in_sample==i,])[1]
  return(c(nS,nNS))
}
sfs <- sapply(c(1:14),generate_sfs)
estim <- least_square_del(qns = 2,syn = sfs[1,],nonsyn = sfs[2,])

1 - gammainc(estim$shape,estim$shape/estim$Smean)/gamma(estim$shape)
1 - gammainc(0.5,0.5/Se)/gamma(0.5)
