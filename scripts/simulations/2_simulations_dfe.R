# Script to analyse SLiM simulation results
# Estimation of the DFE
# Sylvain Glémin (CNRS Rennes, France)
# March 2024

library(expint)
library(ggplot2)


# 1 - Generating the file with DFE results (not needed to redoing it each time) ####
source("/Users/sglemin/Documents/Boulot/Recherche/Thematiques/DFE/FastDFE/fast_dfe.R")

# Function to generate SFS
generate_sfs <- function(i) {
  nS <- dim(df_sfs[df_sfs$selection_coeff==0 & df_sfs$count_in_sample==i,])[1]
  nNS <- dim(df_sfs[df_sfs$selection_coeff<0 & df_sfs$count_in_sample==i,])[1]
  #To avoid division by 0
  if(nS==0) nS <- nS + 1
  if(nNS==0) nNS <- nNS + 0.1
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
    pi_obs_tot <-  mean(data_per_gene$piS)
    Ne_obs_tot <- pi_obs_tot/(4*mu)
    Smean_obs_tot <- 4*Ne_obs_tot*(h + Fis - h*Fis)*smean
    Sbound <- 10 # Proportion: 0 < S < 10
    f0_obs_tot <- 1 - gammainc(estim_tot$shape,Sbound*estim_tot$shape/estim_tot$Sdel)/gamma(estim_tot$shape)
    f0_pred_tot <- 1 - gammainc(shape,Sbound*shape/Smean)/gamma(shape)
    f0_adj_tot <- 1 - gammainc(shape,Sbound*shape/Smean_obs_tot)/gamma(shape)
    omega_tot <- estim_tot$p*estim_tot$Sadv
    
    # High recombination
    df_sfs <- mydata[mydata$rate > quantile(mydata$rate,0.5),]
    sfs <- sapply(c(1:(n-1)),generate_sfs)
    estim_high <- least_square_tot(qns = 1/2,syn = sfs[1,],nonsyn = sfs[2,])
    pi_obs_high <- mean(data_per_gene[data_per_gene$rate > quantile(data_per_gene$rate,0.5),]$piS)
    Ne_obs_high <- pi_obs_high/(4*mu)
    Smean_obs_high <- 4*Ne_obs_high*(h + Fis - h*Fis)*smean
    Sbound <- 10 # Proportion: 0 < S < 10
    f0_obs_high <- 1 - gammainc(estim_high$shape,Sbound*estim_high$shape/estim_high$Sdel)/gamma(estim_high$shape)
    f0_pred_high <- 1 - gammainc(shape,Sbound*shape/Smean)/gamma(shape)
    f0_adj_high <- 1 - gammainc(shape,Sbound*shape/Smean_obs_high)/gamma(shape)
    omega_high <- estim_high$p*estim_high$Sadv
    
    # Low recombination
    df_sfs <- mydata[mydata$rate <= quantile(mydata$rate,0.5),]
    sfs <- sapply(c(1:(n-1)),generate_sfs)
    estim_low <- least_square_tot(qns = 1/2,syn = sfs[1,],nonsyn = sfs[2,])
    pi_obs_low <- mean(data_per_gene[data_per_gene$rate <= quantile(data_per_gene$rate,0.5),]$piS)
    Ne_obs_low <- pi_obs_low/(4*mu)
    Smean_obs_low <- 4*Ne_obs_low*(h + Fis - h*Fis)*smean
    Sbound <- 10 # Proportion: 0 < S < 10
    f0_obs_low <- 1 - gammainc(estim_low$shape,Sbound*estim_low$shape/estim_low$Sdel)/gamma(estim_low$shape)
    f0_pred_low <- 1 - gammainc(shape,Sbound*shape/Smean)/gamma(shape)
    f0_adj_low <- 1 - gammainc(shape,Sbound*shape/Smean_obs_low)/gamma(shape)
    omega_low <- estim_low$p*estim_low$Sadv
    estim_low <- least_square_del(qns = 2,syn = sfs[1,],nonsyn = sfs[2,])
    
    res_tot <- c(s,j,"tot",pi_obs_tot,Ne_obs_tot,Smean_obs_tot,f0_obs_tot,f0_pred_tot,f0_adj_tot,omega_tot)
    res_high <- c(s,j,"high",pi_obs_high,Ne_obs_high,Smean_obs_high,f0_obs_high,f0_pred_high,f0_adj_high,omega_high)
    res_low <- c(s,j,"low",pi_obs_low,Ne_obs_low,Smean_obs_low,f0_obs_low,f0_pred_low,f0_adj_low,omega_low)
    result <- rbind(result,res_tot,res_high,res_low)
  }
  print(i)
}

result <- data.frame(result)
names(result) <- c("self","rep","rec","pi","Ne","Smean","f0_obs","f0_pred","f0_adj","omega")
write.table(result,"outputs/simulations/simul_dfe.txt",row.names = F)


# 2 - Analyzing the results ####

# Function to merge multiple plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mydata <- read.table("outputs/simulations/simul_dfe.txt",header = T)

G1 <- ggplot(data = mydata, aes(x = as.factor(self),y = f0_obs,col=rec)) +
  geom_boxplot(outliers = F) +
  scale_color_discrete(name = "Recombination") +
  geom_boxplot(aes(x = as.factor(self),y = f0_adj,fill=rec),col="black",outliers = F) +
  scale_fill_discrete(name = "Recombination") +
  xlab("Selfing rate") +
  ylab(expression("% of weakly deleterious mutations (-10 < 4"*italic(N[e])*italic(s)*" < 0)"))

G2 <- ggplot(data = mydata, aes(x = as.factor(self),y = omega,col=rec)) +
  geom_boxplot(outliers = F) +
  scale_color_discrete(name = "Recombination") +
  xlab("Selfing rate") +
  ylab(expression("Adaptive substitution rate ("*omega[a]*")"))

pdf("figures/sup_mat/Simulations_DFE_Recombination.pdf",width = 8,height = 12)
multiplot(G1,G2)
dev.off()
