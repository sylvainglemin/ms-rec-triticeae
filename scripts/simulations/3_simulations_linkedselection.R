# Script to analyse SLiM simulation results
# 3 - Fit of the linked
# Sylvain Glémin (CNRS Rennes, France)
# March 2024

# Reading simulation files


data_sim <- read.table("data/simulations/sigma0.0_rep2_data_20N.txt",header = T,sep = ",")
data_rec <- read.table("data/simulations/sigma0.0_rep2_gene_map.txt",header = T,sep = ",")
data_param <- read.table("data/simulations/sigma0.0_rep2_params.txt",header = F,sep = ":")



# simulated parameters
n_sample <- as.numeric(data_param[3,2])
z <- strsplit(data_param[9,2]," ")[[1]]
p_neutral <- as.numeric(z[3])/(as.numeric(z[3])+as.numeric(z[4]))
theta <- 4*as.numeric(data_param[1,2])*as.numeric(data_param[4,2])*p_neutral

data_sim$start <- 1000*floor(data_sim$position/10^3)
data_sim$end <- 1000*ceiling(data_sim$position/10^3)-1

mydata <- merge(data_sim,data_rec,by=c("start","end"))

mydata$pi <- 2*mydata$count_in_sample*(n_sample-mydata$count_in_sample)/(n_sample*(n_sample-1))


# Test of the phenomenological model ####

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


# Grouping SNPs by windows
Ncat <- 100
data_per_gene$window  <- floor(data_per_gene$start/30000)
  
mydataRec <- data.frame(list("RecCat"=c(1:Ncat)))
mydataRec$rec <- 10^6*aggregate(data_per_gene$rate,by=list(data_per_gene$window),FUN = function(x) mean(x,na.rm=T))$x
mydataRec$piS <- aggregate(data_per_gene$piS,by=list(data_per_gene$window),FUN = function(x) mean(x,na.rm=T))$x
mydataRec$piN <- aggregate(data_per_gene$piN,by=list(data_per_gene$window),FUN = function(x) mean(x,na.rm=T))$x
mydataRec$nbSNPs <- aggregate(data_per_gene$nbSNPs,by=list(data_per_gene$window),FUN = function(x) sum(x,na.rm=T))$x
mydataRec$start <- aggregate(data_per_gene$start,by=list(data_per_gene$window),FUN = function(x) sum(x,na.rm=T))$x


plot(mydataRec$start,mydataRec$piS)

