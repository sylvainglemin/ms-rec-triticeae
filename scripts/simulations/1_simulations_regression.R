# Script to analyse SLiM simulation result
# 1 Regression (phenomenological model)
# Sylvain Glémin (CNRS Rennes, France)
# March 2024

library(ggplot2)

# Number of recombination categories
Ncat <- 50

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
  theta <- 4*as.numeric(data_param[1,2])*as.numeric(data_param[4,2])*p_neutral/(1+Fis)

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
    
    # Grouping SNPs by recombination categories
    data_per_gene$CatRec  <- cut(data_per_gene$rate,breaks = unique(quantile(data_per_gene$rate, c(0:Ncat)/Ncat,na.rm=T)))
    mydataRec <- data.frame(list("RecCat"=c(1:Ncat)))
    mydataRec$rec <- 10^6*aggregate(data_per_gene$rate,by=list(data_per_gene$CatRec),FUN = function(x) mean(x,na.rm=T))$x
    mydataRec$piS <- aggregate(data_per_gene$piS,by=list(data_per_gene$CatRec),FUN = function(x) mean(x,na.rm=T))$x
    mydataRec$piN <- aggregate(data_per_gene$piN,by=list(data_per_gene$CatRec),FUN = function(x) mean(x,na.rm=T))$x
    mydataRec$nbSNPs <- aggregate(data_per_gene$nbSNPs,by=list(data_per_gene$CatRec),FUN = function(x) sum(x,na.rm=T))$x
    
    # Fitting the phenomenological model
    x <- mydataRec$rec
    y <- mydataRec$piS
    K <- 1
    repeat{
      ainit <- max(mean(mydataRec$piS),rnorm(1,mean(mydataRec$piS)*2,mean(mydataRec$piS)))
      binit <- max(0.5,rnorm(1,2,2))
      blim <- max(5,2*binit)
      cinit <- rnorm(1,0.02,0.005)
      w <- mydataRec$nbSNPs /(y*(1-y))
      w <- w/mean(w)
      ZERO <- 0.00001
      DONE <- F
      tryCatch({
        model <- 
          nls(y ~ a/(1+exp(b-(4*c*x/a))),weights = w,start = list(a=ainit,b=binit,c=cinit),algorithm = "port",lower=c(ZERO,-blim,ZERO),upper=c(1,blim,100),
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024, printEval = TRUE, warnOnly = FALSE) )
        DONE <- T
      },error=function(e){})
      K <- K+1
      if(DONE | K > 1000) break
    }
    if(exists("model")){
      model_sum <- summary(model)
      aestim <- model_sum$parameters[1,1]
      bestim <- model_sum$parameters[2,1]
      cestim <- model_sum$parameters[3,1]
      result <- rbind(result,c(s,theta,mean(mydataRec$piS),aestim))
      rm(model)
    }
  }
  print(tabself[i])
}

result <- data.frame(result)
names(result) <- c("self","theta","piS","piMax")

z <- c(0:1,0.1)
datapred <- data.frame("x"=z,"y"=theta/(1+z/(2-z)))
ggplot(data=result) + geom_point(aes(x = self,y=piS)) + geom_point(aes(x = self,y=piMax,col="red")) +
  geom_line(data=datapred,aes(x,y,col="red")) + scale_y_log10()
