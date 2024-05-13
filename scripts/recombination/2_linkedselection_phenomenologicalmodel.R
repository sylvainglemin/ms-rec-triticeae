# Script to fit a phenomenological linked-selection model to the data
# Provides the estimate of the slope of piS as a function of recombination and piMax
# Also provides some graphs (not presented in the article)
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022


# packages to install
if (!require("plotrix"))   install.packages("plotrix", dependencies = TRUE)
if (!require("DescTools"))   install.packages("DescTools", dependencies = TRUE)
if (!require("PropCIs"))   install.packages("PropCIs", dependencies = TRUE)


# loading dataset
data_rec <- read.table("data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)

# Formula returning the variance of the product of two random variables
var_prod <- function(m1,m2,s1,s2,r12) (m1*s2)^2 + (m2*s1)^2 + 2*m1*m2*s1*s2*r12 + (1*r12^2)*(s1*s2)^2

species_list <- c(
  "Ae_bicornis",
  "Ae_caudata",
  "Ae_comosa",
  "Ae_longissima",
  "Ae_mutica",
  "Ae_searsii",
  "Ae_sharonensis",
  "Ae_speltoides",
  "Ae_tauschii",
  "Ae_umbellulata",
  "Ae_uniaristata",
  "T_boeticum",
  "T_urartu"
)

# Analysis for piS, Ncat number of categories

result <- c()

Ncat <- 100
Nboot <- 1000

for(SP in c(1:13)) {
  
  # Choice of the species
  # SP <- 11 # An example
  SPECIES <- species_list[SP]
  data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
  data_focal$V1 <- gsub("_simExt","",data_focal$V1)
  data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
  data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
  data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
  names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
  # Filtering out too divergent orthologs. After inspection by eye because 0.35 is an appropriate threshold
  THRESHOLD <- 0.35
  data_focal <- data_focal[data_focal$t<THRESHOLD,]
  FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                 "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
  PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
  data_pol <- read.table(PATH,header=T)
  # To suppress the species name in some variable names 
  names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
  data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
  data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
  data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
  # Preparing the datafile
  mydata <- merge(data_rec,data_hordeum,by="Gene")
  mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
  mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
  mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
  names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
  names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
  mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
  mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
  mydata$PIS <- mydata$piS*mydata$nb_complete_site
  mydata$PIN <- mydata$piN*mydata$nb_complete_site
  mydata$Ldiv <- mydata$S
  mydata$SIM <- mydata$dS*mydata$Ldiv
  mydataWindow <- aggregate(mydata,by=list(mydata$CatRec),FUN = function(x) mean(x,na.rm=T))
  mydataWindow$piSyn <- mydataWindow$PIS/mydataWindow$nb_complete_site
  mydataWindow$piNonSyn <- mydataWindow$PIN/mydataWindow$nb_complete_site
  mydataWindow$f0 <- mydataWindow$piNonSyn/mydataWindow$piSyn
  mydataWindow$div <- mydataWindow$S*mydataWindow$dS/mydataWindow$S
  mydataWindow$Ne <- 12000000*mydataWindow$piSyn/mydataWindow$div
  mydataWindow$reldiv <- mydataWindow$div/mean(mydataWindow$div)
  mydataWindow$piSyn_cor <- mydataWindow$piSyn/mydataWindow$reldiv
  # Filtering
  filter <- which(!is.na(mydataWindow$piSyn_cor) &!is.infinite(mydataWindow$piSyn_cor) & mydataWindow$piSyn_cor>0)
  x <- mydataWindow$Recombination[filter]
  y <- mydataWindow$piSyn_cor[filter]
  n <- mydataWindow$nb_complete_site[filter]
# Bootstrapping the results
  coeffboot <- c()
  i <- 0
  repeat{
    boot <- sample(length(x),Nboot,replace = T)
    xboot <- x[boot]
    yboot <- y[boot]
    nboot <- n[boot]
    piobsboot <- sum(yboot*nboot)/sum(nboot)
    ainit <- mean(yboot)*3
    binit <- sample(x = c(0.1,0.5,1,2,3,4,5),1)
    blim <- 5
    cinit <- sample(x = c(0.001,0.005,0.01,0.05,0.1,0.5),1)
    w <- nboot/(yboot*(1-yboot))
    w <- w/mean(w)
    ZERO <- 0.00001
    tryCatch({
      model <- 
        nls(yboot ~ a/(1+exp(b-(4*c*xboot/a))),weights = w,start = list(a=ainit,b=binit,c=cinit),algorithm = "port",lower=c(ZERO,-blim,ZERO),upper=c(1,blim,100),
            nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024, printEval = TRUE, warnOnly = FALSE) )
      model_sum <- summary(model)
      if(model_sum$parameters[2]<= 0.99*blim) {
        coeffboot <- rbind(coeffboot,c(model_sum$parameters[,1],model_sum$parameters[1,1]/piobsboot))
        i <- i+1
        #print(i)
      }
    }, error=function(e){})
    if(i>=Nboot) break
  }
  # Collecting the results
  meanboot <- apply(coeffboot,2,mean)
  sdboot <- apply(coeffboot,2,sd)
  ci_inf <- apply(coeffboot,2, function(x) sort(x)[round(length(x)*0.025)])
  ci_sup <- apply(coeffboot,2, function(x) sort(x)[round(length(x)*0.975)])
  aestim <- meanboot[1]
  bestim <- meanboot[2]
  cestim <- meanboot[3]
  y_ci <- mapply(function(x,n) exactci(x,n,0.95)$conf.int,x=n*y,n=n)
  # Mean fit
  fit.pi <- function(x) aestim/(1+exp(bestim-(4*cestim*x/aestim)))
  # Plot for each species
  pdf(paste("figures/sup_mat/",SPECIES,"_Ncat=",Ncat,"_piS_withCI.pdf",sep=""),width = 8,height = 6)
    plotCI(x,y,li=y_ci[1,],ui=y_ci[2,],sfrac = 0.002,main=SPECIES,xlab="Recombination (cM/Mb)",ylab=expression(pi[S]),cex.lab=1.2,col="grey")
    points(x,y,pch=16,col="grey",cex=0.8)
    points(x,y)
    points(c(-100:2000)/1000,fit.pi(c(-100:2000)/1000),type="l",lwd=2,col="red")
  dev.off()
  pdf(paste("figures/sup_mat/",SPECIES,"_Ncat=",Ncat,"_piS.pdf",sep=""),width = 8,height = 6)
    plot(x,y,main=SPECIES,xlab="Recombination (cM/Mb)",ylab=expression(pi[S]),cex.lab=1.2)
    points(x,y,pch=16,col="grey",cex=0.8)
    points(c(-100:2000)/1000,fit.pi(c(-100:2000)/1000),type="l",lwd=2,col="red")
  dev.off()
  # Exporting the results
  stat <- data.frame(t(c(SPECIES,meanboot,sdboot,ci_inf,ci_sup)))
  names(stat) <- c("species","pimax","b","slope","ratio","pimax_sd","b_sd","slope_sd","ratio_sd","pimax_inf","b_inf","slope_inf","ratio_inf","pimax_sup","b_sup","slope_sup","ratio_sup")
  result <-rbind(result,stat)
  write.table(result,paste("outputs/recombination/hordeum/stat_recombination_per_species_Ncat=",Ncat,".txt",sep=""),sep = "\t",quote = F,row.names = F)
  
  print(SPECIES)
}