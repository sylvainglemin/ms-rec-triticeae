# Patterns of pi along chromosomes
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022


# packages to install
if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)
if (!require("ggpubr"))   install.packages("ggpubr", dependencies = TRUE)
if (!require("gridExtra"))   install.packages("gridExtra", dependencies = TRUE)


#Loading datasets
data_rec <- read.table("data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
data_rec$Recombination <- ifelse(data_rec$Recombination<0,0,data_rec$Recombination)
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)

# Graph of recombination along chromosomes ####

listsp1 <- c("Ae_bicornis","Ae_caudata","Ae_comosa","Ae_longissima",
             "Ae_mutica","Ae_searsii","Ae_sharonensis","Ae_speltoides",
             "Ae_tauschii","Ae_umbellulata","Ae_uniaristata",
             "T_boeticum","T_urartu")
listsp2 <- c("Aegilops bicornis","Aegilops caudata","Aegilops comosa","Aegilops longissima",
             "Aegilops mutica","Aegilops searsii","Aegilops sharonensis","Aegilops speltoides",
             "Aegilops tauschii","Aegilops umbellutata","Aegilops uniaristata",
             "Triticum monoccocum","Triticum urartu")

# Choosing options ####
SizeMb <- 1
NcatRec <- 50

for(i in c(1:13)) {
  
  SPECIES <- listsp1[i]
  SP_title <- listsp2[i]
  
  # Grooping by windows along chromosomes ####
  
  # Prepare datasets for plot along chromosome ###
  data_focal <- read.table(paste("outputs/orthology/",SPECIES,"_Hordeum_reciprocalbestblast.txt",sep=""),header=T)
  data_focal$Query <- gsub("_simExt","",data_focal$Query)
  data_focal$Query <- gsub("_lgOrf","",data_focal$Query)
  data_focal$Query <- gsub("_ESTScan","",data_focal$Query)
  data_focal <- data_focal[,c("Query","Subject","Coverage.x","Similarity.x")]
  FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                 "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
  PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
  data_pol <- read.table(PATH,header=T)
  # To suppress the species name in some variable names 
  names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
  data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
  data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
  data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
  # Prepare file
  mydata <- merge(data_rec,data_hordeum,by="Gene")
  mydata <- merge(mydata,data_focal,by.x="Gene",by.y="Subject")
  mydata <- merge(mydata,data_pol,by.x="Query",by.y="Contig_name")
  mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
  names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
  names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
  chromlength <- aggregate(mydata$Start,by=list("Chromosome"=mydata$Chromosome),FUN = max)
  names(chromlength)[2] <- "ChromLength"
  mydata <- merge(mydata,chromlength,by="Chromosome")
  mydata$posMb  <- SizeMb*round(mydata$Start/SizeMb)
  mydata$posMbrel  <- mydata$posMb/mydata$ChromLength
  mydata$PIS <- mydata$piS*mydata$nb_complete_site
  mydata$PIN <- mydata$piN*mydata$nb_complete_site
  mydata$Ldiv <- mydata$Coverage.x*mydata$nb_complete_site
  mydata$SIM <- mydata$Similarity.x*mydata$Ldiv
  
  mydataMb <- aggregate(mydata,by=list("posMb"=mydata$posMb,"Chromosome"=mydata$Chromosome),FUN = function(x) mean(x,na.rm=T))
  mydataMb <- mydataMb[-c(3:5)]
  mydataMb$NbGene <- unlist(aggregate(mydata,by=list("posMb"=mydata$posMb,"Chromosome"=mydata$Chromosome),FUN = length)[10])
  mydataMb$piSyn <- mydataMb$PIS/mydataMb$nb_complete_site
  mydataMb$piNonSyn <- mydataMb$PIN/mydataMb$nb_complete_site
  mydataMb$f0 <- mydataMb$piNonSyn/mydataMb$piSyn
  mydataMb$div <- (1 - 0.01*mydataMb$SIM/mydataMb$Ldiv)/60000000
  mydataMb$Ne <- mydataMb$piSyn/mydataMb$div
  
  
  # Plot of recombination and piS along chromosome ###
  G1 <- ggplot(data=mydataMb[mydataMb$Chromosome!="7H",],aes(x=posMb)) + theme_classic() +
    facet_wrap(~Chromosome,scales = "free_x",) + 
    geom_point(data=mydataMb[mydataMb$Chromosome!="7H",],aes(x=posMb,y=piSyn),col="grey") +
    geom_point(data=mydataMb[mydataMb$Chromosome!="7H",],aes(x=posMb,y=piSyn),pch=1,col="grey50") +
    geom_smooth(data=mydataMb[mydataMb$Chromosome!="7H",],aes(y=piSyn),method = "loess", span = 0.1,se = F,col="blue",cex=1.2) +
    geom_line(data=data_rec[data_rec$Chromosome!="7H",],aes(x=Start,y=Recombination/10),col="black",cex=1,lty=2) +
    scale_y_continuous(limits = c(0,0.12),sec.axis = sec_axis(~.*10, name="Recombination rate (in cM/Mb)")) + 
    #theme(axis.text=element_text(size = 14),axis.line = element_line(linetype=1)) +
    theme(axis.text=element_text(size = 10)) +
    theme(axis.title = element_text(size = 12),plot.title = element_text(size = 16,face="italic",hjust=0.5)) + 
    xlab(NULL) + ylab(expression(pi[S])) + ggtitle(SP_title)
  
  G1bis <- ggplot(data=mydataMb[mydataMb$Chromosome=="7H",],aes(x=posMb)) + theme_classic() +
    facet_wrap(~Chromosome,scales = "free_x") + 
    geom_point(data=mydataMb[mydataMb$Chromosome=="7H",],aes(x=posMb,y=piSyn),col="grey") +
    geom_point(data=mydataMb[mydataMb$Chromosome=="7H",],aes(x=posMb,y=piSyn),pch=1,col="grey50") +
    geom_smooth(data=mydataMb[mydataMb$Chromosome=="7H",],aes(y=piSyn),method = "loess", span = 0.1,se = F,col="blue",cex=1.2) +
    geom_line(data=data_rec[data_rec$Chromosome=="7H",],aes(x=Start,y=Recombination/10),col="black",cex=1,lty=2) +
    scale_y_continuous(limits = c(0,0.12),sec.axis = sec_axis(~.*10, name=NULL)) + 
    #theme(axis.text=element_text(size = 14),axis.line = element_line(linetype=1)) +
    theme(axis.text=element_text(size = 10)) +
    theme(axis.title = element_text(size = 12),plot.title = element_text(size = 16,face="italic",hjust=0.5)) + 
    xlab("Chromosome position (in bp)") + ylab(" ") + ggtitle(NULL)
  
  
  # Grouping as a function of recombination rate ####
  
  
  # Prepare dataset ###
  mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:NcatRec)/NcatRec,na.rm=T)))
  mydata$PIS <- mydata$piS*mydata$nb_complete_site
  mydata$PIN <- mydata$piN*mydata$nb_complete_site
  mydata$Ldiv <- mydata$Coverage.x*mydata$nb_complete_site
  mydata$SIM <- mydata$mydata$Ldiv
  mydataWindow <- aggregate(mydata,by=list(mydata$CatRec),FUN = function(x) mean(x,na.rm=T))
  mydataWindow$piSyn <- mydataWindow$PIS/mydataWindow$nb_complete_site
  mydataWindow$piNonSyn <- mydataWindow$PIN/mydataWindow$nb_complete_site
  mydataWindow$f0 <- mydataWindow$piNonSyn/mydataWindow$piSyn
  # To avoid aberrant windows
  filter <- which(!is.na(mydataWindow$piSyn) &!is.infinite(mydataWindow$f0))
  
  # Fitting a logistic model on the data
  x <- mydataWindow$Recombination[filter]
  y <- mydataWindow$piSyn[filter]
  model <- nls(y ~ a/(1+b*exp(-(c*x))),start = list(a=0.002,b=0.1,c=1),algorithm = "port",lower=c(0,0,0),upper=c(100,100,100))
  coeff <- summary(model)$parameters[,1]
  sd <- summary(model)$parameters[,2]
  pimax <- coeff[1]
  pimax.sd <- sd[1]
  pimin <- coeff[1]/(1+coeff[2])
  rec.pi <- coeff[1]*coeff[2]*coeff[3]/(1+coeff[2])^2
  slope.pi <- coeff[1]*coeff[3]/4
  medianpi <- median(y)
  fit.pi <- function(x) coeff[1]/(1+coeff[2]*exp(-(coeff[3]*x)))
  
  
  G2 <- ggplot(mydataWindow,aes(x=Recombination,y=piSyn)) +theme_classic() +
    geom_point(col="grey",cex=2) + geom_point(pch=1,col="grey50",cex=2) +
    geom_line(aes(x=Recombination,y=fit.pi(Recombination)),cex=1.2,col="blue") +
    theme(axis.text=element_text(size = 12),axis.line = element_line(linetype=1)) +
    theme(axis.title = element_text(size = 12),plot.title = element_text(size = 16,face="italic",hjust=0.5)) + 
    theme(plot.margin = margin(l = 100,r = 100, t = 20)) +
    xlab("Recombination rate (in cM/Mb)") + ylab(expression(pi[S])) + ggtitle(NULL)
  
  pdf(paste("figures/sup_mat/",SPECIES,"_piS_alongchrom.pdf",sep = ""),width = 12,height = 10)
  grid.arrange(G1,G1bis,G2,heights=c(2,1),widths=c(2,3.2),
               layout_matrix = rbind(c(1, 1),c(2,3))
  )
  #annotate_figure(fig,fig.lab=SP_title,fig.lab.pos = "top",fig.lab.size = 20,fig.lab.face = "italic")
  dev.off()
}

