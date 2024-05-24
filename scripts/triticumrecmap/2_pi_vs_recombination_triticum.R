# Analysis of the relationship between recombination rate and polymorphism patterns
# Plots with recombination maps of the three T aestivum sub-genomes
# Sylvain Glémin (CNRS Rennes, France)
# April 2024

if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)
if (!require("gridExtra"))   install.packages("gridExtra", dependencies = TRUE)


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



# Analysis for piS Ncat number of categories ####

SizeMb <- 1*10^6
Ncat <- 50

for(SPECIES in species_list) {
  for(GENOME in c("A","B","D")){
    #SPECIES <- species_list[2] # example
    # GENOME  "D"
    # Loading and preparing datasets ####
    data_rec <- read.table(file ="data/recombination/RecombinationRates_AllTriticumGenes.txt",header=T)
    data_rec <- data_rec[substr(data_rec$map,2,2)==GENOME,]
    names(data_rec)[1] <- "Gene"
    mareymap <- read.table("data/recombination/MareyMap_Triticum_aestivum_GutierrezGonzalez2019.txt",header=T)
    # Removing invalid markers from the Marey map
    #mareymap <- mareymap[mareymap$vld==T,]
    data_triticum <- read.table("data/recombination/TriticumGenes.txt",header=T)
    data_triticum <- data_triticum[substr(data_triticum$Chromosome,2,2)==GENOME,]

    data_focal <- read.table(paste("outputs/orthology/triticum/",SPECIES,"_Taestivum",GENOME,"_reciprocalbestblast.txt",sep=""),header=T,sep ="",fill=T)
    data_focal$Query <- gsub("_simExt","",data_focal$Query)
    data_focal$Query <- gsub("_lgORF","",data_focal$Query)
    data_focal$Query <- gsub("_ESTScan","",data_focal$Query)
    data_focal <- data_focal[data_focal$Evalue<10^(-20),]
    data_focal <- data_focal[,c(1,2)]
    names(data_focal) <- c("Gene_focal","Gene")
    
    # Loading polymorphism data
    FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                   "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
    data_pol <- read.table(paste("data/polymorphism/",SPECIES,"/",FILE,sep=""),header=T)
    # To suppress the species name in some variable names 
    names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
    data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
    data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
    data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
    # Final dataset
    mydata <- merge(data_rec,data_triticum,by="Gene")
    mydata <- merge(mydata,data_focal,by="Gene")
    mydata <- merge(mydata,data_pol,by.x="Gene_focal",by.y="Contig_name")
    # Cleaning the data.frame
    mydata <- mydata[-which(names(mydata)==c("map") | names(mydata)==c("start") | names(mydata)==c("end") | names(mydata)==c("set"))]
    mydata$recRate <- ifelse(mydata$recRate<0,0,mydata$recRate)
    mydata$CatRec  <- cut(mydata$recRate,breaks = unique(quantile(mydata$recRate, c(0:Ncat)/Ncat,na.rm=T)))
    length(levels(mydata$CatRec))
    chromlength <- aggregate(mydata$Start,by=list("Chromosome"=mydata$Chromosome),FUN = max)
    names(chromlength)[2] <- "ChromLength"
    mydata <- merge(mydata,chromlength,by="Chromosome")
    mydata$posMb  <- SizeMb*round(mydata$Start/SizeMb)
    mydata$posMbrel  <- mydata$posMb/mydata$ChromLength
    mydata$PIS <- mydata$piS*mydata$nb_complete_site
    mydata$PIN <- mydata$piN*mydata$nb_complete_site
    mydataMb <- aggregate(mydata,by=list("posMb"=mydata$posMb,"Chromosome"=mydata$Chromosome),FUN = function(x) mean(x,na.rm=T))
    mydataMb <- mydataMb[-c(3:5)]
    mydataMb$NbGene <- unlist(aggregate(mydata,by=list("posMb"=mydata$posMb,"Chromosome"=mydata$Chromosome),FUN = length)[10])
    mydataMb$piSyn <- mydataMb$PIS/mydataMb$nb_complete_site
    mydataMb$piNonSyn <- mydataMb$PIN/mydataMb$nb_complete_site
    mydataMb$f0 <- mydataMb$piNonSyn/mydataMb$piSyn
    
    # Plot
    G <- ggplot(data=mydataMb,aes(x=posMb)) + theme_classic() +
      facet_wrap(~Chromosome,scales = "free_x",) + 
      geom_point(data=mydataMb,aes(x=posMb,y=piSyn),col="grey") +
      geom_point(data=mydataMb,aes(x=posMb,y=piSyn),pch=1,col="grey50") +
      geom_smooth(data=mydataMb,aes(y=piSyn),method = "loess", span = 0.1,se = F,col="blue",cex=1.2) +
      geom_line(data=mydataMb,aes(x=Start,y=recRate*3*10^4),col="black",cex=1,lty=1) +
      scale_y_continuous(limits = c(0,0.1),sec.axis = sec_axis(~.*33, name="Recombination rate (in cM/Mb)")) + 
      #theme(axis.text=element_text(size = 14),axis.line = element_line(linetype=1)) +
      theme(axis.text=element_text(size = 10)) +
      theme(axis.title = element_text(size = 12),plot.title = element_text(size = 16,face="italic",hjust=0.5)) + 
      xlab(NULL) + ylab(expression(pi[S]))
    ggsave(filename = paste("figures/additional/triticum_recmap/",SPECIES,"_piS_alongchrom_genome",GENOME,".pdf",sep=""),plot = G,width = 12,height = 10)
  }
}
  
  
  