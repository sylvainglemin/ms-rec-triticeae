# Script to SFS from SNP table
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated April 2024


# Function to project a SFS from n to m categories
reduceSFS <- function(sfs,m){
  n <- length(sfs)
  if(n==m) return(sfs)
  output <- rep(0,m)
  index<-c(1:m)
  for(k in 1:n){
    output <- output + sfs[k]*choose(k,index)*choose(n-k,m-index)/choose(n,m)
  }
  return(output)
}



# To choose the files
SPECIES <- "T_urartu"
NUMBER <- 10
donnees <- read.table(paste("data/polymorphism/",SPECIES,"/",SPECIES,"_T_caputmedusae_snp_table_n",NUMBER,".txt",sep=""),header=T)
donneesATCG <- read.table(paste("data/polymorphism/",SPECIES,"/",SPECIES,"_snp_table_ATGC.txt",sep=""),header=T)
donnees <- merge(donnees,donneesATCG,by=c("gene","position"))
donnees$gene <- gsub(".fasta.fst.clean.fst","",donnees$gene)

data_rec <- read.table("data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)
data_focal <- read.table(paste("outputs/orthology/",SPECIES,"_Hordeum_reciprocalbestblast.txt",sep=""),header=T)
data_focal$Query <- gsub("_simExt","",data_focal$Query)
data_focal$Query <- gsub("_lgOrf","",data_focal$Query)
data_focal$Query <- gsub("_ESTScan","",data_focal$Query)

# Mapped corresponds to genes that can be mapped on Hordeum
mapped <- merge(data_rec,data_hordeum,by="Gene")
mapped <- merge(mapped,data_focal,by.x="Gene",by.y="Subject")
mapped <- merge(mapped,donnees,by.x="Query",by.y="gene")
names(mapped)[c(1,2)] <- c("gene","Hordeum_gene")


# To get the count of syn and non-syn position
contig <- read.table(paste("data/polymorphism/",SPECIES,"/dNdSpiNpiS_output_outgroup",sep=""),header=T)
infocontig <- contig[,c("Contig_name","nS","nN")]


############################ #
# SFS on the whole dataset ###
############################ #

mydata <- merge(donnees,infocontig,by.x="gene",by.y="Contig_name")


# Sample size to project the SFS
# To be chosen as a function of the dataset
hist(mydata$n_cult)
Nmin <- 19
Nchrom <- max(mydata$n_cult)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)


# SFS ####

l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)

sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}

FILENAME <- paste("outputs/dfe/sfs/",SPECIES,"_total_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# Bootstrap files ####

# Number of bootstraps to be choosen
Nboot <- 1000
Ngene <- length(levels(mydata$gene))

for(i in 1:Nboot) {
  sampled_genes <- data.frame(gene=sample(levels(mydata$gene),Ngene,replace = T))
  databoot <- merge(mydata,sampled_genes,by="gene")
  l_S <- sum(aggregate(databoot$nS,by=list(databoot$gene),FUN=mean)$x)
  l_NS <- sum(aggregate(databoot$nN,by=list(databoot$gene),FUN=mean)$x)
  sfs_S <- rep(0,Nmin-1)
  sfs_NS <- rep(0,Nmin-1)
  for(n in Nmin:Nchrom){
    data_S <- databoot[databoot$type=="S", ]
    data_S <- data_S[data_S$n_cult==n,]
    data_NS <- databoot[databoot$type=="NS", ]
    data_NS <- data_NS[data_NS$n_cult==n,]
    if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
    sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
  }
  FILENAME <- paste("outputs/dfe/sfs/",SPECIES,"_total_polyDFE_boot",i,".txt",sep="")
  write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
  write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
  write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)
}



############################ #
# SFS on mapped contigs   ####
############################ #


mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
hist(mydata$n_cult)

Nmin <- 19 # To be chosen as a function of the dataset
Nchrom <- max(mydata$n_cult)

mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)


# SFS ####


l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)

sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n & !is.na(data_S$n_cult),]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n & !is.na(data_S$n_cult),]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}

FILENAME <- paste("outputs/dfe/sfs/",SPECIES,"_mapped_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# Bootstrap files ####

Nboot <- 1000
Ngene <- length(levels(mydata$gene))

for(i in 1:Nboot) {
  sampled_genes <- data.frame(gene=sample(levels(mydata$gene),Ngene,replace = T))
  databoot <- merge(mydata,sampled_genes,by="gene")
  l_S <- sum(aggregate(databoot$nS,by=list(databoot$gene),FUN=mean)$x)
  l_NS <- sum(aggregate(databoot$nN,by=list(databoot$gene),FUN=mean)$x)
  sfs_S <- rep(0,Nmin-1)
  sfs_NS <- rep(0,Nmin-1)
  for(n in Nmin:Nchrom){
    data_S <- databoot[databoot$type=="S", ]
    data_S <- data_S[data_S$n_cult==n & !is.na(data_S$n_cult),]
    data_NS <- databoot[databoot$type=="NS", ]
    data_NS <- data_NS[data_NS$n_cult==n & !is.na(data_S$n_cult),]
    if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
    sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
  }
  FILENAME <- paste("outputs/dfe/sfs",SPECIES,"_mapped_polyDFE_boot",i,".txt",sep="")
  write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
  write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
  write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)
}


##################################### #
# SFS on mapped contigs high rec   ####
##################################### #

# Thereshold for the two categories of recombination
THRESHOLDREC <- 0.5

mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
hist(mydata$Loess2,breaks=50)
mydata <- mydata[mydata$Loess2>THRESHOLDREC,]

hist(mydata$n_cult)

Nmin <- 19 # To be chosen as a function of the dataset

mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)

Nchrom <- max(mydata$n_cult,na.rm = T)

# SFS ####

l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)

sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}

FILENAME <- paste("outputs/dfe/sfs/",SPECIES,"_highrec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# Bootstrap files ####

Nboot <- 1000
Ngene <- length(levels(mydata$gene))

for(i in 1:Nboot) {
  sampled_genes <- data.frame(gene=sample(levels(mydata$gene),Ngene,replace = T))
  databoot <- merge(mydata,sampled_genes,by="gene")
  l_S <- sum(aggregate(databoot$nS,by=list(databoot$gene),FUN=mean)$x)
  l_NS <- sum(aggregate(databoot$nN,by=list(databoot$gene),FUN=mean)$x)
  sfs_S <- rep(0,Nmin-1)
  sfs_NS <- rep(0,Nmin-1)
  for(n in Nmin:Nchrom){
    data_S <- databoot[databoot$type=="S", ]
    data_S <- data_S[data_S$n_cult==n & !is.na(data_S$n_cult),]
    data_NS <- databoot[databoot$type=="NS", ]
    data_NS <- data_NS[data_NS$n_cult==n & !is.na(data_S$n_cult),]
    if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
    sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
  }
  FILENAME <- paste("outputs/dfe/sfs/",SPECIES,"_highrec_polyDFE_boot",i,".txt",sep="")
  write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
  write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
  write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)
}



##################################### #
# SFS on mapped contigs low rec   ####
##################################### #

# Thereshold for the two categories of recombination
THRESHOLDREC <- 0.5

mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2<=THRESHOLDREC,]

hist(mydata$n_cult)

Nmin <- 19 # To be chosen as a function of the dataset

mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)

Nchrom <- max(mydata$n_cult,na.rm = T)

# SFS ####

l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)

sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}

FILENAME <- paste(SPECIES,"_lowrec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# Bootstrap files ####

Nboot <- 1000
Ngene <- length(levels(mydata$gene))

for(i in 1:Nboot) {
  sampled_genes <- data.frame(gene=sample(levels(mydata$gene),Ngene,replace = T))
  databoot <- merge(mydata,sampled_genes,by="gene")
  l_S <- sum(aggregate(databoot$nS,by=list(databoot$gene),FUN=mean)$x)
  l_NS <- sum(aggregate(databoot$nN,by=list(databoot$gene),FUN=mean)$x)
  sfs_S <- rep(0,Nmin-1)
  sfs_NS <- rep(0,Nmin-1)
  for(n in Nmin:Nchrom){
    data_S <- databoot[databoot$type=="S", ]
    data_S <- data_S[data_S$n_cult==n & !is.na(data_S$n_cult),]
    data_NS <- databoot[databoot$type=="NS", ]
    data_NS <- data_NS[data_NS$n_cult==n & !is.na(data_S$n_cult),]
    if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
    sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
    sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
  }
  FILENAME <- paste(SPECIES,"_lowrec_polyDFE_boot",i,".txt",sep="")
  write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
  write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
  write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)
}





########################################### #
# SFS as a function of base composition  ####
########################################### #
Nmin <- 19 # To be chosen as a function of the dataset

# GC conservative
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/T"|mydata$typeATGC=="T/A"|mydata$typeATGC=="C/G"|mydata$typeATGC=="G/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WWSS_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)

# W->S
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/G"|mydata$typeATGC=="A/C"|mydata$typeATGC=="T/G"|mydata$typeATGC=="T/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WS_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)

# S->W
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="G/A"|mydata$typeATGC=="C/A"|mydata$typeATGC=="G/T"|mydata$typeATGC=="C/T"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_SW_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)



###################################################################   #
# SFS on mapped contigs high rec as a function of base compostion  ####
###################################################################   #
Nmin <- 19 # To be chosen as a function of the dataset

# Thereshold for the two categories of recombination
THRESHOLDREC <- 0.5

# GC conservative
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2>THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/T"|mydata$typeATGC=="T/A"|mydata$typeATGC=="C/G"|mydata$typeATGC=="G/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WWSS_HighRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# W->S
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2>THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/G"|mydata$typeATGC=="A/C"|mydata$typeATGC=="T/G"|mydata$typeATGC=="T/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WS_HighRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# S->W
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2>THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="G/A"|mydata$typeATGC=="C/A"|mydata$typeATGC=="G/T"|mydata$typeATGC=="C/T"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_SW_HighRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)



###################################################################   #
# SFS on mapped contigs low rec as a function of base compostion   ####
###################################################################   #
Nmin <- 19 # To be chosen as a function of the dataset

# Thereshold for the two categories of recombination
THRESHOLDREC <- 0.5

# GC conservative
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2<=THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/T"|mydata$typeATGC=="T/A"|mydata$typeATGC=="C/G"|mydata$typeATGC=="G/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WWSS_LowHighRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# W->S
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2<=THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="A/G"|mydata$typeATGC=="A/C"|mydata$typeATGC=="T/G"|mydata$typeATGC=="T/C"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_WS_LowRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)


# S->W
mydata <- merge(mapped,infocontig,by.x="gene",by.y="Contig_name")
mydata <- mydata[mydata$Loess2<=THRESHOLDREC,]
mydata <- mydata[mydata$polarized==0 & (mydata$typeATGC=="G/A"|mydata$typeATGC=="C/A"|mydata$typeATGC=="G/T"|mydata$typeATGC=="C/T"),]
Nchrom <- max(mydata$n_cult,na.rm = T)
mydata <- mydata[mydata$folded==0 & mydata$n_cult>=Nmin,]
mydata$gene <- as.factor(mydata$gene)
l_S <- sum(aggregate(mydata$nS,by=list(mydata$gene),FUN=mean)$x)
l_NS <- sum(aggregate(mydata$nN,by=list(mydata$gene),FUN=mean)$x)
sfs_S <- rep(0,Nmin-1)
sfs_NS <- rep(0,Nmin-1)
for(n in Nmin:Nchrom){
  data_S <- mydata[mydata$type=="S", ]
  data_S <- data_S[data_S$n_cult==n,]
  data_NS <- mydata[mydata$type=="NS", ]
  data_NS <- data_NS[data_NS$n_cult==n,]
  if(dim(data_S)[1]==0) sfs_Sn <- rep(0,Nmin-1) else sfs_Sn <- hist(data_S$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  if(dim(data_NS)[1]==0) sfs_NSn <- rep(0,Nmin-1) else sfs_NSn <- hist(data_NS$x_cult,breaks = (c(1:22)-0.5),plot = F)$counts[1:(n-1)]
  sfs_S <- sfs_S + reduceSFS(sfs_Sn,Nmin-1)
  sfs_NS <- sfs_NS + reduceSFS(sfs_NSn,Nmin-1)
}
FILENAME <- paste(SPECIES,"_SW_LowRec_polyDFE.txt",sep="")
write(x=c(1,1,Nmin),file=FILENAME,sep="\t",ncol=3)
write(x=c(sfs_S,l_S),file=FILENAME,sep="\t",ncol=(length(sfs_S)+1),append = T)
write(x=c(sfs_NS,l_NS),file=FILENAME,sep="\t",ncol=(length(sfs_NS)+1),append = T)



