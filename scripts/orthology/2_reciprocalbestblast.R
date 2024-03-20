# Modified from Pascal Milesi
# All steps are not useful but have been kept from Pascal initial script

#WARNINGS: the output of Blast as to be formated and sorted as follow 
#Blast: -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen gaps'
#Blast ouput has to be sorted with the following command : sort -k1,1 -k2,2 -k7,7n
#meaning: sorted by query name then by subject name then by the position of Query start (increasing order).

# Use the correct PATH here
setwd("PATH")

# Intinial thresholds
INIT_SIM_THRESH <- 70
INIT_EVALUE_THRESH <- 1e-10



############################ #
# Sorting focal 2 Hordeum ####
############################ #


# Best blast thresholds
COV_TRESH <- 0.8
SIM_2nd_THRESH <- 90 # in %
SIM_BEST_THRESH <- 90

LISTFILES <- list.files(pattern = ".fasta_focal2hordeum.blast.txt")

for(FILENAME in LISTFILES[3:19]) {

  OUTPUT <- paste(strsplit(FILENAME,".blast.txt"),".bestblast.txt",sep="")
  dataf <- read.table(FILENAME, sep = "\t", header =F)
  colnames(dataf) <- c("Query","Subject","Similarity","Length","Mismatch","Gap","Qstart","Qend","Sstart","Send","E-value","Bitscore","Querylength","Totgaplentgh")
  dataf2 <- dataf[which(dataf[,3] >= INIT_SIM_THRESH & dataf[,11] <= INIT_EVALUE_THRESH),] # You can chosse the treshold to be more or less convervative the first argument refers to the similarity, the second to the e-value 
  Lqh <- dataf2$Length*dataf2$Querylength
  dataf3 <- cbind.data.frame(dataf2,Lqh)
  dataf4 <- dataf3[order(dataf3$Lqh),]
  # Bit score correction
  mod1 <- lm(dataf4$Bitscore~dataf4$Lqh)
  a <- mod1$coefficients[2]
  b <- mod1$coefficients[1]
  Bitcor <- dataf4$Bitscore/(b+dataf4$Lqh*a)
  dataf5 <- cbind.data.frame(dataf4,Bitcor)
  dataf6 <- dataf5[order(dataf5[,1],dataf5[,2],dataf5[,7]),]
  cov <- dataf6$Length/dataf6$Querylength #computed the coverage by dividing the hit length by query length
  epond <- dataf6[,11]*dataf6[,4] #computed pondered evalue by hit length
  bitpond <- dataf6[,16]*dataf6[,4] #computed ponderd bitscore by hit length
  dataf6 <- cbind.data.frame(dataf6,cov,epond,bitpond)
  # The pattern corresponding to transcript that must be suppress can differ form one ref species to the other -> to be checked
  # To supress transcript redundancy
  pattern <- "[:.:][0-9]*"
  dataf6$Subject <- sub(pattern,"",dataf6$Subject)
  query <- unique(dataf6[,1])

  # Vector initialization
  epondj <- c()
  bpondj <- c()
  sim <- c()
  names <- c()
  #mismatch <- c()
  cov <- c()
  list <- c()
  Qlen <- c()
  m <- 0
  n <- length(dataf6[,1])
  sortiefull <- c()
  #the following loop aims at gathering all hits of a same Subject and compute summary statistique for each pair Query - Subject
  #The percentage of advanecement is printed.
  #unique(dataf2[which(dataf2[,1]==i),2])
  t <- 1 #Counter for printing
  
  for(i in query) {
    
    subsample <- dataf6[which(dataf6$Query==i),]
    subsample <- subsample[order(subsample$Bitcor,decreasing = T),]
    best <- unique(subsample$Subject)
    
    for (j in best) {
      temp <- dataf6[which(dataf6[,1]==i & dataf6[,2]==j),]
      if(length(temp[,1])>1){
        epondj[length(epondj)+1] <- sum(temp[,18])/sum(temp[,4])
        bpondj[length(bpondj)+1] <- sum(temp[,12])/sum(temp[,4])
        totlength <- length(unique(as.vector(unlist(mapply(function(x,y) c(x:y),x=temp[,7],y=temp[,8])))))
        cov[length(cov)+1] <- totlength/temp[1,13]
        mat_sim <- matrix(rep(NA,temp[1,13]*length(temp[,1])),length(temp[,1]),temp[1,13])
        for(jj in 1:length(temp[,1])) {
        mat_sim[jj,c(temp[jj,7]:temp[jj,8])] <- temp[jj,3]
        }
      sim[length(sim)+1] <- mean(apply(mat_sim,2,function(x) mean(x,na.rm = T)),na.rm = T)
      names[length(names)+1] <- as.character(j)
      #mismatch[length(mismatch)+1] <- sum(temp[,5])
      } else {
        epondj[length(epondj)+1] <- temp[,11]
        bpondj[length(bpondj)+1] <- temp[,16]
        cov[length(cov)+1] <- (temp[1,8]-temp[1,7]+1)/temp[1,13]
        sim[length(sim)+1] <- temp[,3]
        names[length(names)+1] <- as.character(j)
        #mismatch[length(mismatch)+1] <- temp[,5]
      }
      list[length(list)+1] <- as.character(i)
      Qlen[length(Qlen)+1] <- temp[1,13]
    }
  
    sortie <- cbind.data.frame(list,names,Qlen,bpondj,epondj,cov,sim)
    # To keep only the best blast corresponding to thresholds criteria
    if(dim(sortie[sortie$cov > COV_TRESH & sortie$sim > SIM_2nd_THRESH,])[1]==1) {
      sortiefiltered <- sortie[sortie$cov > COV_TRESH & sortie$sim > SIM_BEST_THRESH,]
      sortiefull <- rbind(sortiefull,sortiefiltered)
      }
    # Printing progression
    #if(t%/%100 == t/100) print(paste(round(100*t/length(query)),"%",sep=" "))
    t <- t+1
    #re-initialization
    epondj <- c()
    bpondj <- c()
    sim <- c()
    names <- c()
  # mismatch <- c()
    cov <- c()
    list <- c()
    Qlen <- c()
  }
  print(FILENAME)
  names(sortiefull) <- c("Query","Subject","Query_length","BitScore_corrected","Evalue","Coverage","Similarity")
  write.table(x=sortiefull,file=OUTPUT,quote = F,row.names = F,sep="\t")
}




############################ #
# Sorting Hordeum 2 Focal ####
############################ #


# Here the coverage is computed on the
# Best blast thresholds
COV_TRESH <- 0.8
SIM_2nd_THRESH <- 90 # in %
SIM_BEST_THRESH <- 90


LISTFILES <- list.files(pattern = ".fasta_hordeum2focal.blast.txt")

for(FILENAME in LISTFILES[9:19]) {
  
  OUTPUT <- paste(strsplit(FILENAME,".blast.txt"),".bestblast.txt",sep="")
  dataf <- read.table(FILENAME, sep = "\t", header =F)
  colnames(dataf) <- c("Query","Subject","Similarity","Length","Mismatch","Gap","Qstart","Qend","Sstart","Send","E-value","Bitscore","Subjectlength","Totgaplentgh")
  dataf2 <- dataf[which(dataf[,3] >= INIT_SIM_THRESH & dataf[,11] <= INIT_EVALUE_THRESH),] # You can chosse the treshold to be more or less convervative the first argument refers to the similarity, the second to the e-value 
  Lqh <- dataf2$Length*dataf2$Subjectlength
  dataf3 <- cbind.data.frame(dataf2,Lqh)
  dataf4 <- dataf3[order(dataf3$Lqh),]
  # Bit score correction
  mod1 <- lm(dataf4$Bitscore~dataf4$Lqh)
  a <- mod1$coefficients[2]
  b <- mod1$coefficients[1]
  Bitcor <- dataf4$Bitscore/(b+dataf4$Lqh*a)
  dataf5 <- cbind.data.frame(dataf4,Bitcor)
  dataf6 <- dataf5[order(dataf5[,1],dataf5[,2],dataf5[,7]),]
  cov <- dataf6$Length/dataf6$Subjectlength #computed the coverage by dividing the hit length by query length
  epond <- dataf6[,11]*dataf6[,4] #computed pondered evalue by hit length
  bitpond <- dataf6[,16]*dataf6[,4] #computed ponderd bitscore by hit length
  dataf6 <- cbind.data.frame(dataf6,cov,epond,bitpond)
  # The pattern corresponding to transcript that must be suppress can differ form one ref species to the other -> to be checked
  # To supress transcript redundancy
  pattern <- "[:.:][0-9]*"
  dataf6$Query <- sub(pattern,"",dataf6$Query)
  query <- unique(dataf6[,1])
  
  # Vector initialization
  epondj <- c()
  bpondj <- c()
  sim <- c()
  names <- c()
  #mismatch <- c()
  cov <- c()
  list <- c()
  Qlen <- c()
  m <- 0
  n <- length(dataf6[,1])
  sortiefull <- c()
  #the following loop aims at gathering all hits of a same Subject and compute summary statistique for each pair Query - Subject
  #The percentage of advanecement is printed.
  #unique(dataf2[which(dataf2[,1]==i),2])
  t <- 1 #Counter for printing
  
  for(i in query) {
    
    subsample <- dataf6[which(dataf6$Query==i),]
    subsample <- subsample[order(subsample$Bitcor,decreasing = T),]
    best <- unique(subsample$Subject)
    
    for (j in best) {
      temp <- dataf6[which(dataf6[,1]==i & dataf6[,2]==j),]
      if(length(temp[,1])>1){
        epondj[length(epondj)+1] <- sum(temp[,18])/sum(temp[,4])
        bpondj[length(bpondj)+1] <- sum(temp[,12])/sum(temp[,4])
        totlength <- length(unique(as.vector(unlist(mapply(function(x,y) c(x:y),x=temp[,9],y=temp[,10])))))
        cov[length(cov)+1] <- totlength/temp[1,13]
        mat_sim <- matrix(rep(NA,temp[1,13]*length(temp[,1])),length(temp[,1]),temp[1,13])
        for(jj in 1:length(temp[,1])) {
          mat_sim[jj,c(temp[jj,9]:temp[jj,10])] <- temp[jj,3]
        }
        sim[length(sim)+1] <- mean(apply(mat_sim,2,function(x) mean(x,na.rm = T)),na.rm = T)
        names[length(names)+1] <- as.character(j)
        #mismatch[length(mismatch)+1] <- sum(temp[,5])
      } else {
        epondj[length(epondj)+1] <- temp[,11]
        bpondj[length(bpondj)+1] <- temp[,16]
        cov[length(cov)+1] <- (temp[1,8]-temp[1,7]+1)/temp[1,13]
        sim[length(sim)+1] <- temp[,3]
        names[length(names)+1] <- as.character(j)
        #mismatch[length(mismatch)+1] <- temp[,5]
      }
      list[length(list)+1] <- as.character(i)
      Qlen[length(Qlen)+1] <- temp[1,13]
    }
    
    sortie <- cbind.data.frame(list,names,Qlen,bpondj,epondj,cov,sim)
    # To keep only the best blast corresponding to thresholds criteria
    if(dim(sortie[sortie$cov > COV_TRESH & sortie$sim > SIM_2nd_THRESH,])[1]==1) {
      sortiefiltered <- sortie[sortie$cov > COV_TRESH & sortie$sim > SIM_BEST_THRESH,]
      sortiefull <- rbind(sortiefull,sortiefiltered)
    }
    # Printing progression
    # if(t%/%100 == t/100) print(paste(round(100*t/length(query)),"%",sep=" "))
    t <- t+1
    #re-initialization
    epondj <- c()
    bpondj <- c()
    sim <- c()
    names <- c()
    # mismatch <- c()
    cov <- c()
    list <- c()
    Qlen <- c()
  }
  names(sortiefull) <- c("Query","Subject","Query_length","BitScore_corrected","Evalue","Coverage","Similarity")
  write.table(x=sortiefull,file=OUTPUT,quote = F,row.names = F,sep="\t")
}




######################################## #
# Searching for best reciprocal blast ####
######################################## #

# Simply by merging the two file by query and subject

FOCAL <- list.files(pattern = "focal2hordeum.bestblast.txt")
HORDEUM <- list.files(pattern = "hordeum2focal.bestblast.txt")

Nfile <- length(FOCAL)

for(i in 1:Nfile) {
  focal <- read.table(FOCAL[i],header=T)
  hordeum <- read.table(HORDEUM[i],header=T)
  reciprocal <- merge(focal,hordeum,by.x=c("Query","Subject"),by.y=c("Subject","Query"))
  output <- paste(strsplit(FOCAL[i],".fasta_focal2hordeum.bestblast.txt"),"reciprocalbestblast.txt",sep = "")
  write.table(x=reciprocal,file=output,quote = F,row.names = F,sep="\t")
  rm(reciprocal)
  }

