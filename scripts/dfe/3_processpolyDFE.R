library(ggplot2)

# Require the following R script from the polyDFE package
source("postprocessing.R")

require(cubature)


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


#### integration over a DFE to obtain the proportion of mutations within a given range


ZERO <- 10^(-10) ## something almost zero to avoid numerical divergence in 0

# The DFE is composed of a negative gamma with Smean (< 0) and shape (> 0) in proportion 1-p
# and a positive exponential with parameter Spos in proportion p


phi_neg <- function(Smean,shape){
  g <- function(S){
    exp((-S)*(shape)/Smean)*(-S)^(shape-1)*(-Smean/shape)^(-shape)/gamma(shape)
  }
  return(g)
}

phi_pos <- function(Spos){
  g <- function(S){
    exp((-S/Spos))/Spos
  }
  return(g)
}

phi_disp <- function(Smean,shape,Smax){
  g <- function(S){
    exp((-S+Smax)*(shape)/Smean)*(-S+Smax)^(shape-1)*(-Smean/shape)^(-shape)/gamma(shape)
  }
  return(g)
}


classMutation <- function(Smean,shape,Spos,p,Smin,Smax){
  if(abs(Smin)< ZERO) Smin <- ZERO
  if(abs(Smax)<ZERO) Smax <- -ZERO
  if(Smin < ZERO){
    if(Smax <= -ZERO){
      Cneg <- adaptIntegrate(phi_neg(Smean,shape), lowerLimit=Smin, upperLimit=Smax)$integral
      return((1-p)*Cneg)
    }
    else{
      Cneg <- adaptIntegrate(phi_neg(Smean,shape), lowerLimit=Smin, upperLimit=-ZERO)$integral
      Cpos <- adaptIntegrate(phi_pos(Spos), lowerLimit=ZERO, upperLimit=Smax)$integral
      return((1-p)*Cneg+p*Cpos)
    }
  }
  else{
    Cpos <- adaptIntegrate(phi_pos(Spos), lowerLimit=Smin, upperLimit=Smax)$integral
    return(p*Cpos)
  }
}


classMutationDisp <- function(Smean,shape,Smax,Slow,Sup){
  C <- adaptIntegrate(phi_disp(Smean,shape,Smax), lowerLimit=Slow, upperLimit=Sup)$integral
  return(C)
}


w_neg <- function(Smean,shape) {
  g <- function(S) {phi_neg(Smean,shape)(S) * S / (1 - exp(-S))}
  return( g )
}

w_pos <- function(Spos) {
  g <- function(S) {phi_pos(Spos)(S) * S / (1 - exp(-S))}
  return( g )
}

omega <- function(Smean,shape,Spos,p,Smin,Smax){
  if(abs(Smin) < ZERO) Smin <- ZERO
  if(abs(Smax) < ZERO) Smax <- -ZERO
  if(Smin < ZERO){
    if(Smax <= -ZERO){
      wneg <- adaptIntegrate(w_neg(Smean,shape), lowerLimit=Smin, upperLimit=Smax)$integral
      return((1-p)*wneg)
    }
    else{
      wneg <- adaptIntegrate(w_neg(Smean,shape), lowerLimit=Smin, upperLimit=-ZERO)$integral
      wpos <- adaptIntegrate(w_pos(Spos), lowerLimit=ZERO, upperLimit=Smax)$integral
      return((1-p)*wneg+p*wpos)
    }
  }
  else{
    wpos <- adaptIntegrate(w_pos(Spos), lowerLimit=Smin, upperLimit=Smax)$integral
    return(p*wpos)
  }
}




##################### #
# SUMARIZE RESULTS ####
##################### #

# Generate the model averaged estimates for all boostraped datasets for the four focal species and the four data categories (total, mapped, lowre, highrec)
# Create a single data.frame with all data

list_species <- c("Ae_mutica","Ae_speltoides","Ae_tauschii","T_urartu")

dfe_species <- c()
Vnames <- c()

for(SPECIES in list_species) {
  
  setwd(paste("/Users/sylvain/Documents/Boulot/Recherche/Thematiques/SystRepro/TRANS/Triticees/DFE/Old/PolyDFEResults/",SPECIES,"/",sep=""))
  list_type <- c("total","mapped","highrec","lowrec")
  
  for(TYPE in list_type) {
  
    dof <- c(3,2,1,0) #
    Sd_weighted <- vector("numeric")
    b_weighted <- vector("numeric")
    Sb_weighted <- vector("numeric")
    pb_weighted <- vector("numeric")

    for(i in 1:1000){
      PREFIX <- paste(SPECIES,"_",TYPE,"_polyDFE_boot",i,sep="")
      aic <- vector("numeric")
      Sd <- vector("numeric")
      beta <- vector("numeric")
      Sb <- vector("numeric")
      pb <- vector("numeric")
      for(j in c(1:4)){
        FILE <- paste(PREFIX,"_modelC",j,"_res.txt",sep="")
        if(file.exists(FILE)){
          if(length(grep("alpha_dfe =",readLines(FILE)))>0){
            est <- parseOutput(FILE)
            aic <- c(aic,2*dof[j] - 2*est[[1]]$lk)
            Sd <- c(Sd,est[[1]]$values["S_d"])
            beta <- c(beta,est[[1]]$values["b"])
            Sb <- c(Sb,est[[1]]$values["S_b"])
            pb <- c(pb,est[[1]]$values["p_b"]) 
          }
        }
      }
      if(length(aic)==4){
        weight <- exp(-(aic-min(aic))/2)
        #weight <- c(0,0,1,0)
        Sd_weighted <- c(Sd_weighted,sum(Sd*weight)/sum(weight))
        b_weighted <- c(b_weighted,sum(beta*weight)/sum(weight))
        Sb_weighted <- c(Sb_weighted,sum(Sb*weight)/sum(weight))
        pb_weighted <- c(pb_weighted,sum(pb*weight)/sum(weight))
      }
    }
    Cweak <- mapply(function(w,x,y,z) classMutation(w,x,y,z,-10,-1),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    CneutralNeg <- mapply(function(w,x,y,z) classMutation(w,x,y,z,-1,0),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    CneutralPos <- mapply(function(w,x,y,z) classMutation(w,x,y,z,0,1),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    Cpos <- mapply(function(w,x,y,z) classMutation(w,x,y,z,0,1000),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    CposStrong <- mapply(function(w,x,y,z) classMutation(w,x,y,z,1,1000),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    Cdel <- 1 - Cweak - CneutralNeg -CneutralPos - Cpos 
    wa <- mapply(function(w,x,y,z) omega(w,x,y,z,0,1000*y),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    wd <- mapply(function(w,x,y,z) omega(w,x,y,z,-100000,0),Sd_weighted,b_weighted,Sb_weighted,pb_weighted)
    alpha <- wa/(wa+wd)
    species <- rep(SPECIES,length(Cdel))
    one_sp <- cbind(species,TYPE,Sd_weighted,b_weighted,Sb_weighted,pb_weighted,Cdel,Cweak,CneutralNeg,CneutralPos,Cpos,CposStrong,wa,wd,alpha)
    dfe_species <- rbind(dfe_species,one_sp)
    print(c(SPECIES,TYPE))
  }
}

mydata <- as.data.frame(dfe_species)
names(mydata) <- c("Species","Type","Sd","b","Sb","pb","Cdel","Cweak","CneutralNeg","CneutralPos","Cpos","CposStrong","wa","wd","alpha")
mydata$Sd <- as.numeric(as.character(mydata$Sd))
mydata$b <- as.numeric(as.character(mydata$b))
mydata$Sb <- as.numeric(as.character(mydata$Sb))
mydata$pb <- as.numeric(as.character(mydata$pb))
mydata$Cdel <- as.numeric(as.character(mydata$Cdel))
mydata$Cweak <- as.numeric(as.character(mydata$Cweak))
mydata$CneutralNeg <- as.numeric(as.character(mydata$CneutralNeg))
mydata$CneutralPos <- as.numeric(as.character(mydata$CneutralPos))
mydata$Cpos <- as.numeric(as.character(mydata$Cpos))
mydata$CposStrong <- as.numeric(as.character(mydata$CposStrong))
mydata$wa <- as.numeric(as.character(mydata$wa))
mydata$wd <- as.numeric(as.character(mydata$wd))
mydata$alpha <- as.numeric(as.character(mydata$alpha))

mydata$Rec <- as.numeric(sapply(mydata$Type,function(x) strsplit(as.character(x),"_")[[1]][2]))
mydata$MS <- as.factor(ifelse(mydata$Species=="Ae_mutica","Outcrossing",ifelse(mydata$Species=="Ae_speltoides","Outcrossing","Selfing")))

mydata <- droplevels(mydata)
                    
write.table(mydata,"dfe_results_4focal.txt",sep = "\t",quote = F,row.names = F)


############################# #
# ANALYSE AND PLOT RESULTS ####
############################# #


mydata <- read.table("dfe_results_4focal.txt",header=T)

dataplot <- mydata[is.na(mydata$Rec) & mydata$Type!="mapped",]
dataplot$Species <- ifelse(dataplot$Species=="Ae_mutica","Ae. mutica",dataplot$Species)
dataplot$Species <- ifelse(dataplot$Species=="Ae_speltoides","Ae. speltoides",dataplot$Species)
dataplot$Species <- ifelse(dataplot$Species=="Ae_tauschii","Ae. tauschii",dataplot$Species)
dataplot$Species <- ifelse(dataplot$Species=="T_urartu","T. urartu",dataplot$Species)


G1 <- ggplot(data = dataplot,aes(x = Type,y=(Cweak+CneutralNeg)/(1-Cpos),fill=Type)) + geom_boxplot(outlier.shape = NA)
G1 <- G1 + facet_grid(~Species)
G1 <- G1 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G1 <- G1 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G1 <- G1 + xlab(expression("% of weakly deleterious mutations (-10 < 4"*italic(N[e])*italic(s)*" < 0)")) + ylab("")
G1


G2 <- ggplot(data = dataplot,aes(x = Type,y=(Cdel+CneutralPos)/(1-Cpos),fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0.78,0.91)
G2 <- G2 + facet_grid(~Species)
G2 <- G2 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G2 <- G2 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G2 <- G2 + xlab("% of strongly deleterious mutations (2Ns < -10)") + ylab("")
G2

G3 <- ggplot(data = dataplot,aes(x = Type,y=CneutralPos,fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0,0.15)
G3 <- G3 + facet_grid(~Species)
G3 <- G3 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G3 <- G3 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G3 <- G3 + xlab("% of weakly beneficial mutations (0 < 2Ns < 1)") + ylab("")
G3

G4 <- ggplot(data = dataplot,aes(x = Type,y=CposStrong,fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0,0.025)
G4 <- G4 + facet_grid(~Species)
G4 <- G4 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G4 <- G4 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G4 <- G4 + xlab("% of strongly beneficial mutations (2Ns > 1)") + ylab("")
G4

G5 <- ggplot(data = dataplot,aes(x = Type,y=Cpos,fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0,0.15)
G5 <- G5 + facet_grid(~Species)
G5 <- G5 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G5 <- G5 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G5 <- G5 + xlab("% of  beneficial mutations") + ylab("")
G5

G6 <- ggplot(data = dataplot,aes(x = Type,y=wa,fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0,0.25)
G6 <- G6 + facet_grid(~Species)
G6 <- G6 + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
G6 <- G6 + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
G6 <- G6 + xlab(expression("Adaptive substitution rate ("*omega[a]*")")) + ylab("")
G6

# Only graphs 1 and 6 are presented in the main text

pdf("DFE_Species_Recombination.pdf",width = 8,height = 12)
multiplot(G1,G6)
dev.off()



##################################### #
# POINT ESTIMATES FOR W/S results ####
##################################### #


dfeGC <- data.frame()

for(SPECIES in list_species) {
  
  FILE <- paste(SPECIES,"_output_test.txt",sep="")
  dfe <- c()
  est <- parseOutput(FILE)
  for(i in c(1:length(est))){
    param <- est[[i]]$values
    lnL <- est[[i]]$lk
    a <- est[[i]]$alpha
    Sd <- param["S_d"]
    b <- param["b"]
    Sb <- param["S_b"]
    pb <- param["p_b"]
    err <- param["eps_an"]
    dfe <- rbind(dfe,c(lnL,a,Sd,b,Sb,pb,err))
  }
  
  dfe <- data.frame(dfe)
  names(dfe) <- c("lnL","alpha","Sd","beta","Sb","pb","err")
  dfe$ATCG <- rep(c(rep("WWSS",4),rep("WS",4),rep("SW",4)),3)
  dfe$Type <- c(rep("total",12),rep("lowrec",12),rep("highrec",12))
  dfe$dataset <- interaction(dfe$ATCG,dfe$Type)
  dfe$Species <- rep(SPECIES,36)
  dfe$dof <- 3 -(is.na(dfe$pb/dfe$pb) + is.na(dfe$Sb/dfe$Sb) + is.na(dfe$err/dfe$err))
  dfe$AIC <- 2 * (dfe$dof - dfe$lnL)
  
  dfe$AICmin <- vector(length = length(dfe$dataset))
  for(i in levels(dfe$dataset)) {
    dfe[dfe$dataset==i,]$AICmin <- min(dfe[dfe$dataset==i,]$AIC)
  }
  dfe$weight <- exp(-(dfe$AIC-dfe$AICmin)/2)
  
  dfe$Sd_weighted <- vector(length = length(dfe$dataset))
  dfe$b_weighted <- vector(length = length(dfe$dataset))
  dfe$Sb_weighted <- vector(length = length(dfe$dataset))
  dfe$pb_weighted <- vector(length = length(dfe$dataset))
  dfe$omegaA <- vector(length = length(dfe$dataset))
  dfe$Cweak <- vector(length = length(dfe$dataset))
  for(i in levels(dfe$dataset)) {
    mydata <- dfe[dfe$dataset==i,]
    dfe[dfe$dataset==i,]$Sd_weighted <- sum(mydata$Sd*mydata$weight)/sum(mydata$weight)
    dfe[dfe$dataset==i,]$b_weighted <-  sum(mydata$beta*mydata$weight)/sum(mydata$weight)
    dfe[dfe$dataset==i,]$Sb_weighted <- sum(mydata$Sb*mydata$weight)/sum(mydata$weight)
    dfe[dfe$dataset==i,]$pb_weighted <- sum(mydata$pb*mydata$weight)/sum(mydata$weight)
    dfe[dfe$dataset==i,]$omegaA <- omega(dfe[dfe$dataset==i,]$Sd_weighted,
                                         dfe[dfe$dataset==i,]$b_weighted,
                                         dfe[dfe$dataset==i,]$Sb_weighted,
                                         dfe[dfe$dataset==i,]$pb_weighted,ZERO,1000)
    dfe[dfe$dataset==i,]$Cweak <- classMutation(dfe[dfe$dataset==i,]$Sd_weighted,
                                                dfe[dfe$dataset==i,]$b_weighted,
                                                dfe[dfe$dataset==i,]$Sb_weighted,
                                                dfe[dfe$dataset==i,]$pb_weighted,-10,0)
    rm(mydata)
  }
  dfeGC <- rbind(dfeGC,dfe)
}

dfeGC <- dfeGC[dfeGC$AIC==dfeGC$AICmin,]

dfeGC$Species <- ifelse(dfeGC$Species=="Ae_mutica","Ae. mutica",dfeGC$Species)
dfeGC$Species <- ifelse(dfeGC$Species=="Ae_speltoides","Ae. speltoides",dfeGC$Species)
dfeGC$Species <- ifelse(dfeGC$Species=="Ae_tauschii","Ae. tauschii",dfeGC$Species)
dfeGC$Species <- ifelse(dfeGC$Species=="T_urartu","T. urartu",dfeGC$Species)

Ggc <- ggplot(data = dataplot,aes(x = Type,y=wa,fill=Type)) + geom_boxplot(outlier.shape = NA) + ylim(0,0.5)
Ggc <- Ggc + geom_point(data = dfeGC[dfeGC$ATCG=="WWSS",],aes(x=Type,y=omegaA),pch="+",cex=5)
Ggc <- Ggc + geom_point(data = dfeGC[dfeGC$ATCG=="WS",],aes(x=Type,y=omegaA),pch=1,cex=2)
Ggc <- Ggc + geom_point(data = dfeGC[dfeGC$ATCG=="SW",],aes(x=Type,y=omegaA),pch=16,cex=2)
Ggc <- Ggc + facet_grid(~Species)
Ggc <- Ggc + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
Ggc <- Ggc + theme(title = element_text(size = 12),axis.title = element_text(size = 12),legend.position="none",strip.text = element_text(size = 12, face = "italic"))
Ggc <- Ggc + xlab(expression("Adaptive substitution rate ("*omega[a]*")")) + ylab("")
Ggc


pdf("DFE_Species_RecombinationGC.pdf",width = 8,height = 5)
Ggc
dev.off()

