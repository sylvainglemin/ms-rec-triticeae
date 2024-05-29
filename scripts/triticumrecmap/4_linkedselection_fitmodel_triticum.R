# Script to fut the linked selection model
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated March 2024

# Version with only BS (no sweep) and another way of modelling:
# Several categories of selection (instead of one) but with fixed s
# The total mutation rate for each category is estimated

if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)
if (!require("parallel"))   install.packages("parallel", dependencies = TRUE)

PATH <- NULL
# PATH <- "/scratch/sglemin/ms-rec-triticeae/"

# Option to choose to filter the data
FILTER <- 300 # Remove windows with less than FILTER complete position
size <- 1 # Size en cM
# Option to choose for bootsptrap and parallelisation
Nboot <- 100
Ncores <- 96

############## #
# Functions ####
############## #

# Background selection coefficients
BScoef <- function(chrom,pos,sel,fis) {
  posmax <- max(data_triticumcM[which(data_triticumcM$Chromosome==chrom),]$poscM)
  pos <- min(pos,posmax) # To avoid that the rounded position be higher than chromosome length
  r <- ifelse(data_triticumcM$Chromosome==chrom,
              (1 - exp(-2*abs(data_triticumcM$poscM-pos)))/2,
              1/2)
  r <- r*(1-fis)
  sum(data_triticumcM$GenDens*sel/(r+sel)^2)
}

# Juke-Cantor transformation to transform theta=4Neu into piS
JC <- function(x) -3*(-1+exp(-4*x/3))/4

# Set a limit when exact 0 is not allowed
ZERO <- 10^(-50)





# Likelihood function, Fis fixed
# lnLfixfis <- function(par) {
#   pi0 <- 10^par[1]
#   u1 <- 10^par[2]
#   u2 <- 10^(par[2]+par[3])
#   u3 <- 10^(par[2]+par[3]+par[4])
#   s1 <- 10^par[5]
#   s2 <- 10^par[6]
#   s3 <- 10^par[7]
#   #fis <- fis_list[i]
#   mydatacM$BS1 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s1,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#   mydatacM$BS2 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s2,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#   mydatacM$BS3 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s3,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#   theta <- mydatacM$reldiv*pi0*mydatacM$nb_complete_site*exp(-u1*mydatacM$BS1-u2*mydatacM$BS2-u3*mydatacM$BS3) + ZERO
#   k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
#   lik <- sum(k*log(theta)-(k+1)*log(theta+1))
#   return(-lik)
# }



# Function to run the optimization, fis free

run_fit_free <- function(mydatacM,fisinit,init=NULL,PRINT=F) {
  #fisinit <- fis_list[i]
  # Likelihood function, Fis estimated
  lnL <- function(par) {
    pi0 <- 10^par[1]
    u1 <- 10^par[2]
    u2 <- 10^(par[2]+par[3])
    u3 <- 10^(par[2]+par[3]+par[4])
    s1 <- 10^par[5]
    s2 <- 10^par[6]
    s3 <- 10^par[7]
    fis <- par[8]
    #fis <- fis_list[i]
    mydatacM$BS1 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s1,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
    mydatacM$BS2 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s2,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
    mydatacM$BS3 <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s3,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
    theta <- pi0*mydatacM$nb_complete_site*exp(-u1*mydatacM$BS1-u2*mydatacM$BS2-u3*mydatacM$BS3) + ZERO
    k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
    lik <- sum(k*log(theta)-(k+1)*log(theta+1))
    return(-lik)
  }
  if(is.null(init)) {
    pi0init <- log10(2*mean(mydatacM$piSyn)/(1-fisinit))
    u1init <- log10(10^(-9))
    addu2init <- -0.1
    addu3init <- -0.1
    s1init <- log10(10^(-1))
    s2init <- log10(10^(-2))
    s3init <- log10(10^(-3))
    init <- c(pi0init,u1init,addu2init,addu3init,s1init,s2init,s3init,fisinit) 
  }
  # Boundaries for optimization
  inf <- c(-3,-13,-3,-3,-4,-4,-4,0)
  sup <- c(3,-7,0,0,0,0,0,0.99)
  scaling <- c(abs(init[1:7]),1)
  #lnL(init)
  # Optimization of the likelihood function
  ml <- optim(init,lnL,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling,maxit=100,factr=10^7,lmm=5,trace=0))
  #Output
  minLik <- ml$value
  pi_est <- 10^ml$par[1]
  u1_est <- 10^ml$par[2]
  u2_est <- 10^(ml$par[2]+ml$par[3])
  u3_est <- 10^(ml$par[2]+ml$par[3]+ml$par[4])
  s1_est <- 10^ml$par[5]
  s2_est <- 10^ml$par[6]
  s3_est <- 10^ml$par[7]
  f_est <- ml$par[8]
  # Deviance R squared
  k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0) + ZERO
  likNull <- sum(k*log(mean(k))-(k+1)*log(mean(k)+1))
  likSat <- sum(k*log(k)-(k+1)*log(k+1))
  R2deviance <- 1 - (likSat + minLik)/(likSat - likNull)
  # Efron R squared
  mydatacM$pred <-  pi_est*exp(
    -u1_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s1_est,fis = f_est),x=mydatacM$Chromosome,y=mydatacM$poscM)
    -u2_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s2_est,fis = f_est),x=mydatacM$Chromosome,y=mydatacM$poscM)
    -u3_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s3_est,fis = f_est),x=mydatacM$Chromosome,y=mydatacM$poscM)
  )
  R2efron <- summary(lm(mydatacM$piSyn~mydatacM$pred))$r.squared
  res <- c(SPECIES,GENOME,WINDOW,mean(mydatacM$piSyn),pi_est,fisinit,f_est,u1_est,s1_est,u2_est,s2_est,u3_est,s3_est,likNull,-minLik,likSat,R2deviance,R2efron)
  names(res) <- c("Species","Genome","Windows","piS_obs","piS_est","Fis_obs","Fis_est","u1_est","s1_est","u2_est","s2_est","u3_est","s3_est","lnLNull","lnLMax","lnLSat","R2deviance","R2efron")
  # Control plots
  if(PRINT){
    G <- ggplot(data=mydatacM,aes(x=Start,y=piSyn)) + geom_point() + geom_line(aes(x=Start,y=pred,col="red"))
    G <- G +  facet_wrap(~Chromosome,scales = "free_x")
    G <- G + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
    G <- G + theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic"))
    G <- G + xlab("Chromosome position (in bp)") + ylab(expression(pi[S])) 
    ggsave(filename = paste(PATH,"figures/additional/triticum_recmap/",SPECIES,"_BSfit_genome",GENOME,"_",WINDOW,"filter",FILTER,".pdf",sep = ""),plot = G)
    G <- ggplot(data=mydatacM,aes(x=pred,y=piSyn)) + geom_point()
    G <- G + scale_x_log10() + scale_y_log10()
    G <- G + geom_smooth(method = "lm")
    G <- G + xlab(expression(pi[Pred])) + ylab(expression(pi[Obs])) + ggtitle(SPECIES)
    G <- G + theme(title = element_text(size = 14),axis.title = element_text(size = 18), axis.text = element_text(size = 14))
    ggsave(filename = paste(PATH,"figures/additional/triticum_recmap/",SPECIES,"_GoF_genome",GENOME,"_",WINDOW,"filter",FILTER,".pdf",sep = ""),plot = G)
  }

  return(res)
}

################ #
# Estimations ####
################ #

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
fis_list <- c(0.9,
              0.6,
              0.8,
              0.7,
              0.0,
              0.9,
              0.7,
              0.0,
              0.9,
              0.8,
              0.9,
              0.85,
              0.95
)

# Estimation Free Fis ####


for(GENOME in c("A","B","D")){
  # Loop over the 13 species
  # GENOME <- "A"
  output <- paste(PATH,"outputs/recombination/triticum/Results_BSfit_genome",GENOME,"_filter",FILTER,"_window",size,"cM.txt",sep="")
  varnames <- c("Species","Genome","Windows","piS_obs","piS_est","Fis_obs","Fis_est","u1_est","s1_est","u2_est","s2_est","u3_est","s3_est","lnLNull","lnLMax","lnLSat","R2deviance","R2efron")
  write(x = varnames,file = output,ncolumns = length(varnames),append = F)
  
  
  for(i in c(1:length(species_list))) {
    # i <- 5
    
    #Preparing the file
    WINDOW <- paste(size,"cM",sep = "")
    data_triticumcM <- read.table(paste(PATH,"outputs/recombination/triticum/Triticum_genome",GENOME,"_",WINDOW,".txt",sep=""),header = T)
    SPECIES <- species_list[i]
    FILE <- paste(PATH,"outputs/recombination/triticum/",SPECIES,"_genome",GENOME,"_",WINDOW,".txt",sep="")
    mydatacM <- read.table(FILE,header = T)
    mydatacM <- mydatacM[!is.na(mydatacM$piSyn),]
    # Filtering out windows with too few positions
    mydatacM <- mydatacM[mydatacM$nb_complete_site>FILTER,]

    # Runing the model on the orginal data
    result <- run_fit_free(mydatacM = mydatacM,fisinit=fis_list[i],PRINT=T)

    # Exporting the results
    write(x = result,file = output,ncolumns = length(result),append = T)
    
    # Bootstrap
    if(Nboot > 1) {
      n_sample <- dim(mydatacM)[1]
      data_boot <- replicate(n = Nboot,
                             expr = mydatacM[sample(n_sample,n_sample,T),],
                             simplify = F)
      result <- as.list(result)
      names(result) <- varnames
      init_boot <- c(
        log10(as.numeric(result$piS_est)),
        log10(as.numeric(result$u1_est)),
        log10(as.numeric(result$u2_est))-log10(as.numeric(result$u1_est)),
        log10(as.numeric(result$u3_est))-log10(as.numeric(result$u2_est)),
        log10(as.numeric(result$s1_est)),
        log10(as.numeric(result$s2_est)),
        log10(as.numeric(result$s3_est)),
        as.numeric(result$Fis_est)
      )
      ww <- mclapply(
        X = data_boot,
        FUN = function(df) run_fit_free(mydatacM = df,fisinit=fis_list[i],init = init_boot),
        mc.preschedule = T,
        mc.cores = Ncores
      )
      result_boot <- data.frame(matrix(unlist(ww), nrow=length(ww), byrow=TRUE))
      names(result_boot) <- varnames
      out_boot <- paste(PATH,"outputs/recombination/triticum/Bootstrap_",SPECIES,"_genome",GENOME,"_filter",FILTER,"_window",size,"cM.txt",sep="")
      write.table(x = result_boot,file = out_boot,quote = F,row.names = F)
    }
    
    print(c(SPECIES,GENOME))
  }
}



# # Estimation Fixed Fis ####
# 
# output <- paste("outputs/recombination/hordeum/Results_BSfit_FixedFis_filter",FILTER,"_window",size,"cM.txt",sep="")
# varnames <- c("Species","Windows","piS_obs","piS_est","Fis_obs","u1_est","s1_est","u2_est","s2_est","u3_est","s3_est","lnLNull","lnLMax","lnLSat","R2deviance","R2efron")
# write(x = data,file = output,ncolumns = length(varnames),append = F)
# 
# for(i in c(1:length(species_list))) {
# 
# # i <- 1
#   WINDOW <- paste(size,"cM",sep = "")
#   data_hordeumcM <- read.table(paste("outputs/recombination/hordeum/data_hordeum",WINDOW,".txt",sep=""),header = T)
#   
#   SPECIES <- species_list[i]
#   FILE <- paste("outputs/recombination/hordeum/",SPECIES,WINDOW,".txt",sep="")
#   mydatacM <- read.table(FILE,header = T)
#   mydatacM <- mydatacM[!is.na(mydatacM$piSyn),]
#   mydatacM$reldiv <- mydatacM$div/mean(mydatacM$div)
#   
#   # Filtering out windows with too few positions
#   mydatacM <- mydatacM[mydatacM$nb_complete_site>FILTER,]
#   
#   # Initialization of parameters for optimization
#   fis <- fis_list[i]
#   pi0init <- log10(2*mean(mydatacM$piSyn)/(1-fis))
#   u1init <- log10(10^(-9))
#   addu2init <- -0.1
#   addu3init <- -0.1
#   s1init <- log10(10^(-1))
#   s2init <- log10(10^(-2))
#   s3init <- log10(10^(-3))
#   init <- c(pi0init,u1init,addu2init,addu3init,s1init,s2init,s3init)
#   # Boundaries for optimization
#   inf <- c(-3,-13,-3,-3,-4,-4,-4)
#   sup <- c(3,-7,0,0,0,0,0)
#   scaling <- abs(init[1:7])
#   lnLfixfis(init)
#   # Optimization of the likelihood function
#   ml <- optim(init,lnLfixfis,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling,maxit=1000,factr=10^7,lmm=20,trace=0))
#   #Output
#   minLik <- ml$value
#   pi_est <- 10^ml$par[1]
#   u1_est <- 10^ml$par[2]
#   u2_est <- 10^(ml$par[2]+ml$par[3])
#   u3_est <- 10^(ml$par[2]+ml$par[3]+ml$par[4])
#   s1_est <- 10^ml$par[5]
#   s2_est <- 10^ml$par[6]
#   s3_est <- 10^ml$par[7]
#   # Deviance R squared
#   k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0) + ZERO
#   likNull <- sum(k*log(mean(k))-(k+1)*log(mean(k)+1))
#   likSat <- sum(k*log(k)-(k+1)*log(k+1))
#   R2deviance <- 1 - (likSat + minLik)/(likSat - likNull)
#   # Efron R squared
#   mydatacM$pred <-  pi_est*mydatacM$reldiv*exp(
#     -u1_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s1_est,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#     -u2_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s2_est,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#     -u2_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = s3_est,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
#   )
#   R2efron <- summary(lm(mydatacM$piSyn~mydatacM$pred))$r.squared
#   
#   G <- ggplot(data=mydatacM,aes(x=Start,y=piSyn)) + geom_point() + geom_line(aes(x=Start,y=pred,col="red"))
#   G <- G +  facet_wrap(~Chromosome,scales = "free_x")
#   G <- G + theme(panel.background = element_rect(fill = "grey95", colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5))
#   G <- G + theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic"))
#   G <- G + xlab("Chromosome position (in bp)") + ylab(expression(pi[S])) 
#   ggsave(filename = paste("figures/additional/",SPECIES,"_BSfit_",WINDOW,"filter",FILTER,".pdf",sep = ""),plot = G)
#   
#   G <- ggplot(data=mydatacM,aes(x=pred,y=piSyn)) + geom_point()
#   G <- G + scale_x_log10() + scale_y_log10()
#   G <- G + geom_smooth(method = "lm")
#   G <- G + xlab(expression(pi[Pred])) + ylab(expression(pi[Obs])) + ggtitle(SPECIES)
#   G <- G + theme(title = element_text(size = 14),axis.title = element_text(size = 18), axis.text = element_text(size = 14))
#   ggsave(filename = paste("figures/additional/",SPECIES,"_GoF_",WINDOW,"filter",FILTER,".pdf",sep = ""),plot = G)
#   
#   
#   result <- c(SPECIES,WINDOW,mean(mydatacM$piSyn),pi_est,fis,u1_est,s1_est,u2_est,s2_est,u3_est,s3_est,likNull,-minLik,likSat,R2deviance,R2efron)
#   write(x = result,file = output,ncolumns = length(result),append = T)
#   
#   print(c(SPECIES,mean(mydatacM$piSyn),pi_est,fis,round(R2deviance,3)))
#   
# }


