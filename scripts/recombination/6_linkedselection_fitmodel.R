# Script to fut the linked selection model
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)


# Option to choose to filter the data
FILTER <- 1000 #
size <- 1 # Size en cM



# Functions ####


# Background selection coefficients
BScoef <- function(chrom,pos,sel,fis) {
  posmax <- max(data_hordeumcM[which(data_hordeumcM$Chromosome==paste("chr",chrom,"H",sep="")),]$poscM)
  pos <- min(pos,posmax) # To avoid that the rounded position be higher than chromosome length
  r <- ifelse(data_hordeumcM$Chromosome==paste("chr",chrom,"H",sep=""),
              (1 - exp(-2*abs(data_hordeumcM$poscM-pos)))/2,
              1/2)
  r <- r*(1-fis)
  sum(data_hordeumcM$GenDens*sel/(r+sel)^2)
}
# Selective sweep coefficients
SWcoef <- function(chrom,pos,tau,fis) {
  posmax <- max(data_hordeumcM[which(data_hordeumcM$Chromosome==paste("chr",chrom,"H",sep="")),]$poscM)
  pos <- min(pos,posmax) # To avoid that the rounded position be higher than chromosome length
  r <- ifelse(data_hordeumcM$Chromosome==paste("chr",chrom,"H",sep=""),
              (1 - exp(-2*abs(data_hordeumcM$poscM-pos)))/2,
              1/2)
  r <- r*(1-fis)
  sum(data_hordeumcM$GenDens*exp(-r*tau))
}

# Juke-Cantor transformation to transform theta=4Neu into piS
JC <- function(x) -3*(-1+exp(-4*x/3))/4

# Set a limit when exact 0 is not allowed
ZERO <- 10^(-50)

# Likelihood function, Fis estimated
lnL <- function(par) {
  pi0 <- 10^par[1]
  udel <- 10^par[2]
  seldel <- 10^par[3]
  usweep <- 10^par[4]
  tau <- 10^par[5]
  fis <- par[6]
  #fis <- fis_list[i]
  mydatacM$BS <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = seldel,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
  mydatacM$SW <- mapply(function(x,y) SWcoef(chrom = x,pos = y,tau = tau,fis = fis),x=mydatacM$Chromosome,y=mydatacM$poscM)
  x <- mydatacM$reldiv*pi0/(exp(udel*mydatacM$BS) + usweep*mydatacM$SW)
  #x <- mydatacM$reldiv*pi0/(exp(udel*mydatacM$BS))
  x <- JC(x) + ZERO
  k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
  n <- mydatacM$nb_complete_site
  lik <- sum(dbinom(x = k,size = n,prob = x,log = T))
  return(-lik)
}

# Likelihood function, Fis fixed
lnLfixfis <- function(par) {
  pi0 <- 10^par[1]
  udel <- 10^par[2]
  seldel <- 10^par[3]
  usweep <- 10^par[4]
  tau <- 10^par[5]
  #fis <- fis_list[i]
  mydatacM$BS <- mapply(function(x,y) BScoef(chrom = x,pos = y,sel = seldel,fis = fisinit),x=mydatacM$Chromosome,y=mydatacM$poscM)
  mydatacM$SW <- mapply(function(x,y) SWcoef(chrom = x,pos = y,tau = tau,fis = fisinit),x=mydatacM$Chromosome,y=mydatacM$poscM)
  x <- mydatacM$reldiv*pi0/(exp(udel*mydatacM$BS) + usweep*mydatacM$SW)
  #x <- mydatacM$reldiv*pi0/(exp(udel*mydatacM$BS))
  x <- JC(x) + ZERO
  k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
  n <- mydatacM$nb_complete_site
  lik <- sum(dbinom(x = k,size = n,prob = x,log = T))
  return(-lik)
}


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


output <- paste("outputs/recombination/Results_BSfit_filter",FILTER,"_window",size,"cM.txt",sep="")

varnames <- c("Species","Windows","piS_obs","piS_est","Fis_obs","Fis_est","fdel","s_del","usweep","tau","lnLNull","lnLMax","lnLSat","R2deviance","R2efron")
write(x = varnames,file = output,ncolumns = length(varnames),append = F)

for(i in c(1:length(species_list))) {

i <- 1
  
WINDOW <- paste(size,"cM",sep = "")
data_hordeumcM <- read.table(paste("outputs/recombination/data_hordeum",WINDOW,".txt",sep=""),header = T)
  

SPECIES <- species_list[i]
FILE <- paste("outputs/recombination/",SPECIES,WINDOW,".txt",sep="")
mydatacM <- read.table(FILE,header = T)
mydatacM <- mydatacM[!is.na(mydatacM$piSyn),]
mydatacM$reldiv <- mydatacM$div/mean(mydatacM$div)

# Filtering out windows with too few positions
mydatacM <- mydatacM[mydatacM$nb_complete_site>FILTER,]

# Initialization of paramters for optimization
fisinit <- fis_list[i]
pi0init <- log10(2*mean(mydatacM$piSyn)/(1-fisinit))
udelinit <- log10(10^(-11))
seldelinit <- log10(0.1*(1+fisinit))
usweepinit <- log10(10^(-12))
tauinit <- log10(100)
init <- c(pi0init,udelinit,seldelinit,usweepinit,tauinit,fisinit)
# Boundaries for optimization
inf <- c(-3,-15,-5,-15,0,0)
sup <- c(3,-9,0,-9,10,0.999)
scaling <- c(abs(init[1:5]),1)
#lnL(init)
# Optimization of the likelihood function
ml <- optim(init,lnL,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling,maxit=1000,factr=10^7,lmm=20,trace=0))
#Output
minLik <- ml$value
pi_est <- 10^ml$par[1]
udel_est <- 10^ml$par[2]
sdel_est <- 10^(ml$par[3])
usweep_est <- 10^(ml$par[4])
tau_est <- 10^(ml$par[5])
f_est <- ml$par[6]
k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
n <- mydatacM$nb_complete_site
likNull <- sum(dbinom(x = k,size = n,prob = mean(k/n),log = T))
likSat <- sum(dbinom(x = k,size = n,prob = k/n,log = T))
# Deviance R squared
R2deviance <- 1 - (likSat + minLik)/(likSat - likNull)
# Efron R squared
mydatacM$pred <-  pi_est*mydatacM$reldiv /(exp(udel_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = sdel_est,fis = f_est),x=mydatacM$Chromosome,y=mydatacM$poscM))
                                   + usweep_est*mapply(function(x,y) SWcoef(chrom = x,pos = y,tau = tau_est,fis = f_est),x=mydatacM$Chromosome,y=mydatacM$poscM))
R2efron <- summary(lm(mydatacM$piSyn~mydatacM$pred))$r.squared
result <- c(SPECIES,WINDOW,mean(mydatacM$piSyn),pi_est,fisinit,f_est,udel_est,sdel_est,usweep_est,tau_est,likNull,-minLik,likSat,R2deviance,R2efron)
write(x = result,file = output,ncolumns = length(result),append = T)

print(c(SPECIES,mean(mydatacM$piSyn),pi_est,fisinit,round(f_est,3),round(R2deviance,3)))

}



# Estimation Fixed Fis ####

output <- paste("outputs/recombination/Results_BSfit_FixedFis_filter",FILTER,"_window",size,"cM.txt",sep="")
varnames <- c("Species","Windows","piS_obs","piS_est","Fis_obs","fdel","s_del","usweep","tau","lnLNull","lnLMax","lnLSat","R2deviance","R2efron")
write(x = varnames,file = output,ncolumns = length(varnames),append = F)

for(i in c(1:length(species_list))) {
  
  WINDOW <- paste(size,"cM",sep = "")
  data_hordeumcM <- read.table(paste("outputs/recombination/data_hordeum",WINDOW,".txt",sep=""),header = T)
  SPECIES <- species_list[i]
  FILE <- paste(SPECIES,WINDOW,".txt",sep="")
  mydatacM <- read.table(FILE,header = T)
  mydatacM$reldiv <- mydatacM$div/mean(mydatacM$div)

  # Filtering out windows with too few positions
    mydatacM <- mydatacM[mydatacM$nb_complete_site>FILTER,]

# Initialization of paramters for optimization
  fisinit <- fis_list[i]
  pi0init <- log10(2*mean(mydatacM$piSyn)/(1-fisinit))
  udelinit <- log10(10^(-8))
  seldelinit <- log10(0.1*(1+fisinit))
  usweepinit <- log10(10^(-9))
  tauinit <- log10(100)
  init <- c(pi0init,udelinit,seldelinit,usweepinit,tauinit)
# Boundaries for optimization
  inf <- c(-3,-12,-5,-12,0)
  sup <- c(3,-6,0,-6,10)
# Optimization of the likelihood function
  ml <- optim(init,lnLfixfis,lower=inf,upper=sup,method="L-BFGS-B",control=list(maxit=1000,factr=10^7,lmm=20,trace=0))
# Outputs
  minLik <- ml$value
  pi_est <- 10^ml$par[1]
  udel_est <- 10^ml$par[2]
  sdel_est <- 10^(ml$par[3])
  usweep_est <- 10^(ml$par[4])
  tau_est <- 10^(ml$par[5])
  k <- round(mydatacM$piSyn*mydatacM$nb_complete_site,0)
  n <- mydatacM$nb_complete_site
  likNull <- sum(dbinom(x = k,size = n,prob = mean(k/n),log = T))
  likSat <- sum(dbinom(x = k,size = n,prob = k/n,log = T))
  # Deviance R squared
  R2deviance <- 1 - (likSat + minLik)/(likSat - likNull)
  # Efron R squared
  mydatacM$pred <-  pi_est*mydatacM$reldiv /(exp(udel_est*mapply(function(x,y) BScoef(chrom = x,pos = y,sel = sdel_est,fis = fisinit),x=mydatacM$Chromosome,y=mydatacM$poscM))
                                             + usweep_est*mapply(function(x,y) SWcoef(chrom = x,pos = y,tau = tau_est,fis = fisinit),x=mydatacM$Chromosome,y=mydatacM$poscM))
  R2efron <- summary(lm(mydatacM$piSyn~mydatacM$pred))$r.squared
  
  result <- c(SPECIES,WINDOW,mean(mydatacM$piSyn),pi_est,fisinit,udel_est,sdel_est,usweep_est,tau_est,likNull,-minLik,likSat,R2deviance,R2efron)
  write(x = result,file = output,ncolumns = length(result),append = T)
  
  print(c(SPECIES,mean(mydatacM$piSyn),pi_est,fisinit,round(R2deviance,3)))
  
}



