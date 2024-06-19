# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Simulation of the sampling matrix of input parameters of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(deSolve)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(pracma)
library(sensitivity)
library(parallel)
library(ggplot2)

#Functions for running the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021)
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencATload.R",sep=""))
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencATout.R",sep=""))
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencAT.R",sep=""))

#Functions for computing Shapley effects
source(paste(getwd(),"/FunctionsShapleyeffects.R",sep=""))

# ---- Model characteristics ----

Inputs<-c("khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2",
          "k_br","p1","p2","p3","p4","p5","p6") #Input parameters
k<-length(Inputs) #Number of input parameters

# ---- Sampling of input parameters -----

# List of input parameters studied:

# khyd_ndf: Hydrolysis rate constant of neutral detergent fiber, 1/h 
# khyd_nsc: Hydrolysis rate constant of non-structural carbohydrates, 1/h  
# khyd_pro: Hydrolysis rate constant of proteins, 1/h 
# km_su: Maximum specific utilization rate constant of sugars, mol/(mol h) 
# Ks_su: Monod constant associated with the utilization of sugars, mol/L   
# km_aa: Maximum specific utilization rate constant of amino acids, mol/(mol h)
# Ks_aa: Monod constant associated with the utilization of amino acids mol/L 
# km_h2: Maximum specific utilization rate constant of hydrogen, mol/(mol h)
# Ks_h2: Monod constant associated with the utilization of hydrogen, mol/L
# k_br: Kinetic rate constant of bromoform utilization, 1/h
# p1: parameter of the sigmoid function associated with the bromoform inhibition factor (Ibr) 
# p2: parameter of the sigmoid function associated with the bromoform inhibition factor (Ibr)
# p3: parameter of the affine function associated with the flux allocation parameter towards acetate production (lambda1)
# p4: parameter of the affine function associated with the flux allocation parameter towards acetate production (lambda1)
# p5: parameter of the affine function associated with the flux allocation parameter towards propionate production (lambda2)   
# p6: parameter of the affine function associated with the flux allocation parameter towards propionate production (lambda2) 

# -- Range of variation of input parameters (uniform distributions set for exploring the variability of input parameters due to the lack of data available) --

# - Parameters khyd_ndf, khyd_nsc, khyd_pro, km_su, Ks_su, km_aa, Ks_aa, km_h2, Ks_h2 -

#Minimum and maximum values explored for each parameter
#khyd_ndf
Minkhyd_ndf<-0.01;Maxkhyd_ndf<-0.33
#khyd_nsc
Minkhyd_nsc<-0.06;Maxkhyd_nsc<-0.22
#khyd_pro
Minkhyd_pro<-0.05;Maxkhyd_pro<-0.25
#km_su
Minkm_su<-0.94;Maxkm_su<-4.33
#Ks_su
MinKs_su<-1e-04;MaxKs_su<-9e-03
#km_aa
Minkm_aa<-1;Maxkm_aa<-5
#Ks_aa
MinKs_aa<-3e-04;MaxKs_aa<-8e-03
#km_h2
Minkm_h2<-12;Maxkm_h2<-25
#Ks_h2
MinKs_h2<-1e-07;MaxKs_h2<-1e-05
  
# - Parameters of functions representing the bromoform inhibition factor (Ibr; p1 and p2),the flux allocation towards acetate production (lambda1; p3 and p4) and the flux allocation towards propionate production (lambda2; p5 and p6) -

#Percentage of variation considered around the value of Munoz-Tamayo et al., 2021 (similar to Dougherty et al., 2017)
p<-0.10 

#Range of variation explored for each parameter
#k_br
Initialvaluek_br<-0.095;Variabilityk_br<-c(Initialvaluek_br-(p*Initialvaluek_br),Initialvaluek_br+(p*Initialvaluek_br))
#p1
Initialvaluep1<-7.2551e+04;Variabilityp1<-c(Initialvaluep1-(p*Initialvaluep1),Initialvaluep1+(p*Initialvaluep1))
#p2
Initialvaluep2<--1.0837e-04;Variabilityp2<-c(Initialvaluep2+(p*Initialvaluep2),Initialvaluep2-(p*Initialvaluep2))
#p3
Initialvaluep3<-0.3655;Variabilityp3<-c(Initialvaluep3-(p*Initialvaluep3),Initialvaluep3+(p*Initialvaluep3))
#p4
Initialvaluep4<-0.6371;Variabilityp4<-c(Initialvaluep4-(p*Initialvaluep4),Initialvaluep4+(p*Initialvaluep4))
#p5
Initialvaluep5<-0.3787;Variabilityp5<-c(Initialvaluep5-(p*Initialvaluep5),Initialvaluep5+(p*Initialvaluep5))
#p6
Initialvaluep6<-0.1160;Variabilityp6<-c(Initialvaluep6-(p*Initialvaluep6),Initialvaluep6+(p*Initialvaluep6))

#Distribution fitted on input parameters
VarInputparameters<-data.frame(matrix(NA,nrow=k,ncol=3))
row.names(VarInputparameters)<-Inputs
colnames(VarInputparameters)<-c("Distribution","First_est_param","Second_est_param")

VarInputparameters[,"Distribution"]<-c("runif","runif","runif","runif","runif","runif",
                                       "runif","runif","runif","runif","runif","runif",
                                       "runif","runif","runif","runif")
VarInputparameters[,"First_est_param"]<-c(Minkhyd_ndf,Minkhyd_nsc,Minkhyd_pro,
                                          Minkm_su,MinKs_su,
                                          Minkm_aa,MinKs_aa,
                                          Minkm_h2,MinKs_h2,
                                          Variabilityk_br[1],
                                          Variabilityp1[1],Variabilityp2[1],
                                          Variabilityp3[1],Variabilityp4[1],
                                          Variabilityp5[1],Variabilityp6[1])
VarInputparameters[,"Second_est_param"]<-c(Maxkhyd_ndf,Maxkhyd_nsc,Maxkhyd_pro,
                                           Maxkm_su,MaxKs_su,
                                           Maxkm_aa,MaxKs_aa,
                                           Maxkm_h2,MaxKs_h2,
                                           Variabilityk_br[2],
                                           Variabilityp1[2],Variabilityp2[2],Variabilityp3[2],
                                           Variabilityp4[2],Variabilityp5[2],Variabilityp6[2])

# -- Functions used for sampling the input parameters with the random permutation method --

# - Function for sampling the input parameters -

Xall<-function(n){ 
  
  cbind(
    #Parameters of hydrolysis rate and substrate kinetic rate functions 
    runif(n,min=Minkhyd_ndf,max=Maxkhyd_ndf),
    runif(n,min=Minkhyd_nsc,max=Maxkhyd_nsc),
    runif(n,min=Minkhyd_pro,max=Maxkhyd_pro),
    runif(n,min=Minkm_su,max=Maxkm_su),
    runif(n,min=MinKs_su,max=MaxKs_su),
    runif(n,min=Minkm_aa,max=Maxkm_aa),
    runif(n,min=MinKs_aa,max=MaxKs_aa),
    runif(n,min=Minkm_h2,max=Maxkm_h2),
    runif(n,min=MinKs_h2,max=MaxKs_h2),
    #Constants and parameters of functions representing the effect of bromoform on the rumen fermentation
    runif(n,Variabilityk_br[1],Variabilityk_br[2]),
    runif(n,Variabilityp1[1],Variabilityp1[2]),
    runif(n,Variabilityp2[1],Variabilityp2[2]),
    runif(n,Variabilityp3[1],Variabilityp3[2]),
    runif(n,Variabilityp4[1],Variabilityp4[2]),
    runif(n,Variabilityp5[1],Variabilityp5[2]),
    runif(n,Variabilityp6[1],Variabilityp6[2]))
}

# - Function for sampling the input parameters considering the permutations -

Xset<-function(n,Sj,Sjc,xjc){
  
  if(length(Sj)>1){
    
    VarSjInputs<-matrix(NA,ncol=length(Sj),nrow=n)
    for(i in 1:length(Sj)){
      
      DrawDistributionfitted<-paste(VarInputparameters[match(Sj[i],row.names(VarInputparameters)),"Distribution"],
                                    "(",n,",",VarInputparameters[match(Sj[i],row.names(VarInputparameters)),"First_est_param"],
                                    ",",VarInputparameters[match(Sj[i],row.names(VarInputparameters)),"Second_est_param"],")",sep="")
      VarSjInputs[,i]<-eval(parse(text=DrawDistributionfitted))
      
      rm(DrawDistributionfitted)
    }
    
  } else{
    
    DrawDistributionfitted<-paste(VarInputparameters[match(Sj,row.names(VarInputparameters)),"Distribution"],
                                  "(",n,",",VarInputparameters[match(Sj,row.names(VarInputparameters)),"First_est_param"],
                                  ",",VarInputparameters[match(Sj,row.names(VarInputparameters)),"Second_est_param"],")",sep="")
    VarSjInputs<-matrix(eval(parse(text=DrawDistributionfitted)),nc=length(Sj))
    
    rm(DrawDistributionfitted)
    
  }
  
  return(VarSjInputs)
  
}

# -- Definition of the sampling matrix of input parameters --

# - Sampling of input parameters considering the random permutations - 

#Generation of the sampling matrix
SamplingInputparameters<-shapleyPermRand(model=NULL,Xall=Xall,Xset=Xset,d=k,Nv=1e4,m=1e4,No=1,Ni=3,colnames=Inputs)
SamplingMatrix<-SamplingInputparameters$SamplingInputparameters
write.table(data.frame(SamplingMatrix),
            file=paste(getwd(),"/InputparametersSamplingMatrix_txt/InputSamplingMatrix.txt",sep=""),
            row.names=FALSE,col.names=TRUE)

#Splitting of the sampling matrix in 20 parts
for(i in 1:20){
  
  SampMatSplitting<-seq(0,nrow(SamplingMatrix),by=nrow(SamplingMatrix)/20)
  eval(parse(text=paste("SamplingMatrixP",i,"<-SamplingMatrix[(seq(0,nrow(SamplingMatrix),by=nrow(SamplingMatrix)/20)[",i,"]+1):seq(0,nrow(SamplingMatrix),by=nrow(SamplingMatrix)/20)[",i+1,"],]",sep="")))
  write.table(data.frame(eval(parse(text=paste("SamplingMatrixP",i,sep="")))),
              file=paste(getwd(),"/InputparametersSamplingMatrix_txt/InputSamplingMatrixP",i,".txt",sep=""),
              row.names=FALSE,col.names=TRUE)
  
  rm(SampMatSplitting)
  
}

save.image(paste(getwd(),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))