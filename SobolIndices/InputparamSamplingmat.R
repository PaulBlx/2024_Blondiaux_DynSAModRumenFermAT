# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Simulation of the sampling matrix of input parameters of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(sensobol)
library(data.table)
library(ggplot2)
library(ggpubr)
library(boot)
library(rvinecopulib)
library(copula)
library(VineCopula)
library(combinat)

#Functions for running the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021)
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencATload.R",sep=""))
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencATout.R",sep=""))
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/rumencAT.R",sep=""))

#Functions for computing Sobol indices
source(paste(getwd(),"/BootstrapimplementationFunctions.R",sep=""))
source(paste(getwd(),"/SobolindicesCompFunctions.R",sep=""))
source(paste(getwd(),"/PickandFreezeStratFunctions.R",sep=""))

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

VarInputparameters[,"Distribution"]<-c("unif","unif","unif","unif","unif","unif",
                                       "unif","unif","unif","unif","unif","unif",
                                       "unif","unif","unif","unif")
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

# -- Definition of the sampling matrix of input parameters --

# - Sampling of input parameters using the Sobol sequences and considering 3000 simulations - 

#Number of simulations
Nsim<-3e03

#Circular reordering of input parameters
InputsRT1<-c("khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2",
             "k_br","p1","p2","p3","p4","p5","p6") #full=khyd_ndf; independent=p6

InputsRT2<-c("khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1",
             "p2","p3","p4","p5","p6","khyd_ndf") #full=khyd_nsc; independent=khyd_ndf

InputsRT3<-c("khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1","p2",
             "p3","p4","p5","p6","khyd_ndf","khyd_nsc") #full=khyd_pro; independent=khyd_nsc

InputsRT4<-c("km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1","p2","p3",
             "p4","p5","p6","khyd_ndf","khyd_nsc","khyd_pro") #full=km_su; independent=khyd_pro

InputsRT5<-c("Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1","p2","p3","p4",
             "p5","p6","khyd_ndf","khyd_nsc","khyd_pro","km_su") #full=Ks_su; independent=km_su

InputsRT6<-c("km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1","p2","p3","p4","p5",
             "p6","khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su") #full=km_aa; independent=Ks_su

InputsRT7<-c("Ks_aa","km_h2","Ks_h2","k_br","p1","p2","p3","p4","p5","p6",
             "khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa") #full=Ks_aa; independent=km_aa

InputsRT8<-c("km_h2","Ks_h2","k_br","p1","p2","p3","p4","p5","p6","khyd_ndf",
             "khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa") #full=km_h2; independent=Ks_aa

InputsRT9<-c("Ks_h2","k_br","p1","p2","p3","p4","p5","p6","khyd_ndf","khyd_nsc",
             "khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2") #full=Ks_h2; independent=km_h2

InputsRT10<-c("k_br","p1","p2","p3","p4","p5","p6","khyd_ndf","khyd_nsc",
              "khyd_pro","km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2") #full=k_br; independent=Ks_h2

InputsRT11<-c("p1","p2","p3","p4","p5","p6","khyd_ndf","khyd_nsc","khyd_pro",
              "km_su","Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br") #full=p1; independent=k_br

InputsRT12<-c("p2","p3","p4","p5","p6","khyd_ndf","khyd_nsc","khyd_pro","km_su",
              "Ks_su","km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1") #full=p2; independent=p1

InputsRT13<-c("p3","p4","p5","p6","khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su",
              "km_aa","Ks_aa","km_h2","Ks_h2","k_br","p1","p2") #full=p3; independent=p2

InputsRT14<-c("p4","p5","p6","khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa",
              "Ks_aa","km_h2","Ks_h2","k_br","p1","p2","p3") #full=p4; independent=p3

InputsRT15<-c("p5","p6","khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa",
              "km_h2","Ks_h2","k_br","p1","p2","p3","p4") #full=p5; independent=p4

InputsRT16<-c("p6","khyd_ndf","khyd_nsc","khyd_pro","km_su","Ks_su","km_aa","Ks_aa",
              "km_h2","Ks_h2","k_br","p1","p2","p3","p4","p5") #full=p6; independent=p5

#Definition of the sampling matrices for each circular ordering
ListSamplesInputs<-vector("list",length=k)
for(i in 1:k){
  
  #Sampling with the Sobol sequences
  Samplinginputs<-sobol_matrices(matrices=c("A","B","BA"),N=Nsim,params=eval(parse(text=paste("InputsRT",i,sep=""))),
                                 order="first",type="QRN") 
  
  #Selection of the two independent sampling matrices
  SamplinginputsA<-Samplinginputs[1:Nsim,];SamplinginputsB<-Samplinginputs[(Nsim+1):(2*Nsim),]
  
  #Transformation of samples generated in [0,1] with sobol_matrices in the "real" range of variation (defined by VarInputparameters)
  for(j in 1:length(Inputs)){
    
    SamplinginputsA[,Inputs[j]]<-eval(parse(text=paste("q",VarInputparameters[Inputs[j],"Distribution"],
                                                       "(SamplinginputsA[,'",Inputs[j],"'],",VarInputparameters[Inputs[j],"First_est_param"],
                                                       ",",VarInputparameters[Inputs[j],"Second_est_param"],")",sep="")))
    
    SamplinginputsB[,Inputs[j]]<-eval(parse(text=paste("q",VarInputparameters[Inputs[j],"Distribution"],
                                                       "(SamplinginputsB[,'",Inputs[j],"'],",VarInputparameters[Inputs[j],"First_est_param"],
                                                       ",",VarInputparameters[Inputs[j],"Second_est_param"],")",sep="")))
    
  }
  
  ListSamplesInputs[[i]]<-list(A=SamplinginputsA,B=SamplinginputsB)
  rm(Samplinginputs);rm(SamplinginputsA);rm(SamplinginputsB)
  
}

# -- Implementation of the Rosenblatt transformation to have an independent correlation structure of the input parameters --

#Transformation of the sampling matrices using the Rosenblatt transformation
listRTSamplesInputs<-vector("list",length=k)
for(i in 1:k){
  
  #Selection of original sampling matrices
  A<-ListSamplesInputs[[i]]$A;B<-ListSamplesInputs[[i]]$B
  
  #Estimation of a vine distribution model
  FitCopulaSamplesInputsRTiA<-vine(A,copula_controls=list(family_set="par"))
  FitCopulaSamplesInputsRTiB<-vine(B,copula_controls=list(family_set="par"))
  
  #Rosenblatt transformation of the original sampling matrices
  RTSamplesInputsRTiA<-rosenblatt(A,FitCopulaSamplesInputsRTiA)
  RTSamplesInputsRTiB<-rosenblatt(B,FitCopulaSamplesInputsRTiB)
  
  #Transformation of samples generated in [0,1] in the "real" range of variation (defined by VarInputparameters)
  for(j in 1:length(Inputs)){
    
    RTSamplesInputsRTiA[,Inputs[j]]<-eval(parse(text=paste("q",VarInputparameters[Inputs[j],"Distribution"],
                                                           "(RTSamplesInputsRTiA[,'",Inputs[j],"'],",VarInputparameters[Inputs[j],"First_est_param"],
                                                           ",",VarInputparameters[Inputs[j],"Second_est_param"],")",sep="")))
    
    RTSamplesInputsRTiB[,Inputs[j]]<-eval(parse(text=paste("q",VarInputparameters[Inputs[j],"Distribution"],
                                                           "(RTSamplesInputsRTiB[,'",Inputs[j],"'],",VarInputparameters[Inputs[j],"First_est_param"],
                                                           ",",VarInputparameters[Inputs[j],"Second_est_param"],")",sep="")))
    
  }
  
  listRTSamplesInputs[[i]]<-list(A=RTSamplesInputsRTiA,B=RTSamplesInputsRTiB)
  
  rm(A);rm(B);rm(FitCopulaSamplesInputsRTiA);rm(FitCopulaSamplesInputsRTiB);
  rm(RTSamplesInputsRTiA);rm(RTSamplesInputsRTiB)
  
}

# -- Implementation of the "Pick and freeze" strategy to generate the matrix for computing Sobol indices --

#Generation of the sampling matrix of shape (A,B,BA) used for computing Sobol indices
listRTSamplesInputsPickandfreeze<-vector("list",length=k)
for(i in 1:length(listRTSamplesInputs)){
  
  listRTSamplesInputsPickandfreeze[[i]]<-PickandFreeze(A=listRTSamplesInputs[[i]]$A,
                                                       B=listRTSamplesInputs[[i]]$B,
                                                       matrices=c("A","B","BA"),N=Nsim,
                                                       params=eval(parse(text=paste("InputsRT",i,sep=""))),
                                                       order="first")
  
  write.table(listRTSamplesInputsPickandfreeze[[i]],
              paste(getwd(),"/InputparametersSamplingMatrix_txt/RT",i,"/InputSamplingMatrixRT",i,".txt",sep=""),
              sep=" ",row.names=F,col.names=T)
  
}

#Splitting of the sampling matrix in 20 parts
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  SampMatSplitting<-seq(0,nrow(listRTSamplesInputsPickandfreeze[[i]]),by=nrow(listRTSamplesInputsPickandfreeze[[i]])/20)
  
  for(j in 1:20){
    
    eval(parse(text=paste("RT",i,"SamplingMatrixP",j,"<-listRTSamplesInputsPickandfreeze[[",i,"]][(SampMatSplitting[",j,"]+1):SampMatSplitting[",j+1,"],]",sep="")))
    write.table(data.frame(eval(parse(text=paste("RT",i,"SamplingMatrixP",j,sep="")))),
                file=paste(getwd(),"/InputparametersSamplingMatrix_txt/RT",i,"/RT",i,"InputSamplingMatrixP",j,".txt",sep=""),
                row.names=FALSE,col.names=TRUE)
    
  }
  
  rm(SampMatSplitting)
  
}

save.image(paste(getwd(),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))  