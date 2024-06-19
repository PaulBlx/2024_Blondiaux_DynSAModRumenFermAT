# ------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the Shapley effects using the random permutations method (Song et al., 2016) for the high Asparagopsis taxiformis treatment

# ------------------------------------------------------------------------------------------------------------------------------------------

library(deSolve)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(pracma)
library(sensitivity)
library(parallel)

#Recovering of external argument (RData containing the sampling matrix of input parameters and functions for running the mechanistic model and computing the Shapley effects)
load(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))

# -- General parameters of the simulations and of the parallel computing --

#Time scale considered
nd<-4 #number of days simulated
ts<-1/60 #time step, hour
t<-seq(0,24*nd,by=ts) #time in hour
RowSimd4<-seq(length(seq(0,24*(nd-1),by=ts))+1,length(seq(0,24*nd,by=ts)),by=1) #Lines corresponding to the 4th day of simulation

#Rumen volume in Rusitec condition
V_l<-0.74 #Volume in liquid phase of the rumen, L (Belanche et al., 2017)
V_g<-0.06 #Volume in gas phase of the rumen, L (Belanche et al., 2017)

#Parallel computing
#Settings
nbnodes<-detectCores()
cl<-makeCluster(nbnodes)
#Library used for solving the model (for sending the information to cores)
clusterEvalQ(cl,c(library(deSolve),library(gridExtra),library(scales),library(pracma),library(sensitivity)))
#Functions used for solving the model (for sending the information to cores)
clusterExport(cl,c("ode","rumencAT","tDMI","yDMI","tell","SamplingInputparameters","RowSimd4"),envir=environment())

# ---- Estimation of Shapley effects by randomly sampling the permutations of inputs (Method of Song et al., 2016)  ----

#Estimation of Shapley effects using the method of Song et al., 2016 (Using formula of Castro et al., 2009)
ShapleyEffectsComputation<-function(x){
  
  #################################################################################
  # List of inputs to this function:
  # x: vector of the output variable
  #################################################################################
  
  ShapleyResultsTimestep<-tell(x=SamplingInputparameters$sh,y=x)
  
  return(ShapleyResultsTimestep)
  
} 

# ---- Computing the Shapley effects ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h) among the output variables
MethaneDietChagasHigh<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/MethaneDietChagasHigh.txt",sep=""),
                                  header=TRUE)

#Shapley effects computation
ShapleyResultsMethaneDietChagasHigh<-apply(MethaneDietChagasHigh,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsMethaneDietChagasHigh)<-colnames(MethaneDietChagasHigh)

#Estimation
ListShapleyEffectsMethaneDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasHigh,function(x){x$Shapley$original})
ShapleyEffectsMethaneDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsMethaneDietChagasHigh))
row.names(ShapleyEffectsMethaneDietChagasHigh)<-Inputs
colnames(ShapleyEffectsMethaneDietChagasHigh)<-names(ShapleyResultsMethaneDietChagasHigh)

#Standard error
ListShapleyEffectsStdErrMethaneDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasHigh,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrMethaneDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrMethaneDietChagasHigh))
row.names(ShapleyEffectsStdErrMethaneDietChagasHigh)<-Inputs
colnames(ShapleyEffectsStdErrMethaneDietChagasHigh)<-names(ShapleyResultsMethaneDietChagasHigh)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsMethaneDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrMethaneDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsStdErrMethaneDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L) among the output variables
AcetateDietChagasHigh<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/AcetateDietChagasHigh.txt",sep=""),
                                  header=TRUE)

#Shapley effects computation
ShapleyResultsAcetateDietChagasHigh<-apply(AcetateDietChagasHigh,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsAcetateDietChagasHigh)<-colnames(AcetateDietChagasHigh)

#Estimation
ListShapleyEffectsAcetateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasHigh,function(x){x$Shapley$original})
ShapleyEffectsAcetateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsAcetateDietChagasHigh))
row.names(ShapleyEffectsAcetateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsAcetateDietChagasHigh)<-names(ShapleyResultsAcetateDietChagasHigh)

#Standard error
ListShapleyEffectsStdErrAcetateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasHigh,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrAcetateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrAcetateDietChagasHigh))
row.names(ShapleyEffectsStdErrAcetateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsStdErrAcetateDietChagasHigh)<-names(ShapleyResultsAcetateDietChagasHigh)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsAcetateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsAcetateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrAcetateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsStdErrAcetateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L) among the output variables
ButyrateDietChagasHigh<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/ButyrateDietChagasHigh.txt",sep=""),
                                   header=TRUE)

#Shapley effects computation
ShapleyResultsButyrateDietChagasHigh<-apply(ButyrateDietChagasHigh,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHigh)

#Estimation
ListShapleyEffectsButyrateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasHigh,function(x){x$Shapley$original})
ShapleyEffectsButyrateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsButyrateDietChagasHigh))
row.names(ShapleyEffectsButyrateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsButyrateDietChagasHigh)<-names(ShapleyResultsButyrateDietChagasHigh)

#Standard error
ListShapleyEffectsStdErrButyrateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasHigh,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrButyrateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrButyrateDietChagasHigh))
row.names(ShapleyEffectsStdErrButyrateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsStdErrButyrateDietChagasHigh)<-names(ShapleyResultsButyrateDietChagasHigh)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsButyrateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrButyrateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsStdErrButyrateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L) among the output variables
PropionateDietChagasHigh<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/PropionateDietChagasHigh.txt",sep=""),
                                     header=TRUE)

#Shapley effects computation
ShapleyResultsPropionateDietChagasHigh<-apply(PropionateDietChagasHigh,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsPropionateDietChagasHigh)<-colnames(PropionateDietChagasHigh)

#Estimation
ListShapleyEffectsPropionateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasHigh,function(x){x$Shapley$original})
ShapleyEffectsPropionateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsPropionateDietChagasHigh))
row.names(ShapleyEffectsPropionateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsPropionateDietChagasHigh)<-names(ShapleyResultsPropionateDietChagasHigh)

#Standard error
ListShapleyEffectsStdErrPropionateDietChagasHigh<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasHigh,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrPropionateDietChagasHigh<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrPropionateDietChagasHigh))
row.names(ShapleyEffectsStdErrPropionateDietChagasHigh)<-Inputs
colnames(ShapleyEffectsStdErrPropionateDietChagasHigh)<-names(ShapleyResultsPropionateDietChagasHigh)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsPropionateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrPropionateDietChagasHigh,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsStdErrPropionateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl) #Stop the clusters
