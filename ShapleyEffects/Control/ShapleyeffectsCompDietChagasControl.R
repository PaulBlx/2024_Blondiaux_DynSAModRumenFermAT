# -----------------------------------------------------------------------------------------------------------

# Computation of the Shapley effects using the random permutations method (Song et al., 2016) for the control

# -----------------------------------------------------------------------------------------------------------

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
MethaneDietChagasControl<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/MethaneDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Shapley effects computation
ShapleyResultsMethaneDietChagasControl<-apply(MethaneDietChagasControl,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsMethaneDietChagasControl)<-colnames(MethaneDietChagasControl)

#Estimation
ListShapleyEffectsMethaneDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasControl,function(x){x$Shapley$original})
ShapleyEffectsMethaneDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsMethaneDietChagasControl))
row.names(ShapleyEffectsMethaneDietChagasControl)<-Inputs
colnames(ShapleyEffectsMethaneDietChagasControl)<-names(ShapleyResultsMethaneDietChagasControl)

#Standard error
ListShapleyEffectsStdErrMethaneDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasControl,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrMethaneDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrMethaneDietChagasControl))
row.names(ShapleyEffectsStdErrMethaneDietChagasControl)<-Inputs
colnames(ShapleyEffectsStdErrMethaneDietChagasControl)<-names(ShapleyResultsMethaneDietChagasControl)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsMethaneDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrMethaneDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsStdErrMethaneDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L) among the output variables
AcetateDietChagasControl<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/AcetateDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Shapley effects computation
ShapleyResultsAcetateDietChagasControl<-apply(AcetateDietChagasControl,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsAcetateDietChagasControl)<-colnames(AcetateDietChagasControl)

#Estimation
ListShapleyEffectsAcetateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasControl,function(x){x$Shapley$original})
ShapleyEffectsAcetateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsAcetateDietChagasControl))
row.names(ShapleyEffectsAcetateDietChagasControl)<-Inputs
colnames(ShapleyEffectsAcetateDietChagasControl)<-names(ShapleyResultsAcetateDietChagasControl)

#Standard error
ListShapleyEffectsStdErrAcetateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasControl,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrAcetateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrAcetateDietChagasControl))
row.names(ShapleyEffectsStdErrAcetateDietChagasControl)<-Inputs
colnames(ShapleyEffectsStdErrAcetateDietChagasControl)<-names(ShapleyResultsAcetateDietChagasControl)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsAcetateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsAcetateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrAcetateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsStdErrAcetateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L) among the output variables
ButyrateDietChagasControl<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/ButyrateDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Shapley effects computation
ShapleyResultsButyrateDietChagasControl<-apply(ButyrateDietChagasControl,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsButyrateDietChagasControl)<-colnames(ButyrateDietChagasControl)

#Estimation
ListShapleyEffectsButyrateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasControl,function(x){x$Shapley$original})
ShapleyEffectsButyrateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsButyrateDietChagasControl))
row.names(ShapleyEffectsButyrateDietChagasControl)<-Inputs
colnames(ShapleyEffectsButyrateDietChagasControl)<-names(ShapleyResultsButyrateDietChagasControl)

#Standard error
ListShapleyEffectsStdErrButyrateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasControl,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrButyrateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrButyrateDietChagasControl))
row.names(ShapleyEffectsStdErrButyrateDietChagasControl)<-Inputs
colnames(ShapleyEffectsStdErrButyrateDietChagasControl)<-names(ShapleyResultsButyrateDietChagasControl)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsButyrateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrButyrateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsStdErrButyrateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L) among the output variables
PropionateDietChagasControl<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/PropionateDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Shapley effects computation
ShapleyResultsPropionateDietChagasControl<-apply(PropionateDietChagasControl,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsPropionateDietChagasControl)<-colnames(PropionateDietChagasControl)

#Estimation
ListShapleyEffectsPropionateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasControl,function(x){x$Shapley$original})
ShapleyEffectsPropionateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsPropionateDietChagasControl))
row.names(ShapleyEffectsPropionateDietChagasControl)<-Inputs
colnames(ShapleyEffectsPropionateDietChagasControl)<-names(ShapleyResultsPropionateDietChagasControl)

#Standard error
ListShapleyEffectsStdErrPropionateDietChagasControl<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasControl,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrPropionateDietChagasControl<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrPropionateDietChagasControl))
row.names(ShapleyEffectsStdErrPropionateDietChagasControl)<-Inputs
colnames(ShapleyEffectsStdErrPropionateDietChagasControl)<-names(ShapleyResultsPropionateDietChagasControl)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsPropionateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrPropionateDietChagasControl,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsStdErrPropionateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl) #Stop the clusters
