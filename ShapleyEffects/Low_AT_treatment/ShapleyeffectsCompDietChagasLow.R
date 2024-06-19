# -----------------------------------------------------------------------------------------------------------------------------------------

# Computation of the Shapley effects using the random permutations method (Song et al., 2016) for the low Asparagopsis taxiformis treatment

# -----------------------------------------------------------------------------------------------------------------------------------------

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
MethaneDietChagasLow<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/MethaneDietChagasLow.txt",sep=""),
                                 header=TRUE)

#Shapley effects computation
ShapleyResultsMethaneDietChagasLow<-apply(MethaneDietChagasLow,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsMethaneDietChagasLow)<-colnames(MethaneDietChagasLow)

#Estimation
ListShapleyEffectsMethaneDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasLow,function(x){x$Shapley$original})
ShapleyEffectsMethaneDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsMethaneDietChagasLow))
row.names(ShapleyEffectsMethaneDietChagasLow)<-Inputs
colnames(ShapleyEffectsMethaneDietChagasLow)<-names(ShapleyResultsMethaneDietChagasLow)

#Standard error
ListShapleyEffectsStdErrMethaneDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsMethaneDietChagasLow,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrMethaneDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrMethaneDietChagasLow))
row.names(ShapleyEffectsStdErrMethaneDietChagasLow)<-Inputs
colnames(ShapleyEffectsStdErrMethaneDietChagasLow)<-names(ShapleyResultsMethaneDietChagasLow)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsMethaneDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrMethaneDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Methane/ShapleyEffectsStdErrMethaneDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L) among the output variables
AcetateDietChagasLow<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/AcetateDietChagasLow.txt",sep=""),
                                 header=TRUE)

#Shapley effects computation
ShapleyResultsAcetateDietChagasLow<-apply(AcetateDietChagasLow,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsAcetateDietChagasLow)<-colnames(AcetateDietChagasLow)

#Estimation
ListShapleyEffectsAcetateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasLow,function(x){x$Shapley$original})
ShapleyEffectsAcetateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsAcetateDietChagasLow))
row.names(ShapleyEffectsAcetateDietChagasLow)<-Inputs
colnames(ShapleyEffectsAcetateDietChagasLow)<-names(ShapleyResultsAcetateDietChagasLow)

#Standard error
ListShapleyEffectsStdErrAcetateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsAcetateDietChagasLow,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrAcetateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrAcetateDietChagasLow))
row.names(ShapleyEffectsStdErrAcetateDietChagasLow)<-Inputs
colnames(ShapleyEffectsStdErrAcetateDietChagasLow)<-names(ShapleyResultsAcetateDietChagasLow)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsAcetateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsAcetateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrAcetateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Acetate/ShapleyEffectsStdErrAcetateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L) among the output variables
ButyrateDietChagasLow<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/ButyrateDietChagasLow.txt",sep=""),
                                  header=TRUE)

#Shapley effects computation
ShapleyResultsButyrateDietChagasLow<-apply(ButyrateDietChagasLow,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsButyrateDietChagasLow)<-colnames(ButyrateDietChagasLow)

#Estimation
ListShapleyEffectsButyrateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasLow,function(x){x$Shapley$original})
ShapleyEffectsButyrateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsButyrateDietChagasLow))
row.names(ShapleyEffectsButyrateDietChagasLow)<-Inputs
colnames(ShapleyEffectsButyrateDietChagasLow)<-names(ShapleyResultsButyrateDietChagasLow)

#Standard error
ListShapleyEffectsStdErrButyrateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsButyrateDietChagasLow,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrButyrateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrButyrateDietChagasLow))
row.names(ShapleyEffectsStdErrButyrateDietChagasLow)<-Inputs
colnames(ShapleyEffectsStdErrButyrateDietChagasLow)<-names(ShapleyResultsButyrateDietChagasLow)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsButyrateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrButyrateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Butyrate/ShapleyEffectsStdErrButyrateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L) among the output variables
PropionateDietChagasLow<-read.table(paste(getwd(),"/ModRumFermATOutVar_txt/PropionateDietChagasLow.txt",sep=""),
                                    header=TRUE)

#Shapley effects computation
ShapleyResultsPropionateDietChagasLow<-apply(PropionateDietChagasLow,2,FUN=ShapleyEffectsComputation)
names(ShapleyResultsPropionateDietChagasLow)<-colnames(PropionateDietChagasLow)

#Estimation
ListShapleyEffectsPropionateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasLow,function(x){x$Shapley$original})
ShapleyEffectsPropionateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsPropionateDietChagasLow))
row.names(ShapleyEffectsPropionateDietChagasLow)<-Inputs
colnames(ShapleyEffectsPropionateDietChagasLow)<-names(ShapleyResultsPropionateDietChagasLow)

#Standard error
ListShapleyEffectsStdErrPropionateDietChagasLow<-parLapply(cl=cl,X=ShapleyResultsPropionateDietChagasLow,function(x){x$Shapley$`std. error`})
ShapleyEffectsStdErrPropionateDietChagasLow<-as.data.frame(do.call(cbind,ListShapleyEffectsStdErrPropionateDietChagasLow))
row.names(ShapleyEffectsStdErrPropionateDietChagasLow)<-Inputs
colnames(ShapleyEffectsStdErrPropionateDietChagasLow)<-names(ShapleyResultsPropionateDietChagasLow)

#Integration of Shapley effects (estimation + standard error) estimated in the output files
#Estimation
write.table(ShapleyEffectsPropionateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

#Standard error
write.table(ShapleyEffectsStdErrPropionateDietChagasLow,
            file=paste(getwd(),"/ShapEff_txt/Propionate/ShapleyEffectsStdErrPropionateDietChagasLow.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl) #Stop the clusters
