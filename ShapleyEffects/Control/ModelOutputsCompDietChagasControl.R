# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the output variables of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied in the sensitivity analysis analysis using the sampling matrix of input parameters simulated for the control

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# ---- Computation of model outputs (Running the model using parallel computing) -----

# -- Function for running the model (solving the ode system) -- 
ModelRun<-function(x,Cinit,t,auxs){
  
  #################################################################################
  # List of inputs to this function:
  # x: Sampling matrix of inputs
  # Cinit: Initial conditions of the state variables
  # t: vector of time
  # auxs: auxiliary parameters used in the model
  #################################################################################
  
  Inputparams<<-x
  Xm<-data.frame(ode(y=Cinit,times=t,func=rumencAT,
                     parms=auxs,method="lsode"))
  
  return(Xm)
  
}

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

# ---- Running the model ----

#Parameters of the simulations
XDietChagasControl<-Cinit_Control #Initial conditions of the system
auxsDietChagasControl<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrControl) #Diet of Chagas et al., 2019 and treatment control

#Computation of output variables of the model
XmDietChagasControl<-parApply(cl=cl,X=SamplingMatrix,MARGIN=1,FUN=ModelRun,Cinit=XDietChagasControl,t=t,auxs=auxsDietChagasControl)
#Selection of the 4th day of simulation
XmDietChagasControl<-parLapply(cl=cl,X=XmDietChagasControl,fun=function(x){return(x[RowSimd4,])})
#Selection of the output variables of interest (methane and volatile fatty acids)
SAOutVarConcDietChagasControl<-parLapply(cl=cl,X=XmDietChagasControl,fun=rumencATout,V_l=V_l,V_g=V_g)

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h) among the output variables
ListMethaneDietChagasControl<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(x){x[,"q_ch4g_out"]})
MethaneDietChagasControl<-as.data.frame(do.call(rbind,ListMethaneDietChagasControl))
row.names(MethaneDietChagasControl)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(MethaneDietChagasControl)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of methane concentration dynamic simulated in the output files
write.table(MethaneDietChagasControl[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/MethaneDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L) among the output variables
ListAcetateDietChagasControl<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(x){x[,"s_ac"]})
AcetateDietChagasControl<-as.data.frame(do.call(rbind,ListAcetateDietChagasControl))
row.names(AcetateDietChagasControl)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(AcetateDietChagasControl)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of acetate concentration dynamic simulated in the output files
write.table(AcetateDietChagasControl[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/AcetateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L) among the output variables
ListButyrateDietChagasControl<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(x){x[,"s_bu"]})
ButyrateDietChagasControl<-as.data.frame(do.call(rbind,ListButyrateDietChagasControl))
row.names(ButyrateDietChagasControl)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(ButyrateDietChagasControl)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of butyrate concentration dynamic simulated in the output files
write.table(ButyrateDietChagasControl[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/ButyrateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L) among the output variables
ListPropionateDietChagasControl<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(x){x[,"s_pr"]})
PropionateDietChagasControl<-as.data.frame(do.call(rbind,ListPropionateDietChagasControl))
row.names(PropionateDietChagasControl)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(PropionateDietChagasControl)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of propionate concentration dynamic simulated in the output files
write.table(PropionateDietChagasControl[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/PropionateDietChagasControl.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl) #Stop the clusters
