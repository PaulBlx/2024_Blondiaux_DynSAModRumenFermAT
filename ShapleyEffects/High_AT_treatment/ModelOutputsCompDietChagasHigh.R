# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the output variables of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied in the sensitivity analysis analysis using the sampling matrix of input parameters simulated for the high Asparagopsis taxiformis treatment 

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
XDietChagasHigh<-Cinit_Control #Initial conditions of the system
auxsDietChagasHigh<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrHigh) #Diet of Chagas et al., 2019 and treatment high

#Computation of output variables of the model
XmDietChagasHigh<-parApply(cl=cl,X=SamplingMatrix,MARGIN=1,FUN=ModelRun,Cinit=XDietChagasHigh,t=t,auxs=auxsDietChagasHigh)
#Selection of the 4th day of simulation
XmDietChagasHigh<-parLapply(cl=cl,X=XmDietChagasHigh,fun=function(x){return(x[RowSimd4,])})
#Selection of the output variables of interest (methane and volatile fatty acids)
SAOutVarConcDietChagasHigh<-parLapply(cl=cl,X=XmDietChagasHigh,fun=rumencATout,V_l=V_l,V_g=V_g)

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h) among the output variables
ListMethaneDietChagasHigh<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(x){x[,"q_ch4g_out"]})
MethaneDietChagasHigh<-as.data.frame(do.call(rbind,ListMethaneDietChagasHigh))
row.names(MethaneDietChagasHigh)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(MethaneDietChagasHigh)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of methane concentration dynamic simulated in the output files
write.table(MethaneDietChagasHigh[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/MethaneDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L) among the output variables
ListAcetateDietChagasHigh<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(x){x[,"s_ac"]})
AcetateDietChagasHigh<-as.data.frame(do.call(rbind,ListAcetateDietChagasHigh))
row.names(AcetateDietChagasHigh)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(AcetateDietChagasHigh)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of acetate concentration dynamic simulated in the output files
write.table(AcetateDietChagasHigh[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/AcetateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L) among the output variables
ListButyrateDietChagasHigh<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(x){x[,"s_bu"]})
ButyrateDietChagasHigh<-as.data.frame(do.call(rbind,ListButyrateDietChagasHigh))
row.names(ButyrateDietChagasHigh)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(ButyrateDietChagasHigh)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of butyrate concentration dynamic simulated in the output files
write.table(ButyrateDietChagasHigh[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/ButyrateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L) among the output variables
ListPropionateDietChagasHigh<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(x){x[,"s_pr"]})
PropionateDietChagasHigh<-as.data.frame(do.call(rbind,ListPropionateDietChagasHigh))
row.names(PropionateDietChagasHigh)<-paste("Simulation",seq(1,nrow(SamplingMatrix),by=1),sep="")
colnames(PropionateDietChagasHigh)<-paste("t=",round(t[RowSimd4],3),"h",sep="")

#Integration of propionate concentration dynamic simulated in the output files
write.table(PropionateDietChagasHigh[,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
            file=paste(getwd(),"/ModRumFermATOutVar_txt/PropionateDietChagasHigh.txt",sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl) #Stop the clusters
