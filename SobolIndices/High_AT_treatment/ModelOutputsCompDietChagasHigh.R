# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the output variables of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied in the sensitivity analysis analysis using the sampling matrix of input parameters simulated for the high Asparagopsis taxiformis treatment

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(data.table)
library(boot)
library(parallel)
library(deSolve)

#Loading RData containing input parameter sampling matrix and functions for running mechanistic model and computing full and independent Sobol indices
load(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))

# ---- Computation of model outputs (Running the model using parallel computing) -----

# -- Function for computing model output (solving ode system) -- 
ModelRun<-function(x,Cinit,t,auxs){
  
  #################################################################################
  # List of inputs to this function:
  # x: Sampling matrix of inputs
  # Cinit: Initial conditions of state variables
  # t: vector of time
  # auxs:parameters used in the model
  #################################################################################
  
  Inputparams<<-x
  Xm<-data.frame(ode(y=Cinit,times=t,func=rumencAT,
                     parms=auxs,method="lsode"))
  
  return(Xm)
  
}

# -- General parameters of simulations and parallel computing implementation --

#Time scale considered
nd<-4 #number of days simulated
ts<-1/60 #time step, hour
t<-seq(0,24*nd,by=ts) #time in hour
RowSimd4<-seq(length(seq(0,24*(nd-1),by=ts))+1,length(seq(0,24*nd,by=ts)),by=1) #Lines corresponding to the 4th day of simulation

#Rumen volume in Rusitec condition
V_l<-0.74 #Volume in liquid phase of the rumen, L (Communication Peter Moate)
V_g<-0.06 #Volume in gas phase of the rumen, L (Communication Peter Moate)

#Parallel computing
#Settings
nbnodes<-detectCores()
cl<-makeCluster(nbnodes)
#Library used for solving the model (for sending information to clusters)
clusterEvalQ(cl,c(library(deSolve),library(gridExtra),library(pracma),library(scales),library(boot)))
#Functions used for solving the model (for sending information to clusters)
clusterExport(cl,c("ode","rumencAT","tDMI","yDMI","listRTSamplesInputsPickandfreeze","sobol_indices","RowSimd4"),envir=environment())

# ---- Running the model ----

#Computation of output variables of the model
XDietChagasHigh<-Cinit_Control #Initial conditions of the system
auxsDietChagasHigh<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrHigh) #High forage diet and treatment high
XmDietChagasHighRT<-vector("list",length=length(listRTSamplesInputsPickandfreeze))
for(i in 1:length(XmDietChagasHighRT)){
  
  eval(parse(text=paste("SamplingMatrix<-listRTSamplesInputsPickandfreeze[[",i,"]]",sep="")))
  XmDietChagasHighRT[[i]]<-parApply(cl=cl,X=SamplingMatrix,MARGIN=1,FUN=ModelRun,Cinit=XDietChagasHigh,t=t,auxs=auxsDietChagasHigh) #Solving ode system
  
  rm(SamplingMatrix)
  
}

#Selection of the 4th day of simulation
XmDietChagasHighRTD4<-vector("list",length=length(XmDietChagasHighRT))
for(i in 1:length(XmDietChagasHighRTD4)){
  
  XmDietChagasHigh<-XmDietChagasHighRT[[i]]
  XmDietChagasHighRTD4[[i]]<-parLapply(cl=cl,X=XmDietChagasHigh,fun=function(x){return(x[RowSimd4,])})
  
  rm(XmDietChagasHigh)
  
}
XmDietChagasHighRT<-XmDietChagasHighRTD4

#Selection of output variables of interest (methane and volatile fatty acids (VFA) concentrations)
SAOutVarConcDietChagasHighRT<-vector("list",length=length(XmDietChagasHighRT))
for(i in 1:length(SAOutVarConcDietChagasHighRT)){
  
  XmDietChagasHigh<-XmDietChagasHighRT[[i]]
  SAOutVarConcDietChagasHighRT[[i]]<-parLapply(cl=cl,X=XmDietChagasHigh,fun=rumencATout,V_l=V_l,V_g=V_g)
  
  rm(XmDietChagasHigh)
  
}

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
ListMethaneDietChagasHighRT<-vector("list",length=length(SAOutVarConcDietChagasHighRT))
for(i in 1:length(ListMethaneDietChagasHighRT)){
  
  SAOutVarConcDietChagasHigh<-SAOutVarConcDietChagasHighRT[[i]]
  ListMethaneDietChagasHighRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(y){y[,"q_ch4g_out"]})
  
  rm(SAOutVarConcDietChagasHigh)
  
}

MethaneDietChagasHighRT<-parLapply(cl=cl,X=ListMethaneDietChagasHighRT,
                                     function(x){
                                       ListMethaneDietChagasHigh<-x
                                       as.data.frame(do.call(rbind,ListMethaneDietChagasHigh))
                                       }
                                     )
for(i in 1:length(MethaneDietChagasHighRT)){
  
  row.names(MethaneDietChagasHighRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(MethaneDietChagasHighRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of methane concentration dynamic simulated in output files
for(i in 1:length(MethaneDietChagasHighRT)){
  
  write.table(MethaneDietChagasHighRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasHighRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
ListAcetateDietChagasHighRT<-vector("list",length=length(SAOutVarConcDietChagasHighRT))
for(i in 1:length(ListAcetateDietChagasHighRT)){
  
  SAOutVarConcDietChagasHigh<-SAOutVarConcDietChagasHighRT[[i]]
  ListAcetateDietChagasHighRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(y){y[,"s_ac"]})
  
  rm(SAOutVarConcDietChagasHigh)
  
}

AcetateDietChagasHighRT<-parLapply(cl=cl,X=ListAcetateDietChagasHighRT,
                                     function(x){
                                       ListAcetateDietChagasHigh<-x
                                       as.data.frame(do.call(rbind,ListAcetateDietChagasHigh))
                                       }
                                     )
for(i in 1:length(AcetateDietChagasHighRT)){
  
  row.names(AcetateDietChagasHighRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(AcetateDietChagasHighRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of acetate concentration dynamic simulated in output files
for(i in 1:length(AcetateDietChagasHighRT)){
  
  write.table(AcetateDietChagasHighRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasHighRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ListButyrateDietChagasHighRT<-vector("list",length=length(SAOutVarConcDietChagasHighRT))
for(i in 1:length(ListButyrateDietChagasHighRT)){
  
  SAOutVarConcDietChagasHigh<-SAOutVarConcDietChagasHighRT[[i]]
  ListButyrateDietChagasHighRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(y){y[,"s_bu"]})
  
  rm(SAOutVarConcDietChagasHigh)
  
}

ButyrateDietChagasHighRT<-parLapply(cl=cl,X=ListButyrateDietChagasHighRT,
                                      function(x){
                                        ListButyrateDietChagasHigh<-x
                                        as.data.frame(do.call(rbind,ListButyrateDietChagasHigh))
                                        }
                                      )

for(i in 1:length(ButyrateDietChagasHighRT)){
  
  row.names(ButyrateDietChagasHighRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(ButyrateDietChagasHighRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of butyrate concentration dynamic simulated in output files
for(i in 1:length(ButyrateDietChagasHighRT)){
  
  write.table(ButyrateDietChagasHighRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasHighRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
ListPropionateDietChagasHighRT<-vector("list",length=length(SAOutVarConcDietChagasHighRT))
for(i in 1:length(ListPropionateDietChagasHighRT)){
  
  SAOutVarConcDietChagasHigh<-SAOutVarConcDietChagasHighRT[[i]]
  ListPropionateDietChagasHighRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasHigh,function(y){y[,"s_pr"]})
  
  rm(SAOutVarConcDietChagasHigh)
  
}

PropionateDietChagasHighRT<-parLapply(cl=cl,X=ListPropionateDietChagasHighRT,
                                        function(x){
                                          ListPropionateDietChagasHigh<-x
                                          as.data.frame(do.call(rbind,ListPropionateDietChagasHigh))
                                          }
                                        )
for(i in 1:length(PropionateDietChagasHighRT)){
  
  row.names(PropionateDietChagasHighRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(PropionateDietChagasHighRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of Propionate concentration dynamic simulated in output files
for(i in 1:length(PropionateDietChagasHighRT)){
  
  write.table(PropionateDietChagasHighRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasHighRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

stopCluster(cl)