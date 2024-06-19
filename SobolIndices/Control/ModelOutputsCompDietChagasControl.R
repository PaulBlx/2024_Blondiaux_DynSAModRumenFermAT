# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the output variables of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied in the sensitivity analysis analysis using the sampling matrix of input parameters simulated for the control

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
XDietChagasControl<-Cinit_Control #Initial conditions of the system
auxsDietChagasControl<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrControl) #High forage diet and treatment control
XmDietChagasControlRT<-vector("list",length=length(listRTSamplesInputsPickandfreeze))
for(i in 1:length(XmDietChagasControlRT)){
  
  eval(parse(text=paste("SamplingMatrix<-listRTSamplesInputsPickandfreeze[[",i,"]]",sep="")))
  XmDietChagasControlRT[[i]]<-parApply(cl=cl,X=SamplingMatrix,MARGIN=1,FUN=ModelRun,Cinit=XDietChagasControl,t=t,auxs=auxsDietChagasControl) #Solving ode system
  
  rm(SamplingMatrix)
  
}

#Selection of the 4th day of simulation
XmDietChagasControlRTD4<-vector("list",length=length(XmDietChagasControlRT))
for(i in 1:length(XmDietChagasControlRTD4)){
  
  XmDietChagasControl<-XmDietChagasControlRT[[i]]
  XmDietChagasControlRTD4[[i]]<-parLapply(cl=cl,X=XmDietChagasControl,fun=function(x){return(x[RowSimd4,])})
  
  rm(XmDietChagasControl)
  
}
XmDietChagasControlRT<-XmDietChagasControlRTD4

#Selection of output variables of interest (methane and volatile fatty acids (VFA) concentrations)
SAOutVarConcDietChagasControlRT<-vector("list",length=length(XmDietChagasControlRT))
for(i in 1:length(SAOutVarConcDietChagasControlRT)){
  
  XmDietChagasControl<-XmDietChagasControlRT[[i]]
  SAOutVarConcDietChagasControlRT[[i]]<-parLapply(cl=cl,X=XmDietChagasControl,fun=rumencATout,V_l=V_l,V_g=V_g)
  
  rm(XmDietChagasControl)
  
}

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
ListMethaneDietChagasControlRT<-vector("list",length=length(SAOutVarConcDietChagasControlRT))
for(i in 1:length(ListMethaneDietChagasControlRT)){
  
  SAOutVarConcDietChagasControl<-SAOutVarConcDietChagasControlRT[[i]]
  ListMethaneDietChagasControlRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(y){y[,"q_ch4g_out"]})
  
  rm(SAOutVarConcDietChagasControl)
  
}

MethaneDietChagasControlRT<-parLapply(cl=cl,X=ListMethaneDietChagasControlRT,
                                        function(x){
                                          ListMethaneDietChagasControl<-x
                                          as.data.frame(do.call(rbind,ListMethaneDietChagasControl))
                                          }
                                        )
for(i in 1:length(MethaneDietChagasControlRT)){
  
  row.names(MethaneDietChagasControlRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(MethaneDietChagasControlRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of methane concentration dynamic simulated in output files
for(i in 1:length(MethaneDietChagasControlRT)){
  
  write.table(MethaneDietChagasControlRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasControlRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
ListAcetateDietChagasControlRT<-vector("list",length=length(SAOutVarConcDietChagasControlRT))
for(i in 1:length(ListAcetateDietChagasControlRT)){
  
  SAOutVarConcDietChagasControl<-SAOutVarConcDietChagasControlRT[[i]]
  ListAcetateDietChagasControlRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(y){y[,"s_ac"]})
  
  rm(SAOutVarConcDietChagasControl)
  
}

AcetateDietChagasControlRT<-parLapply(cl=cl,X=ListAcetateDietChagasControlRT,
                                        function(x){
                                          ListAcetateDietChagasControl<-x
                                          as.data.frame(do.call(rbind,ListAcetateDietChagasControl))
                                        }
)
for(i in 1:length(AcetateDietChagasControlRT)){
  
  row.names(AcetateDietChagasControlRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(AcetateDietChagasControlRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of acetate concentration dynamic simulated in output files
for(i in 1:length(AcetateDietChagasControlRT)){
  
  write.table(AcetateDietChagasControlRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasControlRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ListButyrateDietChagasControlRT<-vector("list",length=length(SAOutVarConcDietChagasControlRT))
for(i in 1:length(ListButyrateDietChagasControlRT)){
  
  SAOutVarConcDietChagasControl<-SAOutVarConcDietChagasControlRT[[i]]
  ListButyrateDietChagasControlRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(y){y[,"s_bu"]})
  
  rm(SAOutVarConcDietChagasControl)
  
}

ButyrateDietChagasControlRT<-parLapply(cl=cl,X=ListButyrateDietChagasControlRT,
                                         function(x){
                                           ListButyrateDietChagasControl<-x
                                           as.data.frame(do.call(rbind,ListButyrateDietChagasControl))
                                           }
                                         )
for(i in 1:length(ButyrateDietChagasControlRT)){
  
  row.names(ButyrateDietChagasControlRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(ButyrateDietChagasControlRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of butyrate concentration dynamic simulated in output files
for(i in 1:length(ButyrateDietChagasControlRT)){
  
  write.table(ButyrateDietChagasControlRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasControlRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
ListPropionateDietChagasControlRT<-vector("list",length=length(SAOutVarConcDietChagasControlRT))
for(i in 1:length(ListPropionateDietChagasControlRT)){
  
  SAOutVarConcDietChagasControl<-SAOutVarConcDietChagasControlRT[[i]]
  ListPropionateDietChagasControlRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasControl,function(y){y[,"s_pr"]})
  
  rm(SAOutVarConcDietChagasControl)
  
}

PropionateDietChagasControlRT<-parLapply(cl=cl,X=ListPropionateDietChagasControlRT,
                                           function(x){
                                             ListPropionateDietChagasControl<-x
                                             as.data.frame(do.call(rbind,ListPropionateDietChagasControl))
                                             }
                                           )
for(i in 1:length(PropionateDietChagasControlRT)){
  
  row.names(PropionateDietChagasControlRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(PropionateDietChagasControlRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of Propionate concentration dynamic simulated in output files
for(i in 1:length(PropionateDietChagasControlRT)){
  
  write.table(PropionateDietChagasControlRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasControlRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

stopCluster(cl)