# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Computation of the output variables of the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021) studied in the sensitivity analysis analysis using the sampling matrix of input parameters simulated for the low Asparagopsis taxiformis treatment

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
XDietChagasLow<-Cinit_Control #Initial conditions of the system
auxsDietChagasLow<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrLow) #High forage diet and treatment low
XmDietChagasLowRT<-vector("list",length=length(listRTSamplesInputsPickandfreeze))
for(i in 1:length(XmDietChagasLowRT)){
  
  eval(parse(text=paste("SamplingMatrix<-listRTSamplesInputsPickandfreeze[[",i,"]]",sep="")))
  XmDietChagasLowRT[[i]]<-parApply(cl=cl,X=SamplingMatrix,MARGIN=1,FUN=ModelRun,Cinit=XDietChagasLow,t=t,auxs=auxsDietChagasLow) #Solving ode system
  
  rm(SamplingMatrix)
  
}

#Selection of the 4th day of simulation
XmDietChagasLowRTD4<-vector("list",length=length(XmDietChagasLowRT))
for(i in 1:length(XmDietChagasLowRTD4)){
  
  XmDietChagasLow<-XmDietChagasLowRT[[i]]
  XmDietChagasLowRTD4[[i]]<-parLapply(cl=cl,X=XmDietChagasLow,fun=function(x){return(x[RowSimd4,])})
  
  rm(XmDietChagasLow)
  
}
XmDietChagasLowRT<-XmDietChagasLowRTD4

#Selection of output variables of interest (methane and volatile fatty acids (VFA) concentrations)
SAOutVarConcDietChagasLowRT<-vector("list",length=length(XmDietChagasLowRT))
for(i in 1:length(SAOutVarConcDietChagasLowRT)){
  
  XmDietChagasLow<-XmDietChagasLowRT[[i]]
  SAOutVarConcDietChagasLowRT[[i]]<-parLapply(cl=cl,X=XmDietChagasLow,fun=rumencATout,V_l=V_l,V_g=V_g)
  
  rm(XmDietChagasLow)
  
}

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
ListMethaneDietChagasLowRT<-vector("list",length=length(SAOutVarConcDietChagasLowRT))
for(i in 1:length(ListMethaneDietChagasLowRT)){
  
  SAOutVarConcDietChagasLow<-SAOutVarConcDietChagasLowRT[[i]]
  ListMethaneDietChagasLowRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasLow,function(y){y[,"q_ch4g_out"]})
  
  rm(SAOutVarConcDietChagasLow)
  
}

MethaneDietChagasLowRT<-parLapply(cl=cl,X=ListMethaneDietChagasLowRT,
                                        function(x){
                                          ListMethaneDietChagasLow<-x
                                          as.data.frame(do.call(rbind,ListMethaneDietChagasLow))
                                          }
                                        )
for(i in 1:length(MethaneDietChagasLowRT)){
  
  row.names(MethaneDietChagasLowRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(MethaneDietChagasLowRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of methane concentration dynamic simulated in output files
for(i in 1:length(MethaneDietChagasLowRT)){
  
  write.table(MethaneDietChagasLowRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasLowRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
ListAcetateDietChagasLowRT<-vector("list",length=length(SAOutVarConcDietChagasLowRT))
for(i in 1:length(ListAcetateDietChagasLowRT)){
  
  SAOutVarConcDietChagasLow<-SAOutVarConcDietChagasLowRT[[i]]
  ListAcetateDietChagasLowRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasLow,function(y){y[,"s_ac"]})
  
  rm(SAOutVarConcDietChagasLow)
  
}

AcetateDietChagasLowRT<-parLapply(cl=cl,X=ListAcetateDietChagasLowRT,
                                        function(x){
                                          ListAcetateDietChagasLow<-x
                                          as.data.frame(do.call(rbind,ListAcetateDietChagasLow))
                                        }
)
for(i in 1:length(AcetateDietChagasLowRT)){
  
  row.names(AcetateDietChagasLowRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(AcetateDietChagasLowRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of acetate concentration dynamic simulated in output files
for(i in 1:length(AcetateDietChagasLowRT)){
  
  write.table(AcetateDietChagasLowRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasLowRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ListButyrateDietChagasLowRT<-vector("list",length=length(SAOutVarConcDietChagasLowRT))
for(i in 1:length(ListButyrateDietChagasLowRT)){
  
  SAOutVarConcDietChagasLow<-SAOutVarConcDietChagasLowRT[[i]]
  ListButyrateDietChagasLowRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasLow,function(y){y[,"s_bu"]})
  
  rm(SAOutVarConcDietChagasLow)
  
}

ButyrateDietChagasLowRT<-parLapply(cl=cl,X=ListButyrateDietChagasLowRT,
                                         function(x){
                                           ListButyrateDietChagasLow<-x
                                           as.data.frame(do.call(rbind,ListButyrateDietChagasLow))
                                           }
                                         )
for(i in 1:length(ButyrateDietChagasLowRT)){
  
  row.names(ButyrateDietChagasLowRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(ButyrateDietChagasLowRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of butyrate concentration dynamic simulated in output files
for(i in 1:length(ButyrateDietChagasLowRT)){
  
  write.table(ButyrateDietChagasLowRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasLowRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
ListPropionateDietChagasLowRT<-vector("list",length=length(SAOutVarConcDietChagasLowRT))
for(i in 1:length(ListPropionateDietChagasLowRT)){
  
  SAOutVarConcDietChagasLow<-SAOutVarConcDietChagasLowRT[[i]]
  ListPropionateDietChagasLowRT[[i]]<-parLapply(cl=cl,X=SAOutVarConcDietChagasLow,function(y){y[,"s_pr"]})
  
  rm(SAOutVarConcDietChagasLow)
  
}

PropionateDietChagasLowRT<-parLapply(cl=cl,X=ListPropionateDietChagasLowRT,
                                           function(x){
                                             ListPropionateDietChagasLow<-x
                                             as.data.frame(do.call(rbind,ListPropionateDietChagasLow))
                                             }
                                           )
for(i in 1:length(PropionateDietChagasLowRT)){
  
  row.names(PropionateDietChagasLowRT[[i]])<-paste("Simulation",seq(1,eval(parse(text=paste("nrow(listRTSamplesInputsPickandfreeze[[",i,"]])",sep=""))),by=1),sep="")
  colnames(PropionateDietChagasLowRT[[i]])<-paste("t=",round(t[RowSimd4],3),"h",sep="")
  
}

#Integration of Propionate concentration dynamic simulated in output files
for(i in 1:length(PropionateDietChagasLowRT)){
  
  write.table(PropionateDietChagasLowRT[[i]][,c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
              file=paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasLowRT",i,".txt",sep=""),
              row.names=TRUE,col.names=TRUE)
  
}

stopCluster(cl)