# -----------------------------------------------------------------------------------------

# Computation of the full and independent Sobol indices (Mara et al., 2015) for the control

# -----------------------------------------------------------------------------------------

library(data.table)
library(boot)
library(parallel)
library(deSolve)

#Loading RData containing input parameter sampling matrix and functions for running mechanistic model and computing full and independent Sobol indices
load(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))

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

# ---- Estimation of full and independent Sobol indices ----

#Settings of the computation
R<-100 #Number of bootstrap replication
type<-"norm" #Normal method for confidence interval computation
conf<-0.95 #Confidence intervals level

#Estimation of full and independent Sobol indices using sobol (first-order) and saltelli (total) estimators
SobolindicesComputation<-function(x,RT){
  
  #################################################################################
  # List of inputs to this function:
  # x: list of output variable matrix for each permutation 
  # RT: Order of inputs
  #################################################################################
  
  #Indices computed for the reordering considered 
  FullIndicescomputed<-RT[1]
  IndIndicescomputed<-RT[length(RT)]
  
  #List for storing Sobol indice estimates
  ListSobolindices<-vector("list",length=2)
  names(ListSobolindices)<-c(FullIndicescomputed,IndIndicescomputed)
  ListSobolindices[[1]]<-data.frame(matrix(NA,nrow=2,ncol=ncol(x)))
  row.names(ListSobolindices[[1]])<-c("Sfull","Tfull")
  colnames(ListSobolindices[[1]])<-colnames(x)
  ListSobolindices[[2]]<-data.frame(matrix(NA,nrow=2,ncol=ncol(x)))
  row.names(ListSobolindices[[2]])<-c("Sind","Tind")
  colnames(ListSobolindices[[2]])<-colnames(x)
  
  #List for storing Sobol indice standard errors
  ListSobolindicesStdErr<-vector("list",length=2)
  names(ListSobolindicesStdErr)<-c(FullIndicescomputed,IndIndicescomputed)
  ListSobolindicesStdErr[[1]]<-data.frame(matrix(NA,nrow=2,ncol=ncol(x)))
  row.names(ListSobolindicesStdErr[[1]])<-c("Sfull","Tfull")
  colnames(ListSobolindicesStdErr[[1]])<-colnames(x)
  ListSobolindicesStdErr[[2]]<-data.frame(matrix(NA,nrow=2,ncol=ncol(x)))
  row.names(ListSobolindicesStdErr[[2]])<-c("Sind","Tind")
  colnames(ListSobolindicesStdErr[[2]])<-colnames(x)
  
  for(j in 1:ncol(x)){
    
    #Selection of one time step
    OutVarTimestep<-x[,j]
    
    #Computation of Sobol indices
    SobolindicesComputation<-sobol_indices(matrices=c("A","B","BA"),Y=OutVarTimestep,
                                           N=Nsim,params=RT,
                                           first="sobol",total="saltelli",
                                           boot=TRUE,R=R,parallel="no",
                                           conf=conf,type=type,
                                           full=FullIndicescomputed,
                                           indep=IndIndicescomputed)
    
    FullIndices<-SobolindicesComputation$results[c(which(SobolindicesComputation$results[,"sensitivity"]=="Sifull"),
                                                   which(SobolindicesComputation$results[,"sensitivity"]=="Tifull")),]
    IndIndices<-SobolindicesComputation$results[c(which(SobolindicesComputation$results[,"sensitivity"]=="Si-1ind"),
                                                  which(SobolindicesComputation$results[,"sensitivity"]=="Ti-1ind")),]
    
    #Estimate
    ListSobolindices[[FullIndicescomputed]][,j]<-FullIndices[,"Sobolindices.original"]
    ListSobolindices[[IndIndicescomputed]][,j]<-IndIndices[,"Sobolindices.original"]
    #Standard error
    ListSobolindicesStdErr[[FullIndicescomputed]][,j]<-FullIndices[,"Sobolindices.std.error"]
    ListSobolindicesStdErr[[IndIndicescomputed]][,j]<-IndIndices[,"Sobolindices.std.error"]
    
    rm(OutVarTimestep);rm(SobolindicesComputation);rm(FullIndices);rm(IndIndices)
    
  }
  
  return(list(Sobolindices=ListSobolindices,SobolindicesStdErr=ListSobolindicesStdErr))
  
} 

# ---- Computing the full and independent Sobol indices ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
    
    PathMethaneDietChagasControlRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasControlRT",i,".txt",sep="")
    eval(parse(text=paste("MethaneDietChagasControlRT",i,"<-read.table(PathMethaneDietChagasControlRT,header=TRUE)",sep="")))
    
}

#List of the different outputs of the Rosenblatt transformation
MethaneDietChagasControlRT<-list(MethaneDietChagasControlRT1,MethaneDietChagasControlRT2,
                                 MethaneDietChagasControlRT3,MethaneDietChagasControlRT4,
                                 MethaneDietChagasControlRT5,MethaneDietChagasControlRT6,
                                 MethaneDietChagasControlRT7,MethaneDietChagasControlRT8,
                                 MethaneDietChagasControlRT9,MethaneDietChagasControlRT10,
                                 MethaneDietChagasControlRT11,MethaneDietChagasControlRT12,
                                 MethaneDietChagasControlRT13,MethaneDietChagasControlRT14,
                                 MethaneDietChagasControlRT15,MethaneDietChagasControlRT16
                                 )

#Full and independent Sobol indices computation
ListSobolindicesMethaneDietChagasControlRT<-vector("list",length=length(MethaneDietChagasControlRT))
for(i in 1:length(ListSobolindicesMethaneDietChagasControlRT)){
  
  MethaneDietChagasControl<-MethaneDietChagasControlRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesMethaneDietChagasControlRT[[i]]<-SobolindicesComputation(x=MethaneDietChagasControl,RT=InputsRTi)
  
  rm(MethaneDietChagasControl);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullMethaneDietChagasControl)<-Inputs
colnames(SfullMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Total indice
TfullMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullMethaneDietChagasControl)<-Inputs
colnames(TfullMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindMethaneDietChagasControl)<-Inputs
colnames(SindMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Total indice
TindMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindMethaneDietChagasControl)<-Inputs
colnames(TindMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullMethaneDietChagasControl[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullMethaneDietChagasControl[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindMethaneDietChagasControl[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindMethaneDietChagasControl[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrMethaneDietChagasControl)<-Inputs
colnames(SfullStdErrMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Total indice
TfullStdErrMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrMethaneDietChagasControl)<-Inputs
colnames(TfullStdErrMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindStdErrMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrMethaneDietChagasControl)<-Inputs
colnames(SindStdErrMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])
#Total indice
TindStdErrMethaneDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrMethaneDietChagasControl)<-Inputs
colnames(TindStdErrMethaneDietChagasControl)<-colnames(MethaneDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrMethaneDietChagasControl[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrMethaneDietChagasControl[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrMethaneDietChagasControl[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrMethaneDietChagasControl[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullStdErrMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullStdErrMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindStdErrMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrMethaneDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindStdErrMethaneDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathAcetateDietChagasControlRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasControlRT",i,".txt",sep="")
  eval(parse(text=paste("AcetateDietChagasControlRT",i,"<-read.table(PathAcetateDietChagasControlRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
AcetateDietChagasControlRT<-list(AcetateDietChagasControlRT1,AcetateDietChagasControlRT2,
                                 AcetateDietChagasControlRT3,AcetateDietChagasControlRT4,
                                 AcetateDietChagasControlRT5,AcetateDietChagasControlRT6,
                                 AcetateDietChagasControlRT7,AcetateDietChagasControlRT8,
                                 AcetateDietChagasControlRT9,AcetateDietChagasControlRT10,
                                 AcetateDietChagasControlRT11,AcetateDietChagasControlRT12,
                                 AcetateDietChagasControlRT13,AcetateDietChagasControlRT14,
                                 AcetateDietChagasControlRT15,AcetateDietChagasControlRT16
)

#Full and independent Sobol indices computation
ListSobolindicesAcetateDietChagasControlRT<-vector("list",length=length(AcetateDietChagasControlRT))
for(i in 1:length(ListSobolindicesAcetateDietChagasControlRT)){
  
  AcetateDietChagasControl<-AcetateDietChagasControlRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesAcetateDietChagasControlRT[[i]]<-SobolindicesComputation(x=AcetateDietChagasControl,RT=InputsRTi)
  
  rm(AcetateDietChagasControl);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullAcetateDietChagasControl)<-Inputs
colnames(SfullAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Total indice
TfullAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullAcetateDietChagasControl)<-Inputs
colnames(TfullAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindAcetateDietChagasControl)<-Inputs
colnames(SindAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Total indice
TindAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindAcetateDietChagasControl)<-Inputs
colnames(TindAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullAcetateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullAcetateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindAcetateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindAcetateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrAcetateDietChagasControl)<-Inputs
colnames(SfullStdErrAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Total indice
TfullStdErrAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrAcetateDietChagasControl)<-Inputs
colnames(TfullStdErrAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindStdErrAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrAcetateDietChagasControl)<-Inputs
colnames(SindStdErrAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])
#Total indice
TindStdErrAcetateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrAcetateDietChagasControl)<-Inputs
colnames(TindStdErrAcetateDietChagasControl)<-colnames(AcetateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrAcetateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrAcetateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrAcetateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrAcetateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullStdErrAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullStdErrAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindStdErrAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrAcetateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindStdErrAcetateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathButyrateDietChagasControlRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasControlRT",i,".txt",sep="")
  eval(parse(text=paste("ButyrateDietChagasControlRT",i,"<-read.table(PathButyrateDietChagasControlRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
ButyrateDietChagasControlRT<-list(ButyrateDietChagasControlRT1,ButyrateDietChagasControlRT2,
                                  ButyrateDietChagasControlRT3,ButyrateDietChagasControlRT4,
                                  ButyrateDietChagasControlRT5,ButyrateDietChagasControlRT6,
                                  ButyrateDietChagasControlRT7,ButyrateDietChagasControlRT8,
                                  ButyrateDietChagasControlRT9,ButyrateDietChagasControlRT10,
                                  ButyrateDietChagasControlRT11,ButyrateDietChagasControlRT12,
                                  ButyrateDietChagasControlRT13,ButyrateDietChagasControlRT14,
                                  ButyrateDietChagasControlRT15,ButyrateDietChagasControlRT16
                                  )

#Full and independent Sobol indices computation
ListSobolindicesButyrateDietChagasControlRT<-vector("list",length=length(ButyrateDietChagasControlRT))
for(i in 1:length(ListSobolindicesButyrateDietChagasControlRT)){
  
  ButyrateDietChagasControl<-ButyrateDietChagasControlRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesButyrateDietChagasControlRT[[i]]<-SobolindicesComputation(x=ButyrateDietChagasControl,RT=InputsRTi)
  
  rm(ButyrateDietChagasControl);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullButyrateDietChagasControl)<-Inputs
colnames(SfullButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Total indice
TfullButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullButyrateDietChagasControl)<-Inputs
colnames(TfullButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindButyrateDietChagasControl)<-Inputs
colnames(SindButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Total indice
TindButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindButyrateDietChagasControl)<-Inputs
colnames(TindButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullButyrateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullButyrateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindButyrateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindButyrateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrButyrateDietChagasControl)<-Inputs
colnames(SfullStdErrButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Total indice
TfullStdErrButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrButyrateDietChagasControl)<-Inputs
colnames(TfullStdErrButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindStdErrButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrButyrateDietChagasControl)<-Inputs
colnames(SindStdErrButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])
#Total indice
TindStdErrButyrateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrButyrateDietChagasControl)<-Inputs
colnames(TindStdErrButyrateDietChagasControl)<-colnames(ButyrateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrButyrateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrButyrateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrButyrateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrButyrateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullStdErrButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullStdErrButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindStdErrButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrButyrateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindStdErrButyrateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathPropionateDietChagasControlRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasControlRT",i,".txt",sep="")
  eval(parse(text=paste("PropionateDietChagasControlRT",i,"<-read.table(PathPropionateDietChagasControlRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
PropionateDietChagasControlRT<-list(PropionateDietChagasControlRT1,PropionateDietChagasControlRT2,
                                    PropionateDietChagasControlRT3,PropionateDietChagasControlRT4,
                                    PropionateDietChagasControlRT5,PropionateDietChagasControlRT6,
                                    PropionateDietChagasControlRT7,PropionateDietChagasControlRT8,
                                    PropionateDietChagasControlRT9,PropionateDietChagasControlRT10,
                                    PropionateDietChagasControlRT11,PropionateDietChagasControlRT12,
                                    PropionateDietChagasControlRT13,PropionateDietChagasControlRT14,
                                    PropionateDietChagasControlRT15,PropionateDietChagasControlRT16
                                    )

#Full and independent Sobol indices computation
ListSobolindicesPropionateDietChagasControlRT<-vector("list",length=length(PropionateDietChagasControlRT))
for(i in 1:length(ListSobolindicesPropionateDietChagasControlRT)){
  
  PropionateDietChagasControl<-PropionateDietChagasControlRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesPropionateDietChagasControlRT[[i]]<-SobolindicesComputation(x=PropionateDietChagasControl,RT=InputsRTi)
  
  rm(PropionateDietChagasControl);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullPropionateDietChagasControl)<-Inputs
colnames(SfullPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Total indice
TfullPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullPropionateDietChagasControl)<-Inputs
colnames(TfullPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindPropionateDietChagasControl)<-Inputs
colnames(SindPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Total indice
TindPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindPropionateDietChagasControl)<-Inputs
colnames(TindPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullPropionateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullPropionateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindPropionateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindPropionateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrPropionateDietChagasControl)<-Inputs
colnames(SfullStdErrPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Total indice
TfullStdErrPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrPropionateDietChagasControl)<-Inputs
colnames(TfullStdErrPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Independent indices
#First-order indice
SindStdErrPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrPropionateDietChagasControl)<-Inputs
colnames(SindStdErrPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])
#Total indice
TindStdErrPropionateDietChagasControl<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrPropionateDietChagasControl)<-Inputs
colnames(TindStdErrPropionateDietChagasControl)<-colnames(PropionateDietChagasControlRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasControlRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrPropionateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrPropionateDietChagasControl[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrPropionateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrPropionateDietChagasControl[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasControlRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullStdErrPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullStdErrPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindStdErrPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrPropionateDietChagasControl,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindStdErrPropionateDietChagasControl.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl)