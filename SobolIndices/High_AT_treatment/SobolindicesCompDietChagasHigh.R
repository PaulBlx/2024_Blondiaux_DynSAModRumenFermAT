# ------------------------------------------------------------------------------------------------------------------------

# Computation of the full and independent Sobol indices (Mara et al., 2015) for the high Asparagopsis taxiformis treatment

# ------------------------------------------------------------------------------------------------------------------------

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
  
  PathMethaneDietChagasHighRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasHighRT",i,".txt",sep="")
  eval(parse(text=paste("MethaneDietChagasHighRT",i,"<-read.table(PathMethaneDietChagasHighRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
MethaneDietChagasHighRT<-list(MethaneDietChagasHighRT1,MethaneDietChagasHighRT2,
                              MethaneDietChagasHighRT3,MethaneDietChagasHighRT4,
                              MethaneDietChagasHighRT5,MethaneDietChagasHighRT6,
                              MethaneDietChagasHighRT7,MethaneDietChagasHighRT8,
                              MethaneDietChagasHighRT9,MethaneDietChagasHighRT10,
                              MethaneDietChagasHighRT11,MethaneDietChagasHighRT12,
                              MethaneDietChagasHighRT13,MethaneDietChagasHighRT14,
                              MethaneDietChagasHighRT15,MethaneDietChagasHighRT16
                              )

#Full and independent Sobol indices computation
ListSobolindicesMethaneDietChagasHighRT<-vector("list",length=length(MethaneDietChagasHighRT))
for(i in 1:length(ListSobolindicesMethaneDietChagasHighRT)){
  
  MethaneDietChagasHigh<-MethaneDietChagasHighRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesMethaneDietChagasHighRT[[i]]<-SobolindicesComputation(x=MethaneDietChagasHigh,RT=InputsRTi)
  
  rm(MethaneDietChagasHigh);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullMethaneDietChagasHigh)<-Inputs
colnames(SfullMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Total indice
TfullMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullMethaneDietChagasHigh)<-Inputs
colnames(TfullMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindMethaneDietChagasHigh)<-Inputs
colnames(SindMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Total indice
TindMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindMethaneDietChagasHigh)<-Inputs
colnames(TindMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullMethaneDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullMethaneDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindMethaneDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindMethaneDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrMethaneDietChagasHigh)<-Inputs
colnames(SfullStdErrMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Total indice
TfullStdErrMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrMethaneDietChagasHigh)<-Inputs
colnames(TfullStdErrMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindStdErrMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrMethaneDietChagasHigh)<-Inputs
colnames(SindStdErrMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])
#Total indice
TindStdErrMethaneDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrMethaneDietChagasHigh)<-Inputs
colnames(TindStdErrMethaneDietChagasHigh)<-colnames(MethaneDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrMethaneDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrMethaneDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrMethaneDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrMethaneDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullStdErrMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullStdErrMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindStdErrMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrMethaneDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindStdErrMethaneDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathAcetateDietChagasHighRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasHighRT",i,".txt",sep="")
  eval(parse(text=paste("AcetateDietChagasHighRT",i,"<-read.table(PathAcetateDietChagasHighRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
AcetateDietChagasHighRT<-list(AcetateDietChagasHighRT1,AcetateDietChagasHighRT2,
                              AcetateDietChagasHighRT3,AcetateDietChagasHighRT4,
                              AcetateDietChagasHighRT5,AcetateDietChagasHighRT6,
                              AcetateDietChagasHighRT7,AcetateDietChagasHighRT8,
                              AcetateDietChagasHighRT9,AcetateDietChagasHighRT10,
                              AcetateDietChagasHighRT11,AcetateDietChagasHighRT12,
                              AcetateDietChagasHighRT13,AcetateDietChagasHighRT14,
                              AcetateDietChagasHighRT15,AcetateDietChagasHighRT16
                              )

#Full and independent Sobol indices computation
ListSobolindicesAcetateDietChagasHighRT<-vector("list",length=length(AcetateDietChagasHighRT))
for(i in 1:length(ListSobolindicesAcetateDietChagasHighRT)){
  
  AcetateDietChagasHigh<-AcetateDietChagasHighRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesAcetateDietChagasHighRT[[i]]<-SobolindicesComputation(x=AcetateDietChagasHigh,RT=InputsRTi)
  
  rm(AcetateDietChagasHigh);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullAcetateDietChagasHigh)<-Inputs
colnames(SfullAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Total indice
TfullAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullAcetateDietChagasHigh)<-Inputs
colnames(TfullAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindAcetateDietChagasHigh)<-Inputs
colnames(SindAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Total indice
TindAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindAcetateDietChagasHigh)<-Inputs
colnames(TindAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullAcetateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullAcetateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindAcetateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindAcetateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrAcetateDietChagasHigh)<-Inputs
colnames(SfullStdErrAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Total indice
TfullStdErrAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrAcetateDietChagasHigh)<-Inputs
colnames(TfullStdErrAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindStdErrAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrAcetateDietChagasHigh)<-Inputs
colnames(SindStdErrAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])
#Total indice
TindStdErrAcetateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrAcetateDietChagasHigh)<-Inputs
colnames(TindStdErrAcetateDietChagasHigh)<-colnames(AcetateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrAcetateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrAcetateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrAcetateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrAcetateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullStdErrAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullStdErrAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindStdErrAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrAcetateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindStdErrAcetateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathButyrateDietChagasHighRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasHighRT",i,".txt",sep="")
  eval(parse(text=paste("ButyrateDietChagasHighRT",i,"<-read.table(PathButyrateDietChagasHighRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
ButyrateDietChagasHighRT<-list(ButyrateDietChagasHighRT1,ButyrateDietChagasHighRT2,
                               ButyrateDietChagasHighRT3,ButyrateDietChagasHighRT4,
                               ButyrateDietChagasHighRT5,ButyrateDietChagasHighRT6,
                               ButyrateDietChagasHighRT7,ButyrateDietChagasHighRT8,
                               ButyrateDietChagasHighRT9,ButyrateDietChagasHighRT10,
                               ButyrateDietChagasHighRT11,ButyrateDietChagasHighRT12,
                               ButyrateDietChagasHighRT13,ButyrateDietChagasHighRT14,
                               ButyrateDietChagasHighRT15,ButyrateDietChagasHighRT16
                               )

#Full and independent Sobol indices computation
ListSobolindicesButyrateDietChagasHighRT<-vector("list",length=length(ButyrateDietChagasHighRT))
for(i in 1:length(ListSobolindicesButyrateDietChagasHighRT)){
  
  ButyrateDietChagasHigh<-ButyrateDietChagasHighRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesButyrateDietChagasHighRT[[i]]<-SobolindicesComputation(x=ButyrateDietChagasHigh,RT=InputsRTi)
  
  rm(ButyrateDietChagasHigh);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullButyrateDietChagasHigh)<-Inputs
colnames(SfullButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Total indice
TfullButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullButyrateDietChagasHigh)<-Inputs
colnames(TfullButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindButyrateDietChagasHigh)<-Inputs
colnames(SindButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Total indice
TindButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindButyrateDietChagasHigh)<-Inputs
colnames(TindButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullButyrateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullButyrateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindButyrateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindButyrateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrButyrateDietChagasHigh)<-Inputs
colnames(SfullStdErrButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Total indice
TfullStdErrButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrButyrateDietChagasHigh)<-Inputs
colnames(TfullStdErrButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindStdErrButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrButyrateDietChagasHigh)<-Inputs
colnames(SindStdErrButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])
#Total indice
TindStdErrButyrateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrButyrateDietChagasHigh)<-Inputs
colnames(TindStdErrButyrateDietChagasHigh)<-colnames(ButyrateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrButyrateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrButyrateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrButyrateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrButyrateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullStdErrButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullStdErrButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindStdErrButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrButyrateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindStdErrButyrateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathPropionateDietChagasHighRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasHighRT",i,".txt",sep="")
  eval(parse(text=paste("PropionateDietChagasHighRT",i,"<-read.table(PathPropionateDietChagasHighRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
PropionateDietChagasHighRT<-list(PropionateDietChagasHighRT1,PropionateDietChagasHighRT2,
                                 PropionateDietChagasHighRT3,PropionateDietChagasHighRT4,
                                 PropionateDietChagasHighRT5,PropionateDietChagasHighRT6,
                                 PropionateDietChagasHighRT7,PropionateDietChagasHighRT8,
                                 PropionateDietChagasHighRT9,PropionateDietChagasHighRT10,
                                 PropionateDietChagasHighRT11,PropionateDietChagasHighRT12,
                                 PropionateDietChagasHighRT13,PropionateDietChagasHighRT14,
                                 PropionateDietChagasHighRT15,PropionateDietChagasHighRT16
                                 )

#Full and independent Sobol indices computation
ListSobolindicesPropionateDietChagasHighRT<-vector("list",length=length(PropionateDietChagasHighRT))
for(i in 1:length(ListSobolindicesPropionateDietChagasHighRT)){
  
  PropionateDietChagasHigh<-PropionateDietChagasHighRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesPropionateDietChagasHighRT[[i]]<-SobolindicesComputation(x=PropionateDietChagasHigh,RT=InputsRTi)
  
  rm(PropionateDietChagasHigh);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullPropionateDietChagasHigh)<-Inputs
colnames(SfullPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Total indice
TfullPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullPropionateDietChagasHigh)<-Inputs
colnames(TfullPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindPropionateDietChagasHigh)<-Inputs
colnames(SindPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Total indice
TindPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindPropionateDietChagasHigh)<-Inputs
colnames(TindPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullPropionateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullPropionateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindPropionateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindPropionateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrPropionateDietChagasHigh)<-Inputs
colnames(SfullStdErrPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Total indice
TfullStdErrPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrPropionateDietChagasHigh)<-Inputs
colnames(TfullStdErrPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Independent indices
#First-order indice
SindStdErrPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrPropionateDietChagasHigh)<-Inputs
colnames(SindStdErrPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])
#Total indice
TindStdErrPropionateDietChagasHigh<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrPropionateDietChagasHigh)<-Inputs
colnames(TindStdErrPropionateDietChagasHigh)<-colnames(PropionateDietChagasHighRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasHighRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrPropionateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrPropionateDietChagasHigh[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrPropionateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrPropionateDietChagasHigh[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasHighRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullStdErrPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullStdErrPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindStdErrPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrPropionateDietChagasHigh,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindStdErrPropionateDietChagasHigh.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl)