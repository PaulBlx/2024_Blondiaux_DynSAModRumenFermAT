# -----------------------------------------------------------------------------------------------------------------------

# Computation of the full and independent Sobol indices (Mara et al., 2015) for the low Asparagopsis taxiformis treatment

# -----------------------------------------------------------------------------------------------------------------------

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
  
  PathMethaneDietChagasLowRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Methane/MethaneDietChagasLowRT",i,".txt",sep="")
  eval(parse(text=paste("MethaneDietChagasLowRT",i,"<-read.table(PathMethaneDietChagasLowRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
MethaneDietChagasLowRT<-list(MethaneDietChagasLowRT1,MethaneDietChagasLowRT2,
                             MethaneDietChagasLowRT3,MethaneDietChagasLowRT4,
                             MethaneDietChagasLowRT5,MethaneDietChagasLowRT6,
                             MethaneDietChagasLowRT7,MethaneDietChagasLowRT8,
                             MethaneDietChagasLowRT9,MethaneDietChagasLowRT10,
                             MethaneDietChagasLowRT11,MethaneDietChagasLowRT12,
                             MethaneDietChagasLowRT13,MethaneDietChagasLowRT14,
                             MethaneDietChagasLowRT15,MethaneDietChagasLowRT16
                             )

#Full and independent Sobol indices computation
ListSobolindicesMethaneDietChagasLowRT<-vector("list",length=length(MethaneDietChagasLowRT))
for(i in 1:length(ListSobolindicesMethaneDietChagasLowRT)){
  
  MethaneDietChagasLow<-MethaneDietChagasLowRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesMethaneDietChagasLowRT[[i]]<-SobolindicesComputation(x=MethaneDietChagasLow,RT=InputsRTi)
  
  rm(MethaneDietChagasLow);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullMethaneDietChagasLow)<-Inputs
colnames(SfullMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Total indice
TfullMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullMethaneDietChagasLow)<-Inputs
colnames(TfullMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindMethaneDietChagasLow)<-Inputs
colnames(SindMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Total indice
TindMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindMethaneDietChagasLow)<-Inputs
colnames(TindMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullMethaneDietChagasLow[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullMethaneDietChagasLow[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindMethaneDietChagasLow[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindMethaneDietChagasLow[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrMethaneDietChagasLow)<-Inputs
colnames(SfullStdErrMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Total indice
TfullStdErrMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrMethaneDietChagasLow)<-Inputs
colnames(TfullStdErrMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindStdErrMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrMethaneDietChagasLow)<-Inputs
colnames(SindStdErrMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])
#Total indice
TindStdErrMethaneDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrMethaneDietChagasLow)<-Inputs
colnames(TindStdErrMethaneDietChagasLow)<-colnames(MethaneDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesMethaneDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrMethaneDietChagasLow[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrMethaneDietChagasLow[FullIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrMethaneDietChagasLow[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrMethaneDietChagasLow[IndIndicescomputed,]<-ListSobolindicesMethaneDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SfullStdErrMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TfullStdErrMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/SindStdErrMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrMethaneDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Methane/TindStdErrMethaneDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathAcetateDietChagasLowRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Acetate/AcetateDietChagasLowRT",i,".txt",sep="")
  eval(parse(text=paste("AcetateDietChagasLowRT",i,"<-read.table(PathAcetateDietChagasLowRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
AcetateDietChagasLowRT<-list(AcetateDietChagasLowRT1,AcetateDietChagasLowRT2,
                             AcetateDietChagasLowRT3,AcetateDietChagasLowRT4,
                             AcetateDietChagasLowRT5,AcetateDietChagasLowRT6,
                             AcetateDietChagasLowRT7,AcetateDietChagasLowRT8,
                             AcetateDietChagasLowRT9,AcetateDietChagasLowRT10,
                             AcetateDietChagasLowRT11,AcetateDietChagasLowRT12,
                             AcetateDietChagasLowRT13,AcetateDietChagasLowRT14,
                             AcetateDietChagasLowRT15,AcetateDietChagasLowRT16
                             )

#Full and independent Sobol indices computation
ListSobolindicesAcetateDietChagasLowRT<-vector("list",length=length(AcetateDietChagasLowRT))
for(i in 1:length(ListSobolindicesAcetateDietChagasLowRT)){
  
  AcetateDietChagasLow<-AcetateDietChagasLowRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesAcetateDietChagasLowRT[[i]]<-SobolindicesComputation(x=AcetateDietChagasLow,RT=InputsRTi)
  
  rm(AcetateDietChagasLow);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullAcetateDietChagasLow)<-Inputs
colnames(SfullAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Total indice
TfullAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullAcetateDietChagasLow)<-Inputs
colnames(TfullAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindAcetateDietChagasLow)<-Inputs
colnames(SindAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Total indice
TindAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindAcetateDietChagasLow)<-Inputs
colnames(TindAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullAcetateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullAcetateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindAcetateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindAcetateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrAcetateDietChagasLow)<-Inputs
colnames(SfullStdErrAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Total indice
TfullStdErrAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrAcetateDietChagasLow)<-Inputs
colnames(TfullStdErrAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindStdErrAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrAcetateDietChagasLow)<-Inputs
colnames(SindStdErrAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])
#Total indice
TindStdErrAcetateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrAcetateDietChagasLow)<-Inputs
colnames(TindStdErrAcetateDietChagasLow)<-colnames(AcetateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesAcetateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrAcetateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrAcetateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrAcetateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrAcetateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesAcetateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SfullStdErrAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TfullStdErrAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/SindStdErrAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrAcetateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Acetate/TindStdErrAcetateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathButyrateDietChagasLowRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Butyrate/ButyrateDietChagasLowRT",i,".txt",sep="")
  eval(parse(text=paste("ButyrateDietChagasLowRT",i,"<-read.table(PathButyrateDietChagasLowRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
ButyrateDietChagasLowRT<-list(ButyrateDietChagasLowRT1,ButyrateDietChagasLowRT2,
                              ButyrateDietChagasLowRT3,ButyrateDietChagasLowRT4,
                              ButyrateDietChagasLowRT5,ButyrateDietChagasLowRT6,
                              ButyrateDietChagasLowRT7,ButyrateDietChagasLowRT8,
                              ButyrateDietChagasLowRT9,ButyrateDietChagasLowRT10,
                              ButyrateDietChagasLowRT11,ButyrateDietChagasLowRT12,
                              ButyrateDietChagasLowRT13,ButyrateDietChagasLowRT14,
                              ButyrateDietChagasLowRT15,ButyrateDietChagasLowRT16
                              )

#Full and independent Sobol indices computation
ListSobolindicesButyrateDietChagasLowRT<-vector("list",length=length(ButyrateDietChagasLowRT))
for(i in 1:length(ListSobolindicesButyrateDietChagasLowRT)){
  
  ButyrateDietChagasLow<-ButyrateDietChagasLowRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesButyrateDietChagasLowRT[[i]]<-SobolindicesComputation(x=ButyrateDietChagasLow,RT=InputsRTi)
  
  rm(ButyrateDietChagasLow);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullButyrateDietChagasLow)<-Inputs
colnames(SfullButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Total indice
TfullButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullButyrateDietChagasLow)<-Inputs
colnames(TfullButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindButyrateDietChagasLow)<-Inputs
colnames(SindButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Total indice
TindButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindButyrateDietChagasLow)<-Inputs
colnames(TindButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullButyrateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullButyrateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindButyrateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindButyrateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrButyrateDietChagasLow)<-Inputs
colnames(SfullStdErrButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Total indice
TfullStdErrButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrButyrateDietChagasLow)<-Inputs
colnames(TfullStdErrButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindStdErrButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrButyrateDietChagasLow)<-Inputs
colnames(SindStdErrButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])
#Total indice
TindStdErrButyrateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrButyrateDietChagasLow)<-Inputs
colnames(TindStdErrButyrateDietChagasLow)<-colnames(ButyrateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesButyrateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrButyrateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrButyrateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrButyrateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrButyrateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesButyrateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SfullStdErrButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TfullStdErrButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/SindStdErrButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrButyrateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Butyrate/TindStdErrButyrateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
for(i in 1:length(listRTSamplesInputsPickandfreeze)){
  
  PathPropionateDietChagasLowRT<-paste(getwd(),"/ModRumFermATOutVar_txt/Propionate/PropionateDietChagasLowRT",i,".txt",sep="")
  eval(parse(text=paste("PropionateDietChagasLowRT",i,"<-read.table(PathPropionateDietChagasLowRT,header=TRUE)",sep="")))
  
}

#List of the different outputs of the Rosenblatt transformation
PropionateDietChagasLowRT<-list(PropionateDietChagasLowRT1,PropionateDietChagasLowRT2,
                                PropionateDietChagasLowRT3,PropionateDietChagasLowRT4,
                                PropionateDietChagasLowRT5,PropionateDietChagasLowRT6,
                                PropionateDietChagasLowRT7,PropionateDietChagasLowRT8,
                                PropionateDietChagasLowRT9,PropionateDietChagasLowRT10,
                                PropionateDietChagasLowRT11,PropionateDietChagasLowRT12,
                                PropionateDietChagasLowRT13,PropionateDietChagasLowRT14,
                                PropionateDietChagasLowRT15,PropionateDietChagasLowRT16
                                )

#Full and independent Sobol indices computation
ListSobolindicesPropionateDietChagasLowRT<-vector("list",length=length(PropionateDietChagasLowRT))
for(i in 1:length(ListSobolindicesPropionateDietChagasLowRT)){
  
  PropionateDietChagasLow<-PropionateDietChagasLowRT[[i]]
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  ListSobolindicesPropionateDietChagasLowRT[[i]]<-SobolindicesComputation(x=PropionateDietChagasLow,RT=InputsRTi)
  
  rm(PropionateDietChagasLow);rm(InputsRTi)
  
}

#Recovering of Sobol indices and standard deviation
#Estimate
#Full indices
#First-order indice
SfullPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullPropionateDietChagasLow)<-Inputs
colnames(SfullPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Total indice
TfullPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullPropionateDietChagasLow)<-Inputs
colnames(TfullPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindPropionateDietChagasLow)<-Inputs
colnames(SindPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Total indice
TindPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindPropionateDietChagasLow)<-Inputs
colnames(TindPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullPropionateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Sfull",]
  TfullPropionateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$Sobolindices[[FullIndicescomputed]]["Tfull",]
  
  SindPropionateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Sind",]
  TindPropionateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$Sobolindices[[IndIndicescomputed]]["Tind",]
  
}
#Standard error
#Full indices
#First-order indice
SfullStdErrPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SfullStdErrPropionateDietChagasLow)<-Inputs
colnames(SfullStdErrPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Total indice
TfullStdErrPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TfullStdErrPropionateDietChagasLow)<-Inputs
colnames(TfullStdErrPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Independent indices
#First-order indice
SindStdErrPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(SindStdErrPropionateDietChagasLow)<-Inputs
colnames(SindStdErrPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])
#Total indice
TindStdErrPropionateDietChagasLow<-data.frame(matrix(NA,nrow=k,ncol=length(t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])])))
row.names(TindStdErrPropionateDietChagasLow)<-Inputs
colnames(TindStdErrPropionateDietChagasLow)<-colnames(PropionateDietChagasLowRT[[1]])

for(i in 1:length(ListSobolindicesPropionateDietChagasLowRT)){
  
  InputsRTi<-eval(parse(text=paste("InputsRT",i,sep="")))
  FullIndicescomputed<-InputsRTi[1]
  IndIndicescomputed<-InputsRTi[length(InputsRTi)]
  
  SfullStdErrPropionateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Sfull",]
  TfullStdErrPropionateDietChagasLow[FullIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$SobolindicesStdErr[[FullIndicescomputed]]["Tfull",]
  
  SindStdErrPropionateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Sind",]
  TindStdErrPropionateDietChagasLow[IndIndicescomputed,]<-ListSobolindicesPropionateDietChagasLowRT[[i]]$SobolindicesStdErr[[IndIndicescomputed]]["Tind",]
  
}

#Integration of Sobol indices (estimate + standard error) estimated in output files
#Estimate
#Full indices
#First-order indice
write.table(SfullPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Standard error
#Full indices
#First-order indice
write.table(SfullStdErrPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SfullStdErrPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TfullStdErrPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TfullStdErrPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Independent indices
#First-order indice
write.table(SindStdErrPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/SindStdErrPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)
#Total indice
write.table(TindStdErrPropionateDietChagasLow,
            file=paste(getwd(),'/Sobolindices_txt/Propionate/TindStdErrPropionateDietChagasLow.txt',sep=""),
            row.names=TRUE,col.names=TRUE)

stopCluster(cl)