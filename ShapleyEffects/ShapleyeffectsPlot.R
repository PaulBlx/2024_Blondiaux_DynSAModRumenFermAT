# ---------------------------------------------------------------------------------------------------------

# Plotting of the dynamic of output variable simulations and of Shapley effects for the 3 dietary scenarios

# ---------------------------------------------------------------------------------------------------------

#Packages used
library(deSolve)
library(ggplot2)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(pracma)
library(sensitivity)
library(ggpubr)
library(parallel)

# ---- 1. Input parameters sampling matrix -----

# List of input parameters studied:

# khyd_ndf: Hydrolysis rate constant of neutral detergent fiber, 1/h 
# khyd_nsc: Hydrolysis rate constant of non-structural carbohydrates, 1/h  
# khyd_pro: Hydrolysis rate constant of proteins, 1/h 
# km_su: Maximum specific utilization rate constant of sugars, mol/(mol h) 
# Ks_su: Monod constant associated with the utilization of sugars, mol/L   
# km_aa: Maximum specific utilization rate constant of amino acids, mol/(mol h)
# Ks_aa: Monod constant associated with the utilization of amino acids mol/L 
# km_h2: Maximum specific utilization rate constant of hydrogen, mol/(mol h)
# Ks_h2: Monod constant associated with the utilization of hydrogen, mol/L
# k_br: Kinetic rate constant of bromoform utilization, 1/h
# p1: parameter of the sigmoid function associated with the bromoform inhibition factor (Ibr) 
# p2: parameter of the sigmoid function associated with the bromoform inhibition factor (Ibr)
# p3: parameter of the affine function associated with the flux allocation parameter towards acetate production (lambda1)
# p4: parameter of the affine function associated with the flux allocation parameter towards acetate production (lambda1)
# p5: parameter of the affine function associated with the flux allocation parameter towards propionate production (lambda2)   
# p6: parameter of the affine function associated with the flux allocation parameter towards propionate production (lambda2) 

#Sampling matrix of input parameters
load(paste(getwd(),"/ModRumenFermentationATIPsSamplingMat.RData",sep=""))

# ---- 2. Uncertainty analysis of output variable simulations  ----

# -- General parameters of simulations and parallel computing implementation --

#Time scale considered
nd<-4 #number of days simulated
ts<-1/60 #step time, hour
t<-seq(0,24*nd,by=ts)#time in hour
RowSimd4<-seq(length(seq(0,24*(nd-1),by=ts))+1,length(seq(0,24*nd,by=ts)),by=1) #Lines corresponding to the 4th day of simulation

#Rumen volume in Rusitec condition
V_l<-0.74 #Volume in liquid phase of the rumen, L (Communication Peter Moate)
V_g<-0.06 #Volume in gas phase of the rumen, L (Communication Peter Moate)

#Parallel computing
#Settings
nbnodes<-detectCores()
cl<-makeCluster(nbnodes)
#Library used for solving the model (for sending information to clusters) 
clusterEvalQ(cl,c(library(deSolve),library(gridExtra),library(scales),library(pracma),library(sensitivity)))
#Functions used for solving the model (for sending information to clusters) 
clusterExport(cl,c("ode","rumencdynamicDMIRusiteccond","tDMI","yDMI","tell","RowSimd4"),envir=environment())

# ---- 2.1. Control (fraction of Asparagopsis taxiformis in the feed = 0%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
MethaneDietChagasControl<-read.table(paste(getwd(),"/Control/ModRumFermATOutVar_txt/MethaneDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Description of methane output flow of gas phase (mol/h) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianMethaneDietChagasControl<-parApply(cl=cl,X=MethaneDietChagasControl,MARGIN=2,median)
Quantile0.1MethaneDietChagasControl<-parApply(cl=cl,X=MethaneDietChagasControl,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9MethaneDietChagasControl<-parApply(cl=cl,X=MethaneDietChagasControl,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01MethaneDietChagasControl<-parApply(cl=cl,X=MethaneDietChagasControl,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99MethaneDietChagasControl<-parApply(cl=cl,X=MethaneDietChagasControl,MARGIN=2,FUN=quantile,0.99)

DistributionMethaneDietChagasControl<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianMethaneDietChagasControl,
                                                       Quantile0.1=Quantile0.1MethaneDietChagasControl,
                                                       Quantile0.9=Quantile0.9MethaneDietChagasControl,
                                                       Quantile0.01=Quantile0.01MethaneDietChagasControl,
                                                       Quantile0.99=Quantile0.99MethaneDietChagasControl)

#Methane output flow of gas phase (mol/h) distribution over time (excluding simulations 10% more extreme)
PlotDistributionMethaneDietChagasControlQuantile0.10.9<-ggplot(DistributionMethaneDietChagasControl,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,4e-03,by=1e-03),limits=c(0,4e-03))+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Control",x="Time (h)",y=bquote(~q[CH[4]~",g,out"]~"(mol/h)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionMethaneDietChagasControlQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsMethaneDietChagasControl<-apply(MethaneDietChagasControl,2,
                                        function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                               Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsMethaneDietChagasControl<-as.data.frame(do.call(rbind,SumstatsMethaneDietChagasControl))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasControl$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasControl$CV),linetype="solid",linewidth=1)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"Coefficient of variation -"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=30),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListMethaneDietChagasControltsRemExtSim<-vector("list",length=ncol(MethaneDietChagasControl))
names(ListMethaneDietChagasControltsRemExtSim)<-colnames(MethaneDietChagasControl)
for(i in 1:ncol(MethaneDietChagasControl)){
  
  if(length(which(MethaneDietChagasControl[,i]<DistributionMethaneDietChagasControl$Quantile0.1[i]))==0){
    
    MethaneDietChagasControltsSupQ90<-MethaneDietChagasControl[,i][-which(MethaneDietChagasControl[,i]>DistributionMethaneDietChagasControl$Quantile0.9[i])]
    ListMethaneDietChagasControltsRemExtSim[[i]]<-MethaneDietChagasControltsSupQ90
    
  } else{
    
    MethaneDietChagasControltsInfQ10<-MethaneDietChagasControl[,i][-which(MethaneDietChagasControl[,i]<DistributionMethaneDietChagasControl$Quantile0.1[i])]
    MethaneDietChagasControltsInfQ10SupQ90<-MethaneDietChagasControltsInfQ10[-which(MethaneDietChagasControltsInfQ10>DistributionMethaneDietChagasControl$Quantile0.9[i])]
    ListMethaneDietChagasControltsRemExtSim[[i]]<-MethaneDietChagasControltsInfQ10SupQ90
    
  }
  
}
median(unlist(lapply(ListMethaneDietChagasControltsRemExtSim,function(x){sd(x)/mean(x)})))

ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasControltsRemExtSim,function(x){sd(x)}))),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasControltsRemExtSim,function(x){sd(x)}))),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,8e-04,by=2e-04),limits=c(0,8e-04))+
  labs(title="Control",x="Time (h)",y=bquote(~"SD --"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
AcetateDietChagasControl<-read.table(paste(getwd(),"/Control/ModRumFermATOutVar_txt/AcetateDietChagasControl.txt",sep=""),
                                     header=TRUE)

#Description of acetate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianAcetateDietChagasControl<-parApply(cl=cl,X=AcetateDietChagasControl,MARGIN=2,median)
Quantile0.1AcetateDietChagasControl<-parApply(cl=cl,X=AcetateDietChagasControl,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9AcetateDietChagasControl<-parApply(cl=cl,X=AcetateDietChagasControl,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01AcetateDietChagasControl<-parApply(cl=cl,X=AcetateDietChagasControl,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99AcetateDietChagasControl<-parApply(cl=cl,X=AcetateDietChagasControl,MARGIN=2,FUN=quantile,0.99)

DistributionAcetateDietChagasControl<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianAcetateDietChagasControl,
                                                       Quantile0.1=Quantile0.1AcetateDietChagasControl,
                                                       Quantile0.9=Quantile0.9AcetateDietChagasControl,
                                                       Quantile0.01=Quantile0.01AcetateDietChagasControl,
                                                       Quantile0.99=Quantile0.99AcetateDietChagasControl)

#Acetate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionAcetateDietChagasControlQuantile0.10.9<-ggplot(DistributionAcetateDietChagasControl,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.03,0.08)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Control",x="Time (h)",y=bquote(~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionAcetateDietChagasControlQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsAcetateDietChagasControl<-apply(AcetateDietChagasControl,2,
                                        function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                               Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsAcetateDietChagasControl<-as.data.frame(do.call(rbind,SumstatsAcetateDietChagasControl))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasControl$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasControl$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListAcetateDietChagasControltsRemExtSim<-vector("list",length=ncol(AcetateDietChagasControl))
names(ListAcetateDietChagasControltsRemExtSim)<-colnames(AcetateDietChagasControl)
for(i in 1:ncol(AcetateDietChagasControl)){
  
  if(length(which(AcetateDietChagasControl[,i]<DistributionAcetateDietChagasControl$Quantile0.1[i]))==0){
    
    AcetateDietChagasControltsSupQ90<-AcetateDietChagasControl[,i][-which(AcetateDietChagasControl[,i]>DistributionAcetateDietChagasControl$Quantile0.9[i])]
    ListAcetateDietChagasControltsRemExtSim[[i]]<-AcetateDietChagasControltsSupQ90
    
  } else{
    
    AcetateDietChagasControltsInfQ10<-AcetateDietChagasControl[,i][-which(AcetateDietChagasControl[,i]<DistributionAcetateDietChagasControl$Quantile0.1[i])]
    AcetateDietChagasControltsInfQ10SupQ90<-AcetateDietChagasControltsInfQ10[-which(AcetateDietChagasControltsInfQ10>DistributionAcetateDietChagasControl$Quantile0.9[i])]
    ListAcetateDietChagasControltsRemExtSim[[i]]<-AcetateDietChagasControltsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListAcetateDietChagasControltsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/ModRumFermATOutVar_txt/ButyrateDietChagasControl.txt",sep=""),
                                      header=TRUE)

#Description of butyrate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianButyrateDietChagasControl<-parApply(cl=cl,X=ButyrateDietChagasControl,MARGIN=2,median)
Quantile0.1ButyrateDietChagasControl<-parApply(cl=cl,X=ButyrateDietChagasControl,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9ButyrateDietChagasControl<-parApply(cl=cl,X=ButyrateDietChagasControl,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01ButyrateDietChagasControl<-parApply(cl=cl,X=ButyrateDietChagasControl,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99ButyrateDietChagasControl<-parApply(cl=cl,X=ButyrateDietChagasControl,MARGIN=2,FUN=quantile,0.99)

DistributionButyrateDietChagasControl<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianButyrateDietChagasControl,
                                                        Quantile0.1=Quantile0.1ButyrateDietChagasControl,
                                                        Quantile0.9=Quantile0.9ButyrateDietChagasControl,
                                                        Quantile0.01=Quantile0.01ButyrateDietChagasControl,
                                                        Quantile0.99=Quantile0.99ButyrateDietChagasControl)

#Butyrate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionButyrateDietChagasControlQuantile0.10.9<-ggplot(DistributionButyrateDietChagasControl,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.01,0.025)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Control",x="Time (h)",y=bquote(~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionButyrateDietChagasControlQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsButyrateDietChagasControl<-apply(ButyrateDietChagasControl,2,
                                         function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                                Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsButyrateDietChagasControl<-as.data.frame(do.call(rbind,SumstatsButyrateDietChagasControl))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasControl$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasControl$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListButyrateDietChagasControltsRemExtSim<-vector("list",length=ncol(ButyrateDietChagasControl))
names(ListButyrateDietChagasControltsRemExtSim)<-colnames(ButyrateDietChagasControl)
for(i in 1:ncol(ButyrateDietChagasControl)){
  
  if(length(which(ButyrateDietChagasControl[,i]<DistributionButyrateDietChagasControl$Quantile0.1[i]))==0){
    
    ButyrateDietChagasControltsSupQ90<-ButyrateDietChagasControl[,i][-which(ButyrateDietChagasControl[,i]>DistributionButyrateDietChagasControl$Quantile0.9[i])]
    ListButyrateDietChagasControltsRemExtSim[[i]]<-ButyrateDietChagasControltsSupQ90
    
  } else{
    
    ButyrateDietChagasControltsInfQ10<-ButyrateDietChagasControl[,i][-which(ButyrateDietChagasControl[,i]<DistributionButyrateDietChagasControl$Quantile0.1[i])]
    ButyrateDietChagasControltsInfQ10SupQ90<-ButyrateDietChagasControltsInfQ10[-which(ButyrateDietChagasControltsInfQ10>DistributionButyrateDietChagasControl$Quantile0.9[i])]
    ListButyrateDietChagasControltsRemExtSim[[i]]<-ButyrateDietChagasControltsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListButyrateDietChagasControltsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
PropionateDietChagasControl<-read.table(paste(getwd(),"/Control/ModRumFermATOutVar_txt/PropionateDietChagasControl.txt",sep=""),
                                        header=TRUE)

#Description of propionate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianPropionateDietChagasControl<-parApply(cl=cl,X=PropionateDietChagasControl,MARGIN=2,median)
Quantile0.1PropionateDietChagasControl<-parApply(cl=cl,X=PropionateDietChagasControl,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9PropionateDietChagasControl<-parApply(cl=cl,X=PropionateDietChagasControl,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01PropionateDietChagasControl<-parApply(cl=cl,X=PropionateDietChagasControl,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99PropionateDietChagasControl<-parApply(cl=cl,X=PropionateDietChagasControl,MARGIN=2,FUN=quantile,0.99)

DistributionPropionateDietChagasControl<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianPropionateDietChagasControl,
                                                          Quantile0.1=Quantile0.1PropionateDietChagasControl,
                                                          Quantile0.9=Quantile0.9PropionateDietChagasControl,
                                                          Quantile0.01=Quantile0.01PropionateDietChagasControl,
                                                          Quantile0.99=Quantile0.99PropionateDietChagasControl)

#Propionate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionPropionateDietChagasControlQuantile0.10.9<-ggplot(DistributionPropionateDietChagasControl,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.015,0.04)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Control",x="Time (h)",y=bquote(~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionPropionateDietChagasControlQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsPropionateDietChagasControl<-apply(PropionateDietChagasControl,2,
                                           function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                                  Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsPropionateDietChagasControl<-as.data.frame(do.call(rbind,SumstatsPropionateDietChagasControl))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasControl$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasControl$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListPropionateDietChagasControltsRemExtSim<-vector("list",length=ncol(PropionateDietChagasControl))
names(ListPropionateDietChagasControltsRemExtSim)<-colnames(PropionateDietChagasControl)
for(i in 1:ncol(PropionateDietChagasControl)){
  
  if(length(which(PropionateDietChagasControl[,i]<DistributionPropionateDietChagasControl$Quantile0.1[i]))==0){
    
    PropionateDietChagasControltsSupQ90<-PropionateDietChagasControl[,i][-which(PropionateDietChagasControl[,i]>DistributionPropionateDietChagasControl$Quantile0.9[i])]
    ListPropionateDietChagasControltsRemExtSim[[i]]<-PropionateDietChagasControltsSupQ90
    
  } else{
    
    PropionateDietChagasControltsInfQ10<-PropionateDietChagasControl[,i][-which(PropionateDietChagasControl[,i]<DistributionPropionateDietChagasControl$Quantile0.1[i])]
    PropionateDietChagasControltsInfQ10SupQ90<-PropionateDietChagasControltsInfQ10[-which(PropionateDietChagasControltsInfQ10>DistributionPropionateDietChagasControl$Quantile0.9[i])]
    ListPropionateDietChagasControltsRemExtSim[[i]]<-PropionateDietChagasControltsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListPropionateDietChagasControltsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- 2.2. Low treatment (fraction of Asparagopsis taxiformis in the feed = 0.25%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
MethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ModRumFermATOutVar_txt/MethaneDietChagasLow.txt",sep=""),
                                 header=TRUE)

#Description of methane output flow of gas phase (mol/h) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianMethaneDietChagasLow<-parApply(cl=cl,X=MethaneDietChagasLow,MARGIN=2,median)
Quantile0.1MethaneDietChagasLow<-parApply(cl=cl,X=MethaneDietChagasLow,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9MethaneDietChagasLow<-parApply(cl=cl,X=MethaneDietChagasLow,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01MethaneDietChagasLow<-parApply(cl=cl,X=MethaneDietChagasLow,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99MethaneDietChagasLow<-parApply(cl=cl,X=MethaneDietChagasLow,MARGIN=2,FUN=quantile,0.99)

DistributionMethaneDietChagasLow<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianMethaneDietChagasLow,
                                                   Quantile0.1=Quantile0.1MethaneDietChagasLow,
                                                   Quantile0.9=Quantile0.9MethaneDietChagasLow,
                                                   Quantile0.01=Quantile0.01MethaneDietChagasLow,
                                                   Quantile0.99=Quantile0.99MethaneDietChagasLow)

#Methane output flow of gas phase (mol/h) distribution over time (excluding simulations 10% more extreme)
PlotDistributionMethaneDietChagasLowQuantile0.10.9<-ggplot(DistributionMethaneDietChagasLow,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,4e-03,by=1e-03),limits=c(0,4e-03))+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Low",x="Time (h)",y=bquote(~q[CH[4]~",g,out"]~"(mol/h)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionMethaneDietChagasLowQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsMethaneDietChagasLow<-apply(MethaneDietChagasLow,2,
                                    function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                           Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsMethaneDietChagasLow<-as.data.frame(do.call(rbind,SumstatsMethaneDietChagasLow))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasLow$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasLow$CV),linetype="solid",linewidth=1)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"Coefficient of variation -"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=30),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListMethaneDietChagasLowtsRemExtSim<-vector("list",length=ncol(MethaneDietChagasLow))
names(ListMethaneDietChagasLowtsRemExtSim)<-colnames(MethaneDietChagasLow)
for(i in 1:ncol(MethaneDietChagasLow)){
  
  if(length(which(MethaneDietChagasLow[,i]<DistributionMethaneDietChagasLow$Quantile0.1[i]))==0){
    
    MethaneDietChagasLowtsSupQ90<-MethaneDietChagasLow[,i][-which(MethaneDietChagasLow[,i]>DistributionMethaneDietChagasLow$Quantile0.9[i])]
    ListMethaneDietChagasLowtsRemExtSim[[i]]<-MethaneDietChagasLowtsSupQ90
    
  } else{
    
    MethaneDietChagasLowtsInfQ10<-MethaneDietChagasLow[,i][-which(MethaneDietChagasLow[,i]<DistributionMethaneDietChagasLow$Quantile0.1[i])]
    MethaneDietChagasLowtsInfQ10SupQ90<-MethaneDietChagasLowtsInfQ10[-which(MethaneDietChagasLowtsInfQ10>DistributionMethaneDietChagasLow$Quantile0.9[i])]
    ListMethaneDietChagasLowtsRemExtSim[[i]]<-MethaneDietChagasLowtsInfQ10SupQ90
    
  }
  
}
median(unlist(lapply(ListMethaneDietChagasLowtsRemExtSim,function(x){sd(x)/mean(x)})))

ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasLowtsRemExtSim,function(x){sd(x)}))),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasLowtsRemExtSim,function(x){sd(x)}))),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,8e-04,by=2e-04),limits=c(0,8e-04))+
  labs(title="Low",x="Time (h)",y=bquote(~"SD --"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
AcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ModRumFermATOutVar_txt/AcetateDietChagasLow.txt",sep=""),
                                 header=TRUE)

#Description of acetate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianAcetateDietChagasLow<-parApply(cl=cl,X=AcetateDietChagasLow,MARGIN=2,median)
Quantile0.1AcetateDietChagasLow<-parApply(cl=cl,X=AcetateDietChagasLow,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9AcetateDietChagasLow<-parApply(cl=cl,X=AcetateDietChagasLow,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01AcetateDietChagasLow<-parApply(cl=cl,X=AcetateDietChagasLow,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99AcetateDietChagasLow<-parApply(cl=cl,X=AcetateDietChagasLow,MARGIN=2,FUN=quantile,0.99)

DistributionAcetateDietChagasLow<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianAcetateDietChagasLow,
                                                   Quantile0.1=Quantile0.1AcetateDietChagasLow,
                                                   Quantile0.9=Quantile0.9AcetateDietChagasLow,
                                                   Quantile0.01=Quantile0.01AcetateDietChagasLow,
                                                   Quantile0.99=Quantile0.99AcetateDietChagasLow)

#Acetate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionAcetateDietChagasLowQuantile0.10.9<-ggplot(DistributionAcetateDietChagasLow,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.03,0.07)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Low",x="Time (h)",y=bquote(~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionAcetateDietChagasLowQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsAcetateDietChagasLow<-apply(AcetateDietChagasLow,2,
                                    function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                           Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsAcetateDietChagasLow<-as.data.frame(do.call(rbind,SumstatsAcetateDietChagasLow))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasLow$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasLow$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListAcetateDietChagasLowtsRemExtSim<-vector("list",length=ncol(AcetateDietChagasLow))
names(ListAcetateDietChagasLowtsRemExtSim)<-colnames(AcetateDietChagasLow)
for(i in 1:ncol(AcetateDietChagasLow)){
  
  if(length(which(AcetateDietChagasLow[,i]<DistributionAcetateDietChagasLow$Quantile0.1[i]))==0){
    
    AcetateDietChagasLowtsSupQ90<-AcetateDietChagasLow[,i][-which(AcetateDietChagasLow[,i]>DistributionAcetateDietChagasLow$Quantile0.9[i])]
    ListAcetateDietChagasLowtsRemExtSim[[i]]<-AcetateDietChagasLowtsSupQ90
    
  } else{
    
    AcetateDietChagasLowtsInfQ10<-AcetateDietChagasLow[,i][-which(AcetateDietChagasLow[,i]<DistributionAcetateDietChagasLow$Quantile0.1[i])]
    AcetateDietChagasLowtsInfQ10SupQ90<-AcetateDietChagasLowtsInfQ10[-which(AcetateDietChagasLowtsInfQ10>DistributionAcetateDietChagasLow$Quantile0.9[i])]
    ListAcetateDietChagasLowtsRemExtSim[[i]]<-AcetateDietChagasLowtsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListAcetateDietChagasLowtsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ModRumFermATOutVar_txt/ButyrateDietChagasLow.txt",sep=""),
                                  header=TRUE)

#Description of butyrate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianButyrateDietChagasLow<-parApply(cl=cl,X=ButyrateDietChagasLow,MARGIN=2,median)
Quantile0.1ButyrateDietChagasLow<-parApply(cl=cl,X=ButyrateDietChagasLow,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9ButyrateDietChagasLow<-parApply(cl=cl,X=ButyrateDietChagasLow,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01ButyrateDietChagasLow<-parApply(cl=cl,X=ButyrateDietChagasLow,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99ButyrateDietChagasLow<-parApply(cl=cl,X=ButyrateDietChagasLow,MARGIN=2,FUN=quantile,0.99)

DistributionButyrateDietChagasLow<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianButyrateDietChagasLow,
                                                    Quantile0.1=Quantile0.1ButyrateDietChagasLow,
                                                    Quantile0.9=Quantile0.9ButyrateDietChagasLow,
                                                    Quantile0.01=Quantile0.01ButyrateDietChagasLow,
                                                    Quantile0.99=Quantile0.99ButyrateDietChagasLow)

#Butyrate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionButyrateDietChagasLowQuantile0.10.9<-ggplot(DistributionButyrateDietChagasLow,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.01,0.035)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Low",x="Time (h)",y=bquote(~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionButyrateDietChagasLowQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsButyrateDietChagasLow<-apply(ButyrateDietChagasLow,2,
                                     function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                            Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsButyrateDietChagasLow<-as.data.frame(do.call(rbind,SumstatsButyrateDietChagasLow))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasLow$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasLow$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListButyrateDietChagasLowtsRemExtSim<-vector("list",length=ncol(ButyrateDietChagasLow))
names(ListButyrateDietChagasLowtsRemExtSim)<-colnames(ButyrateDietChagasLow)
for(i in 1:ncol(ButyrateDietChagasLow)){
  
  if(length(which(ButyrateDietChagasLow[,i]<DistributionButyrateDietChagasLow$Quantile0.1[i]))==0){
    
    ButyrateDietChagasLowtsSupQ90<-ButyrateDietChagasLow[,i][-which(ButyrateDietChagasLow[,i]>DistributionButyrateDietChagasLow$Quantile0.9[i])]
    ListButyrateDietChagasLowtsRemExtSim[[i]]<-ButyrateDietChagasLowtsSupQ90
    
  } else{
    
    ButyrateDietChagasLowtsInfQ10<-ButyrateDietChagasLow[,i][-which(ButyrateDietChagasLow[,i]<DistributionButyrateDietChagasLow$Quantile0.1[i])]
    ButyrateDietChagasLowtsInfQ10SupQ90<-ButyrateDietChagasLowtsInfQ10[-which(ButyrateDietChagasLowtsInfQ10>DistributionButyrateDietChagasLow$Quantile0.9[i])]
    ListButyrateDietChagasLowtsRemExtSim[[i]]<-ButyrateDietChagasLowtsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListButyrateDietChagasLowtsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
PropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ModRumFermATOutVar_txt/PropionateDietChagasLow.txt",sep=""),
                                    header=TRUE)

#Description of propionate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianPropionateDietChagasLow<-parApply(cl=cl,X=PropionateDietChagasLow,MARGIN=2,median)
Quantile0.1PropionateDietChagasLow<-parApply(cl=cl,X=PropionateDietChagasLow,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9PropionateDietChagasLow<-parApply(cl=cl,X=PropionateDietChagasLow,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01PropionateDietChagasLow<-parApply(cl=cl,X=PropionateDietChagasLow,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99PropionateDietChagasLow<-parApply(cl=cl,X=PropionateDietChagasLow,MARGIN=2,FUN=quantile,0.99)

DistributionPropionateDietChagasLow<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianPropionateDietChagasLow,
                                                      Quantile0.1=Quantile0.1PropionateDietChagasLow,
                                                      Quantile0.9=Quantile0.9PropionateDietChagasLow,
                                                      Quantile0.01=Quantile0.01PropionateDietChagasLow,
                                                      Quantile0.99=Quantile0.99PropionateDietChagasLow)

#Propionate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionPropionateDietChagasLowQuantile0.10.9<-ggplot(DistributionPropionateDietChagasLow,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.015,0.04)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="Low",x="Time (h)",y=bquote(~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionPropionateDietChagasLowQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsPropionateDietChagasLow<-apply(PropionateDietChagasLow,2,
                                       function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                              Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsPropionateDietChagasLow<-as.data.frame(do.call(rbind,SumstatsPropionateDietChagasLow))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasLow$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasLow$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListPropionateDietChagasLowtsRemExtSim<-vector("list",length=ncol(PropionateDietChagasLow))
names(ListPropionateDietChagasLowtsRemExtSim)<-colnames(PropionateDietChagasLow)
for(i in 1:ncol(PropionateDietChagasLow)){
  
  if(length(which(PropionateDietChagasLow[,i]<DistributionPropionateDietChagasLow$Quantile0.1[i]))==0){
    
    PropionateDietChagasLowtsSupQ90<-PropionateDietChagasLow[,i][-which(PropionateDietChagasLow[,i]>DistributionPropionateDietChagasLow$Quantile0.9[i])]
    ListPropionateDietChagasLowtsRemExtSim[[i]]<-PropionateDietChagasLowtsSupQ90
    
  } else{
    
    PropionateDietChagasLowtsInfQ10<-PropionateDietChagasLow[,i][-which(PropionateDietChagasLow[,i]<DistributionPropionateDietChagasLow$Quantile0.1[i])]
    PropionateDietChagasLowtsInfQ10SupQ90<-PropionateDietChagasLowtsInfQ10[-which(PropionateDietChagasLowtsInfQ10>DistributionPropionateDietChagasLow$Quantile0.9[i])]
    ListPropionateDietChagasLowtsRemExtSim[[i]]<-PropionateDietChagasLowtsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListPropionateDietChagasLowtsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- 2.3. High treatment (fraction of Asparagopsis taxiformis in the feed = 0.50%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of methane output flow of gas phase (mol/h)
MethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ModRumFermATOutVar_txt/MethaneDietChagasHigh.txt",sep=""),
                                  header=TRUE)

#Description of methane output flow of gas phase (mol/h) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianMethaneDietChagasHigh<-parApply(cl=cl,X=MethaneDietChagasHigh,MARGIN=2,median)
Quantile0.1MethaneDietChagasHigh<-parApply(cl=cl,X=MethaneDietChagasHigh,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9MethaneDietChagasHigh<-parApply(cl=cl,X=MethaneDietChagasHigh,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01MethaneDietChagasHigh<-parApply(cl=cl,X=MethaneDietChagasHigh,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99MethaneDietChagasHigh<-parApply(cl=cl,X=MethaneDietChagasHigh,MARGIN=2,FUN=quantile,0.99)

DistributionMethaneDietChagasHigh<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianMethaneDietChagasHigh,
                                                    Quantile0.1=Quantile0.1MethaneDietChagasHigh,
                                                    Quantile0.9=Quantile0.9MethaneDietChagasHigh,
                                                    Quantile0.01=Quantile0.01MethaneDietChagasHigh,
                                                    Quantile0.99=Quantile0.99MethaneDietChagasHigh)

#Methane output flow of gas phase (mol/h) distribution over time (excluding simulations 10% more extreme)
PlotDistributionMethaneDietChagasHighQuantile0.10.9<-ggplot(DistributionMethaneDietChagasHigh,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,4e-03,by=1e-03),limits=c(0,4e-03))+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="High",x="Time (h)",y=bquote(~q[CH[4]~",g,out"]~"(mol/h)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionMethaneDietChagasHighQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsMethaneDietChagasHigh<-apply(MethaneDietChagasHigh,2,
                                     function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                            Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsMethaneDietChagasHigh<-as.data.frame(do.call(rbind,SumstatsMethaneDietChagasHigh))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasHigh$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsMethaneDietChagasHigh$CV),linetype="solid",linewidth=1)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"Coefficient of variation -"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=30),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListMethaneDietChagasHightsRemExtSim<-vector("list",length=ncol(MethaneDietChagasHigh))
names(ListMethaneDietChagasHightsRemExtSim)<-colnames(MethaneDietChagasHigh)
for(i in 1:ncol(MethaneDietChagasHigh)){
  
  if(length(which(MethaneDietChagasHigh[,i]<DistributionMethaneDietChagasHigh$Quantile0.1[i]))==0){
    
    MethaneDietChagasHightsSupQ90<-MethaneDietChagasHigh[,i][-which(MethaneDietChagasHigh[,i]>DistributionMethaneDietChagasHigh$Quantile0.9[i])]
    ListMethaneDietChagasHightsRemExtSim[[i]]<-MethaneDietChagasHightsSupQ90
    
  } else{
    
    MethaneDietChagasHightsInfQ10<-MethaneDietChagasHigh[,i][-which(MethaneDietChagasHigh[,i]<DistributionMethaneDietChagasHigh$Quantile0.1[i])]
    MethaneDietChagasHightsInfQ10SupQ90<-MethaneDietChagasHightsInfQ10[-which(MethaneDietChagasHightsInfQ10>DistributionMethaneDietChagasHigh$Quantile0.9[i])]
    ListMethaneDietChagasHightsRemExtSim[[i]]<-MethaneDietChagasHightsInfQ10SupQ90
    
  }
  
}
median(unlist(lapply(ListMethaneDietChagasHightsRemExtSim,function(x){sd(x)/mean(x)})))
which.max(unlist(lapply(ListMethaneDietChagasHightsRemExtSim,function(x){sd(x)/mean(x)})))

ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasHightsRemExtSim,function(x){sd(x)}))),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=unlist(lapply(ListMethaneDietChagasHightsRemExtSim,function(x){sd(x)}))),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=seq(0,8e-04,by=2e-04),limits=c(0,8e-04))+
  labs(title="High",x="Time (h)",y=bquote(~"SD --"~q[CH[4]~",g,out"]~"simulations"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of acetate concentration (mol/L)
AcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ModRumFermATOutVar_txt/AcetateDietChagasHigh.txt",sep=""),
                                  header=TRUE)

#Description of acetate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianAcetateDietChagasHigh<-parApply(cl=cl,X=AcetateDietChagasHigh,MARGIN=2,median)
Quantile0.1AcetateDietChagasHigh<-parApply(cl=cl,X=AcetateDietChagasHigh,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9AcetateDietChagasHigh<-parApply(cl=cl,X=AcetateDietChagasHigh,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01AcetateDietChagasHigh<-parApply(cl=cl,X=AcetateDietChagasHigh,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99AcetateDietChagasHigh<-parApply(cl=cl,X=AcetateDietChagasHigh,MARGIN=2,FUN=quantile,0.99)

DistributionAcetateDietChagasHigh<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianAcetateDietChagasHigh,
                                                    Quantile0.1=Quantile0.1AcetateDietChagasHigh,
                                                    Quantile0.9=Quantile0.9AcetateDietChagasHigh,
                                                    Quantile0.01=Quantile0.01AcetateDietChagasHigh,
                                                    Quantile0.99=Quantile0.99AcetateDietChagasHigh)

#Acetate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionAcetateDietChagasHighQuantile0.10.9<-ggplot(DistributionAcetateDietChagasHigh,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.02,0.05)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="High",x="Time (h)",y=bquote(~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionAcetateDietChagasHighQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsAcetateDietChagasHigh<-apply(AcetateDietChagasHigh,2,
                                     function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                            Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsAcetateDietChagasHigh<-as.data.frame(do.call(rbind,SumstatsAcetateDietChagasHigh))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasHigh$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsAcetateDietChagasHigh$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[ac]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListAcetateDietChagasHightsRemExtSim<-vector("list",length=ncol(AcetateDietChagasHigh))
names(ListAcetateDietChagasHightsRemExtSim)<-colnames(AcetateDietChagasHigh)
for(i in 1:ncol(AcetateDietChagasHigh)){
  
  if(length(which(AcetateDietChagasHigh[,i]<DistributionAcetateDietChagasHigh$Quantile0.1[i]))==0){
    
    AcetateDietChagasHightsSupQ90<-AcetateDietChagasHigh[,i][-which(AcetateDietChagasHigh[,i]>DistributionAcetateDietChagasHigh$Quantile0.9[i])]
    ListAcetateDietChagasHightsRemExtSim[[i]]<-AcetateDietChagasHightsSupQ90
    
  } else{
    
    AcetateDietChagasHightsInfQ10<-AcetateDietChagasHigh[,i][-which(AcetateDietChagasHigh[,i]<DistributionAcetateDietChagasHigh$Quantile0.1[i])]
    AcetateDietChagasHightsInfQ10SupQ90<-AcetateDietChagasHightsInfQ10[-which(AcetateDietChagasHightsInfQ10>DistributionAcetateDietChagasHigh$Quantile0.9[i])]
    ListAcetateDietChagasHightsRemExtSim[[i]]<-AcetateDietChagasHightsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListAcetateDietChagasHightsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of butyrate concentration (mol/L)
ButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ModRumFermATOutVar_txt/ButyrateDietChagasHigh.txt",sep=""),
                                   header=TRUE)

#Description of butyrate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianButyrateDietChagasHigh<-parApply(cl=cl,X=ButyrateDietChagasHigh,MARGIN=2,median)
Quantile0.1ButyrateDietChagasHigh<-parApply(cl=cl,X=ButyrateDietChagasHigh,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9ButyrateDietChagasHigh<-parApply(cl=cl,X=ButyrateDietChagasHigh,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01ButyrateDietChagasHigh<-parApply(cl=cl,X=ButyrateDietChagasHigh,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99ButyrateDietChagasHigh<-parApply(cl=cl,X=ButyrateDietChagasHigh,MARGIN=2,FUN=quantile,0.99)

DistributionButyrateDietChagasHigh<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianButyrateDietChagasHigh,
                                                     Quantile0.1=Quantile0.1ButyrateDietChagasHigh,
                                                     Quantile0.9=Quantile0.9ButyrateDietChagasHigh,
                                                     Quantile0.01=Quantile0.01ButyrateDietChagasHigh,
                                                     Quantile0.99=Quantile0.99ButyrateDietChagasHigh)

#Butyrate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionButyrateDietChagasHighQuantile0.10.9<-ggplot(DistributionButyrateDietChagasHigh,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.01,0.04)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="High",x="Time (h)",y=bquote(~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionButyrateDietChagasHighQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsButyrateDietChagasHigh<-apply(ButyrateDietChagasHigh,2,
                                      function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                             Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsButyrateDietChagasHigh<-as.data.frame(do.call(rbind,SumstatsButyrateDietChagasHigh))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasHigh$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsButyrateDietChagasHigh$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[bu]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListButyrateDietChagasHightsRemExtSim<-vector("list",length=ncol(ButyrateDietChagasHigh))
names(ListButyrateDietChagasHightsRemExtSim)<-colnames(ButyrateDietChagasHigh)
for(i in 1:ncol(ButyrateDietChagasHigh)){
  
  if(length(which(ButyrateDietChagasHigh[,i]<DistributionButyrateDietChagasHigh$Quantile0.1[i]))==0){
    
    ButyrateDietChagasHightsSupQ90<-ButyrateDietChagasHigh[,i][-which(ButyrateDietChagasHigh[,i]>DistributionButyrateDietChagasHigh$Quantile0.9[i])]
    ListButyrateDietChagasHightsRemExtSim[[i]]<-ButyrateDietChagasHightsSupQ90
    
  } else{
    
    ButyrateDietChagasHightsInfQ10<-ButyrateDietChagasHigh[,i][-which(ButyrateDietChagasHigh[,i]<DistributionButyrateDietChagasHigh$Quantile0.1[i])]
    ButyrateDietChagasHightsInfQ10SupQ90<-ButyrateDietChagasHightsInfQ10[-which(ButyrateDietChagasHightsInfQ10>DistributionButyrateDietChagasHigh$Quantile0.9[i])]
    ListButyrateDietChagasHightsRemExtSim[[i]]<-ButyrateDietChagasHightsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListButyrateDietChagasHightsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of propionate concentration (mol/L)
PropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ModRumFermATOutVar_txt/PropionateDietChagasHigh.txt",sep=""),
                                     header=TRUE)

#Description of propionate concentration (mol/L) distribution over time (computation of median, 0.1 and 0.9 quantiles, and 0.01 and 0.99 quantiles at each time step)
MedianPropionateDietChagasHigh<-parApply(cl=cl,X=PropionateDietChagasHigh,MARGIN=2,median)
Quantile0.1PropionateDietChagasHigh<-parApply(cl=cl,X=PropionateDietChagasHigh,MARGIN=2,FUN=quantile,0.1) 
Quantile0.9PropionateDietChagasHigh<-parApply(cl=cl,X=PropionateDietChagasHigh,MARGIN=2,FUN=quantile,0.9)  
Quantile0.01PropionateDietChagasHigh<-parApply(cl=cl,X=PropionateDietChagasHigh,MARGIN=2,FUN=quantile,0.01)  
Quantile0.99PropionateDietChagasHigh<-parApply(cl=cl,X=PropionateDietChagasHigh,MARGIN=2,FUN=quantile,0.99)

DistributionPropionateDietChagasHigh<-cbind.data.frame(time=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],Median=MedianPropionateDietChagasHigh,
                                                       Quantile0.1=Quantile0.1PropionateDietChagasHigh,
                                                       Quantile0.9=Quantile0.9PropionateDietChagasHigh,
                                                       Quantile0.01=Quantile0.01PropionateDietChagasHigh,
                                                       Quantile0.99=Quantile0.99PropionateDietChagasHigh)

#Propionate concentration (mol/L) distribution over time (excluding simulations 10% more extreme)
PlotDistributionPropionateDietChagasHighQuantile0.10.9<-ggplot(DistributionPropionateDietChagasHigh,aes(x=time))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0.015,0.04)+
  geom_line(aes(y=Median),linetype="solid",linewidth=2)+
  geom_line(aes(y=Quantile0.1),linetype="dashed",linewidth=2)+
  geom_line(aes(y=Quantile0.9),linetype="dashed",linewidth=2)+
  labs(title="High",x="Time (h)",y=bquote(~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotDistributionPropionateDietChagasHighQuantile0.10.9

#Uncertainty analysis of the simulations
#Computation of the summary statistics for each time step
SumstatsPropionateDietChagasHigh<-apply(PropionateDietChagasHigh,2,
                                        function(x){data.frame(Min=min(x),Q10=quantile(x,0.1),Mean=mean(x),Median=median(x),
                                                               Q90=quantile(x,0.9),Max=max(x),SD=sd(x),CV=sd(x)/mean(x))})
SumstatsPropionateDietChagasHigh<-as.data.frame(do.call(rbind,SumstatsPropionateDietChagasHigh))

#Dynamic of the coefficient of variation over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasHigh$CV),size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=SumstatsPropionateDietChagasHigh$CV),linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"Coefficient of variation  -- "~s[pr]~"(mol/L)"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

#Variability of the simulations within the 10 and 90% quantiles
ListPropionateDietChagasHightsRemExtSim<-vector("list",length=ncol(PropionateDietChagasHigh))
names(ListPropionateDietChagasHightsRemExtSim)<-colnames(PropionateDietChagasHigh)
for(i in 1:ncol(PropionateDietChagasHigh)){
  
  if(length(which(PropionateDietChagasHigh[,i]<DistributionPropionateDietChagasHigh$Quantile0.1[i]))==0){
    
    PropionateDietChagasHightsSupQ90<-PropionateDietChagasHigh[,i][-which(PropionateDietChagasHigh[,i]>DistributionPropionateDietChagasHigh$Quantile0.9[i])]
    ListPropionateDietChagasHightsRemExtSim[[i]]<-PropionateDietChagasHightsSupQ90
    
  } else{
    
    PropionateDietChagasHightsInfQ10<-PropionateDietChagasHigh[,i][-which(PropionateDietChagasHigh[,i]<DistributionPropionateDietChagasHigh$Quantile0.1[i])]
    PropionateDietChagasHightsInfQ10SupQ90<-PropionateDietChagasHightsInfQ10[-which(PropionateDietChagasHightsInfQ10>DistributionPropionateDietChagasHigh$Quantile0.9[i])]
    ListPropionateDietChagasHightsRemExtSim[[i]]<-PropionateDietChagasHightsInfQ10SupQ90
    
  }
  
}
max(unlist(lapply(ListPropionateDietChagasHightsRemExtSim,function(x){sd(x)/mean(x)})))

# ---- 3. Dynamic of the Shapley effects  ----

# ---- 3.1. Control (fraction of Asparagopsis taxiformis in the feed = 0%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsMethaneDietChagasControl<-read.table(paste(getwd(),"/Control/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasControl.txt",sep=""),
                                                   header=TRUE)

colors<-c("khyd_ndf"="red","khyd_nsc"="yellow","khyd_pro"="black",
          "km_su"="orange","Ks_su"="darkorange3","km_aa"="green",
          "Ks_aa"="darkgreen","km_h2"="cyan","Ks_h2"="darkblue",
          "k_br"="chocolate4",
          "p1"="hotpink","p2"="purple","p3"="azure3",
          "p4"="burlywood3","p5"="darkred","p6"="khaki",
          "p7"="cyan4","p8"="coral1","otherinputs"="darkslategrey")

PlotShapleyEffectsMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                    t(ShapleyEffectsMethaneDietChagasControl)),
                                                   aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Control",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasControl

#Removing of time steps with poor estimates (t = 74,75 h)
ShapleyEffectsMethaneDietChagasControl[,"t.74h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasControl[,"t.75h"]<-rep(0,k)

PlotShapleyEffectsMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                    t(ShapleyEffectsMethaneDietChagasControl)),
                                                   aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Control",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasControl

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsMethaneDietChagasControlInfluentialInputs<-ShapleyEffectsMethaneDietChagasControl[which(apply(ShapleyEffectsMethaneDietChagasControl,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasControlInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasControlInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1MethaneDietChagasControl<-rep(NA,ncol(ShapleyEffectsMethaneDietChagasControlInfluentialInputs))
names(TimestepSHsup1MethaneDietChagasControl)<-colnames(ShapleyEffectsMethaneDietChagasControlInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsMethaneDietChagasControlInfluentialInputs)){
  
  TimestepSHsup1MethaneDietChagasControl[i]<-identical(which(ShapleyEffectsMethaneDietChagasControlInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsMethaneDietChagasControlInfluentialInputs[which(ShapleyEffectsMethaneDietChagasControlInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsMethaneDietChagasControlInfluentialInputs[,which(TimestepSHsup1MethaneDietChagasControl=="FALSE")]<-0

ShapleyEffectsMethaneDietChagasControlInfluentialInputs<-ShapleyEffectsMethaneDietChagasControlInfluentialInputs[which(apply(ShapleyEffectsMethaneDietChagasControlInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasControlInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasControlInfluentialInputs)))

PlotShapleyEffectsMethaneDietChagasControlInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                     t(ShapleyEffectsMethaneDietChagasControlInfluentialInputs)),
                                                                    aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,0.7,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Control",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","Ks_h2","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~K[S~","~H[2]]),"Other inputs"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsMethaneDietChagasControlInfluentialInputs

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsAcetateDietChagasControl<-read.table(paste(getwd(),"/Control/ShapEff_txt/Acetate/ShapleyEffectsAcetateDietChagasControl.txt",sep=""),
                                                   header=TRUE)

PlotShapleyEffectsAcetateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                    t(ShapleyEffectsAcetateDietChagasControl)),
                                                   aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Control",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsAcetateDietChagasControl

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsAcetateDietChagasControlInfluentialInputs<-ShapleyEffectsAcetateDietChagasControl[which(apply(ShapleyEffectsAcetateDietChagasControl,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasControlInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasControlInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1AcetateDietChagasControl<-rep(NA,ncol(ShapleyEffectsAcetateDietChagasControlInfluentialInputs))
names(TimestepSHsup1AcetateDietChagasControl)<-colnames(ShapleyEffectsAcetateDietChagasControlInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsAcetateDietChagasControlInfluentialInputs)){
  
  TimestepSHsup1AcetateDietChagasControl[i]<-identical(which(ShapleyEffectsAcetateDietChagasControlInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsAcetateDietChagasControlInfluentialInputs[which(ShapleyEffectsAcetateDietChagasControlInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsAcetateDietChagasControlInfluentialInputs[,which(TimestepSHsup1AcetateDietChagasControl=="FALSE")]<-0

ShapleyEffectsAcetateDietChagasControlInfluentialInputs<-ShapleyEffectsAcetateDietChagasControlInfluentialInputs[which(apply(ShapleyEffectsAcetateDietChagasControlInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasControlInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasControlInfluentialInputs)))

PlotShapleyEffectsAcetateDietChagasControlInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                     t(ShapleyEffectsAcetateDietChagasControlInfluentialInputs)),
                                                                    aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,0.7,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Control",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[ac]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","p3","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),"Other inputs",bquote(~p[3])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsAcetateDietChagasControlInfluentialInputs

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasControl.txt",sep=""),
                                                    header=TRUE)

PlotShapleyEffectsButyrateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                     t(ShapleyEffectsButyrateDietChagasControl)),
                                                    aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Control",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsButyrateDietChagasControl

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsButyrateDietChagasControlInfluentialInputs<-ShapleyEffectsButyrateDietChagasControl[which(apply(ShapleyEffectsButyrateDietChagasControl,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasControlInfluentialInputs,
                                                                           OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasControlInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1ButyrateDietChagasControl<-rep(NA,ncol(ShapleyEffectsButyrateDietChagasControlInfluentialInputs))
names(TimestepSHsup1ButyrateDietChagasControl)<-colnames(ShapleyEffectsButyrateDietChagasControlInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsButyrateDietChagasControlInfluentialInputs)){
  
  TimestepSHsup1ButyrateDietChagasControl[i]<-identical(which(ShapleyEffectsButyrateDietChagasControlInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsButyrateDietChagasControlInfluentialInputs[which(ShapleyEffectsButyrateDietChagasControlInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsButyrateDietChagasControlInfluentialInputs[,which(TimestepSHsup1ButyrateDietChagasControl=="FALSE")]<-0

ShapleyEffectsButyrateDietChagasControlInfluentialInputs<-ShapleyEffectsButyrateDietChagasControlInfluentialInputs[which(apply(ShapleyEffectsButyrateDietChagasControlInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasControlInfluentialInputs,
                                                                           OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasControlInfluentialInputs)))

PlotShapleyEffectsButyrateDietChagasControlInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                      t(ShapleyEffectsButyrateDietChagasControlInfluentialInputs)),
                                                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Control",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[bu]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","p3","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),"Other inputs",bquote(~p[3]),bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsButyrateDietChagasControlInfluentialInputs

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsPropionateDietChagasControl<-read.table(paste(getwd(),"/Control/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasControl.txt",sep=""),
                                                      header=TRUE)

PlotShapleyEffectsPropionateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                       t(ShapleyEffectsPropionateDietChagasControl)),
                                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Control",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsPropionateDietChagasControl

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsPropionateDietChagasControlInfluentialInputs<-ShapleyEffectsPropionateDietChagasControl[which(apply(ShapleyEffectsPropionateDietChagasControl,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasControlInfluentialInputs,
                                                                             OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasControlInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1PropionateDietChagasControl<-rep(NA,ncol(ShapleyEffectsPropionateDietChagasControlInfluentialInputs))
names(TimestepSHsup1PropionateDietChagasControl)<-colnames(ShapleyEffectsPropionateDietChagasControlInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsPropionateDietChagasControlInfluentialInputs)){
  
  TimestepSHsup1PropionateDietChagasControl[i]<-identical(which(ShapleyEffectsPropionateDietChagasControlInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsPropionateDietChagasControlInfluentialInputs[which(ShapleyEffectsPropionateDietChagasControlInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsPropionateDietChagasControlInfluentialInputs[,which(TimestepSHsup1PropionateDietChagasControl=="FALSE")]<-0

ShapleyEffectsPropionateDietChagasControlInfluentialInputs<-ShapleyEffectsPropionateDietChagasControlInfluentialInputs[which(apply(ShapleyEffectsPropionateDietChagasControlInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasControlInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasControlInfluentialInputs,
                                                                             OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasControlInfluentialInputs)))

PlotShapleyEffectsPropionateDietChagasControlInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                        t(ShapleyEffectsPropionateDietChagasControlInfluentialInputs)),
                                                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Control",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[pr]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),"Other inputs",bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsPropionateDietChagasControlInfluentialInputs

# ---- 3.2. Low treatment (fraction of Asparagopsis taxiformis in the feed = 0.25%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsMethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasLow.txt",sep=""),
                                               header=TRUE)

PlotShapleyEffectsMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                t(ShapleyEffectsMethaneDietChagasLow)),
                                               aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=1))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Low",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasLow

#Removing of time steps with poor estimates (t = 73,74,75,76,79 h)
ShapleyEffectsMethaneDietChagasLow[,"t.73h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasLow[,"t.74h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasLow[,"t.75h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasLow[,"t.76h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasLow[,"t.79h"]<-rep(0,k)

PlotShapleyEffectsMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                t(ShapleyEffectsMethaneDietChagasLow)),
                                               aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=1))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Low",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasLow

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsMethaneDietChagasLowInfluentialInputs<-ShapleyEffectsMethaneDietChagasLow[which(apply(ShapleyEffectsMethaneDietChagasLow,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasLowInfluentialInputs,
                                                                      OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasLowInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1MethaneDietChagasLow<-rep(NA,ncol(ShapleyEffectsMethaneDietChagasLowInfluentialInputs))
names(TimestepSHsup1MethaneDietChagasLow)<-colnames(ShapleyEffectsMethaneDietChagasLowInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsMethaneDietChagasLowInfluentialInputs)){
  
  TimestepSHsup1MethaneDietChagasLow[i]<-identical(which(ShapleyEffectsMethaneDietChagasLowInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsMethaneDietChagasLowInfluentialInputs[which(ShapleyEffectsMethaneDietChagasLowInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsMethaneDietChagasLowInfluentialInputs[,which(TimestepSHsup1MethaneDietChagasLow=="FALSE")]<-0

ShapleyEffectsMethaneDietChagasLowInfluentialInputs<-ShapleyEffectsMethaneDietChagasLowInfluentialInputs[which(apply(ShapleyEffectsMethaneDietChagasLowInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasLowInfluentialInputs,
                                                                      OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasLowInfluentialInputs)))

PlotShapleyEffectsMethaneDietChagasLowInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                 t(ShapleyEffectsMethaneDietChagasLowInfluentialInputs)),
                                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Low",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","Ks_su","km_h2","Ks_h2","k_br","p2","otherinputs")],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[m~","~H[2]]),
                              bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),"Other inputs",bquote(~p[2])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsMethaneDietChagasLowInfluentialInputs

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsAcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Acetate/ShapEff_txt/ShapleyEffectsAcetateDietChagasLow.txt",sep=""),
                                               header=TRUE)

PlotShapleyEffectsAcetateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                t(ShapleyEffectsAcetateDietChagasLow)),
                                               aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Low",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsAcetateDietChagasLow

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsAcetateDietChagasLowInfluentialInputs<-ShapleyEffectsAcetateDietChagasLow[which(apply(ShapleyEffectsAcetateDietChagasLow,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasLowInfluentialInputs,
                                                                      OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasLowInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1AcetateDietChagasLow<-rep(NA,ncol(ShapleyEffectsAcetateDietChagasLowInfluentialInputs))
names(TimestepSHsup1AcetateDietChagasLow)<-colnames(ShapleyEffectsAcetateDietChagasLowInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsAcetateDietChagasLowInfluentialInputs)){
  
  TimestepSHsup1AcetateDietChagasLow[i]<-identical(which(ShapleyEffectsAcetateDietChagasLowInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsAcetateDietChagasLowInfluentialInputs[which(ShapleyEffectsAcetateDietChagasLowInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsAcetateDietChagasLowInfluentialInputs[,which(TimestepSHsup1AcetateDietChagasLow=="FALSE")]<-0

ShapleyEffectsAcetateDietChagasLowInfluentialInputs<-ShapleyEffectsAcetateDietChagasLowInfluentialInputs[which(apply(ShapleyEffectsAcetateDietChagasLowInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasLowInfluentialInputs,
                                                                      OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasLowInfluentialInputs)))

PlotShapleyEffectsAcetateDietChagasLowInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                 t(ShapleyEffectsAcetateDietChagasLowInfluentialInputs)),
                                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Low",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[ac]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","km_h2","Ks_h2","p3","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[m~","~H[2]]),bquote(~K[S~","~H[2]]),"Other inputs",bquote(~p[3])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsAcetateDietChagasLowInfluentialInputs

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasLow.txt",sep=""),
                                                header=TRUE)

PlotShapleyEffectsButyrateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                 t(ShapleyEffectsButyrateDietChagasLow)),
                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Low",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsButyrateDietChagasLow

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsButyrateDietChagasLowInfluentialInputs<-ShapleyEffectsButyrateDietChagasLow[which(apply(ShapleyEffectsButyrateDietChagasLow,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasLowInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasLowInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1ButyrateDietChagasLow<-rep(NA,ncol(ShapleyEffectsButyrateDietChagasLowInfluentialInputs))
names(TimestepSHsup1ButyrateDietChagasLow)<-colnames(ShapleyEffectsButyrateDietChagasLowInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsButyrateDietChagasLowInfluentialInputs)){
  
  TimestepSHsup1ButyrateDietChagasLow[i]<-identical(which(ShapleyEffectsButyrateDietChagasLowInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsButyrateDietChagasLowInfluentialInputs[which(ShapleyEffectsButyrateDietChagasLowInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsButyrateDietChagasLowInfluentialInputs[,which(TimestepSHsup1ButyrateDietChagasLow=="FALSE")]<-0

ShapleyEffectsButyrateDietChagasLowInfluentialInputs<-ShapleyEffectsButyrateDietChagasLowInfluentialInputs[which(apply(ShapleyEffectsButyrateDietChagasLowInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasLowInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasLowInfluentialInputs)))

PlotShapleyEffectsButyrateDietChagasLowInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                  t(ShapleyEffectsButyrateDietChagasLowInfluentialInputs)),
                                                                 aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Low",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[bu]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","km_h2","Ks_h2","p3","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[m~","~H[2]]),bquote(~K[S~","~H[2]]),
                              "Other inputs",bquote(~p[3]),bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsButyrateDietChagasLowInfluentialInputs

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsPropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasLow.txt",sep=""),
                                                  header=TRUE)

PlotShapleyEffectsPropionateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                   t(ShapleyEffectsPropionateDietChagasLow)),
                                                  aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="Low",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsPropionateDietChagasLow

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsPropionateDietChagasLowInfluentialInputs<-ShapleyEffectsPropionateDietChagasLow[which(apply(ShapleyEffectsPropionateDietChagasLow,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasLowInfluentialInputs,
                                                                         OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasLowInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1PropionateDietChagasLow<-rep(NA,ncol(ShapleyEffectsPropionateDietChagasLowInfluentialInputs))
names(TimestepSHsup1PropionateDietChagasLow)<-colnames(ShapleyEffectsPropionateDietChagasLowInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsPropionateDietChagasLowInfluentialInputs)){
  
  TimestepSHsup1PropionateDietChagasLow[i]<-identical(which(ShapleyEffectsPropionateDietChagasLowInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsPropionateDietChagasLowInfluentialInputs[which(ShapleyEffectsPropionateDietChagasLowInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsPropionateDietChagasLowInfluentialInputs[,which(TimestepSHsup1PropionateDietChagasLow=="FALSE")]<-0

ShapleyEffectsPropionateDietChagasLowInfluentialInputs<-ShapleyEffectsPropionateDietChagasLowInfluentialInputs[which(apply(ShapleyEffectsPropionateDietChagasLowInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasLowInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasLowInfluentialInputs,
                                                                         OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasLowInfluentialInputs)))

PlotShapleyEffectsPropionateDietChagasLowInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                    t(ShapleyEffectsPropionateDietChagasLowInfluentialInputs)),
                                                                   aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="Low",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[pr]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),"Other inputs",bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsPropionateDietChagasLowInfluentialInputs

# ---- 3.3. High treatment (fraction of Asparagopsis taxiformis in the feed = 0.50%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsMethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ShapEff_txt/Methane/ShapleyEffectsMethaneDietChagasHigh.txt",sep=""),
                                                header=TRUE)

PlotShapleyEffectsMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                 t(ShapleyEffectsMethaneDietChagasHigh)),
                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=1))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="High",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasHigh

#Removing of time steps with poor estimates (t = 72,77,78,79,80,84,86 h)
ShapleyEffectsMethaneDietChagasHigh[,"t.72.017h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.77h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.78h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.79h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.80h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.84h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.85h"]<-rep(0,k)
ShapleyEffectsMethaneDietChagasHigh[,"t.86h"]<-rep(0,k)

PlotShapleyEffectsMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                 t(ShapleyEffectsMethaneDietChagasHigh)),
                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=1))+scale_y_continuous(breaks=round(seq(-0.4,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="High",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsMethaneDietChagasHigh

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsMethaneDietChagasHighInfluentialInputs<-ShapleyEffectsMethaneDietChagasHigh[which(apply(ShapleyEffectsMethaneDietChagasHigh,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasHighInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasHighInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1MethaneDietChagasHigh<-rep(NA,ncol(ShapleyEffectsMethaneDietChagasHighInfluentialInputs))
names(TimestepSHsup1MethaneDietChagasHigh)<-colnames(ShapleyEffectsMethaneDietChagasHighInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsMethaneDietChagasHighInfluentialInputs)){
  
  TimestepSHsup1MethaneDietChagasHigh[i]<-identical(which(ShapleyEffectsMethaneDietChagasHighInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsMethaneDietChagasHighInfluentialInputs[which(ShapleyEffectsMethaneDietChagasHighInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsMethaneDietChagasHighInfluentialInputs[,which(TimestepSHsup1MethaneDietChagasHigh=="FALSE")]<-0

ShapleyEffectsMethaneDietChagasHighInfluentialInputs<-ShapleyEffectsMethaneDietChagasHighInfluentialInputs[which(apply(ShapleyEffectsMethaneDietChagasHighInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsMethaneDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsMethaneDietChagasHighInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsMethaneDietChagasHighInfluentialInputs)))

PlotShapleyEffectsMethaneDietChagasHighInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                  t(ShapleyEffectsMethaneDietChagasHighInfluentialInputs)),
                                                                 aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="High",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","khyd_pro","Ks_su","km_h2",
                                     "Ks_h2","k_br","p1","p2","p4","p5","p6","otherinputs")],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~H[2]]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),"Other inputs",
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsMethaneDietChagasHighInfluentialInputs

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsAcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ShapEff_txt/Acetate/ShapleyEffectsAcetateDietChagasHigh.txt",sep=""),
                                                header=TRUE)

PlotShapleyEffectsAcetateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                 t(ShapleyEffectsAcetateDietChagasHigh)),
                                                aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="High",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsAcetateDietChagasHigh

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsAcetateDietChagasHighInfluentialInputs<-ShapleyEffectsAcetateDietChagasHigh[which(apply(ShapleyEffectsAcetateDietChagasHigh,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasHighInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasHighInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1AcetateDietChagasHigh<-rep(NA,ncol(ShapleyEffectsAcetateDietChagasHighInfluentialInputs))
names(TimestepSHsup1AcetateDietChagasHigh)<-colnames(ShapleyEffectsAcetateDietChagasHighInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsAcetateDietChagasHighInfluentialInputs)){
  
  TimestepSHsup1AcetateDietChagasHigh[i]<-identical(which(ShapleyEffectsAcetateDietChagasHighInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsAcetateDietChagasHighInfluentialInputs[which(ShapleyEffectsAcetateDietChagasHighInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsAcetateDietChagasHighInfluentialInputs[,which(TimestepSHsup1AcetateDietChagasHigh=="FALSE")]<-0

ShapleyEffectsAcetateDietChagasHighInfluentialInputs<-ShapleyEffectsAcetateDietChagasHighInfluentialInputs[which(apply(ShapleyEffectsAcetateDietChagasHighInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsAcetateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsAcetateDietChagasHighInfluentialInputs,
                                                                       OtherInputs=rep(0,ncol(ShapleyEffectsAcetateDietChagasHighInfluentialInputs)))

PlotShapleyEffectsAcetateDietChagasHighInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                  t(ShapleyEffectsAcetateDietChagasHighInfluentialInputs)),
                                                                 aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="High",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[ac]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","p3","p4","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),"Other inputs",bquote(~p[3]),bquote(~p[4])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsAcetateDietChagasHighInfluentialInputs

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ShapEff_txt/Butyrate/ShapleyEffectsButyrateDietChagasHigh.txt",sep=""),
                                                 header=TRUE)

PlotShapleyEffectsButyrateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                  t(ShapleyEffectsButyrateDietChagasHigh)),
                                                 aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="High",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsButyrateDietChagasHigh

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsButyrateDietChagasHighInfluentialInputs<-ShapleyEffectsButyrateDietChagasHigh[which(apply(ShapleyEffectsButyrateDietChagasHigh,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasHighInfluentialInputs,
                                                                        OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasHighInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1ButyrateDietChagasHigh<-rep(NA,ncol(ShapleyEffectsButyrateDietChagasHighInfluentialInputs))
names(TimestepSHsup1ButyrateDietChagasHigh)<-colnames(ShapleyEffectsButyrateDietChagasHighInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsButyrateDietChagasHighInfluentialInputs)){
  
  TimestepSHsup1ButyrateDietChagasHigh[i]<-identical(which(ShapleyEffectsButyrateDietChagasHighInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsButyrateDietChagasHighInfluentialInputs[which(ShapleyEffectsButyrateDietChagasHighInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsButyrateDietChagasHighInfluentialInputs[,which(TimestepSHsup1ButyrateDietChagasHigh=="FALSE")]<-0

ShapleyEffectsButyrateDietChagasHighInfluentialInputs<-ShapleyEffectsButyrateDietChagasHighInfluentialInputs[which(apply(ShapleyEffectsButyrateDietChagasHighInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsButyrateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsButyrateDietChagasHighInfluentialInputs,
                                                                        OtherInputs=rep(0,ncol(ShapleyEffectsButyrateDietChagasHighInfluentialInputs)))

PlotShapleyEffectsButyrateDietChagasHighInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                   t(ShapleyEffectsButyrateDietChagasHighInfluentialInputs)),
                                                                  aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="High",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[bu]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","p3","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),"Other inputs",bquote(~p[3]),bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsButyrateDietChagasHighInfluentialInputs

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Estimation of Shapley effects using method of Song et al., 2016 (Using formula of Castro et al., 2009)
#Estimation
ShapleyEffectsPropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/ShapEff_txt/Propionate/ShapleyEffectsPropionateDietChagasHigh.txt",sep=""),
                                                   header=TRUE)

PlotShapleyEffectsPropionateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                    t(ShapleyEffectsPropionateDietChagasHigh)),
                                                   aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+scale_y_continuous(breaks=round(seq(-0.1,0.8,by=0.1),1))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=khyd_pro,color="khyd_pro"),size=5,shape=16)+
  geom_line(aes(y=khyd_pro,color="khyd_pro"),linewidth=2)+
  geom_point(aes(y=km_su,color="km_su"),size=5,shape=16)+
  geom_line(aes(y=km_su,color="km_su"),linewidth=2)+
  geom_point(aes(y=Ks_su,color="Ks_su"),size=5,shape=16)+
  geom_line(aes(y=Ks_su,color="Ks_su"),linewidth=2)+
  geom_point(aes(y=km_aa,color="km_aa"),size=5,shape=16)+
  geom_line(aes(y=km_aa,color="km_aa"),linewidth=2)+
  geom_point(aes(y=Ks_aa,color="Ks_aa"),size=5,shape=16)+
  geom_line(aes(y=Ks_aa,color="Ks_aa"),linewidth=2)+
  geom_point(aes(y=km_h2,color="km_h2"),size=5,shape=16)+
  geom_line(aes(y=km_h2,color="km_h2"),linewidth=2)+
  geom_point(aes(y=Ks_h2,color="Ks_h2"),size=5,shape=16)+
  geom_line(aes(y=Ks_h2,color="Ks_h2"),linewidth=2)+
  geom_point(aes(y=k_br,color="k_br"),size=5,shape=16)+
  geom_line(aes(y=k_br,color="k_br"),linewidth=2)+
  geom_point(aes(y=p1,color="p1"),size=5,shape=16)+
  geom_line(aes(y=p1,color="p1"),linewidth=2)+
  geom_point(aes(y=p2,color="p2"),size=5,shape=16)+
  geom_line(aes(y=p2,color="p2"),linewidth=2)+
  geom_point(aes(y=p3,color="p3"),size=5,shape=16)+
  geom_line(aes(y=p3,color="p3"),linewidth=2)+
  geom_point(aes(y=p4,color="p4"),size=5,shape=16)+
  geom_line(aes(y=p4,color="p4"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_point(aes(y=p6,color="p6"),size=5,shape=16)+
  geom_line(aes(y=p6,color="p6"),linewidth=2)+
  labs(title="High",
       x="Time (h)",y="Shapley effects",
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=70,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=32,colour="black"),
        axis.text.y=element_text(size=32,colour="black"))
PlotShapleyEffectsPropionateDietChagasHigh

#Selection of influential input parameters (at least one time step with a contribution >= 0.10)
ShapleyEffectsPropionateDietChagasHighInfluentialInputs<-ShapleyEffectsPropionateDietChagasHigh[which(apply(ShapleyEffectsPropionateDietChagasHigh,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasHighInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasHighInfluentialInputs)))

#Set values <0 at 0 and values >1 at 1 (numerical approcimation)
TimestepSHsup1PropionateDietChagasHigh<-rep(NA,ncol(ShapleyEffectsPropionateDietChagasHighInfluentialInputs))
names(TimestepSHsup1PropionateDietChagasHigh)<-colnames(ShapleyEffectsPropionateDietChagasHighInfluentialInputs)
for(i in 1:ncol(ShapleyEffectsPropionateDietChagasHighInfluentialInputs)){
  
  TimestepSHsup1PropionateDietChagasHigh[i]<-identical(which(ShapleyEffectsPropionateDietChagasHighInfluentialInputs[,i]>1.0),integer(0))
  ShapleyEffectsPropionateDietChagasHighInfluentialInputs[which(ShapleyEffectsPropionateDietChagasHighInfluentialInputs[,i]<0),i]<-0
  
}
ShapleyEffectsPropionateDietChagasHighInfluentialInputs[,which(TimestepSHsup1PropionateDietChagasHigh=="FALSE")]<-0

ShapleyEffectsPropionateDietChagasHighInfluentialInputs<-ShapleyEffectsPropionateDietChagasHighInfluentialInputs[which(apply(ShapleyEffectsPropionateDietChagasHighInfluentialInputs,1,max)>=0.10),]
ShapleyEffectsPropionateDietChagasHighInfluentialInputs<-rbind.data.frame(ShapleyEffectsPropionateDietChagasHighInfluentialInputs,
                                                                          OtherInputs=rep(0,ncol(ShapleyEffectsPropionateDietChagasHighInfluentialInputs)))

PlotShapleyEffectsPropionateDietChagasHighInfluentialInputs<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                                                     t(ShapleyEffectsPropionateDietChagasHighInfluentialInputs)),
                                                                    aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=3))+scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.75))+
  geom_point(aes(y=khyd_ndf,color="khyd_ndf"),size=5,shape=16)+
  geom_line(aes(y=khyd_ndf,color="khyd_ndf"),linewidth=2)+
  geom_point(aes(y=khyd_nsc,color="khyd_nsc"),size=5,shape=16)+
  geom_line(aes(y=khyd_nsc,color="khyd_nsc"),linewidth=2)+
  geom_point(aes(y=p5,color="p5"),size=5,shape=16)+
  geom_line(aes(y=p5,color="p5"),linewidth=2)+
  geom_line(aes(y=OtherInputs,color="otherinputs"),linewidth=2,linetype="dashed")+
  labs(title="High",
       x="Time (h)",y=bquote(~"Shapley effects  -- "~s[pr]),
       color="Input parameters")+
  scale_color_manual(values=colors[c("khyd_ndf","khyd_nsc","p5","otherinputs")],
                     labels=c(bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),"Other inputs",bquote(~p[5])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))
PlotShapleyEffectsPropionateDietChagasHighInfluentialInputs

stopCluster(cl) #Stop the clusters
