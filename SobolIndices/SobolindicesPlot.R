# -----------------------------------------------------------------------------------------

# Plotting of the dynamic of full and independent Sobol indices for the 3 dietary scenarios

# -----------------------------------------------------------------------------------------

#Packages used
library(ggplot2)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(ggpubr)
library(parallel)
library(deSolve)

# ---- 1. Sampling of input parameters -----

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
load("C:/Users/Paul Blondiaux/Seafile/Ma bibliothèque/PhD/Sensitivity_analysis/Implementation/InvitrocontinuousconditionModel/Dependentinputs/FullandIndependentSobolindices/DietChagas/Meso/MechanisticmodFunctionsSamplingMat3000Sim.RData")
SamplingMatrix<-listRTSamplesInputsPickandfreeze[[1]]

# ---- 2. Dynamic of the full and independent Sobol indices  ----

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
clusterEvalQ(cl,c(library(deSolve),library(gridExtra),library(pracma),library(scales),library(boot)))
#Functions used for solving the model (for sending information to clusters)
clusterExport(cl,c("ode","rumencdynamicDMIRusiteccond","tDMI","yDMI","listRTSamplesInputsPickandfreeze","sobol_indices","RowSimd4"),envir=environment())

# ---- 2.1. Control treatment (fraction of Asparagopsis taxiformis in the feed = 0%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullMethaneDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Methane/SfullMethaneDietChagasControl.txt",sep=""),
                                          header=T)
#Total indices
TfullMethaneDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Methane/TfullMethaneDietChagasControl.txt",sep=""),
                                          header=T)
#Independent indices
#First-order indices
SindMethaneDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Methane/SindMethaneDietChagasControl.txt",sep=""),
                                         header=T)
#Total indices
TindMethaneDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Methane/TindMethaneDietChagasControl.txt",sep=""),
                                         header=T)

#Interpretation
#First step: Comparison of Tfull and Tind

colors<-c("khyd_ndf"="red","khyd_nsc"="yellow","khyd_pro"="black",
          "km_su"="orange","Ks_su"="darkorange3","km_aa"="green",
          "Ks_aa"="darkgreen","km_h2"="cyan","Ks_h2"="darkblue",
          "k_br"="chocolate4",
          "p1"="hotpink","p2"="purple","p3"="azure3",
          "p4"="burlywood3","p5"="darkred","p6"="khaki",
          "p7"="cyan4","p8"="coral1","otherinputs"="darkslategrey")

PlotTfullMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TfullMethaneDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(-0.40,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindMethaneDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(-0.40,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTfullTindMethaneDietChagasControl<-ggpubr::ggarrange(PlotTfullMethaneDietChagasControl,
                                                         PlotTindMethaneDietChagasControl,
                                                         nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindMethaneDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=40,hjust=0.5))

#Removing of time steps with poor estimates (t = 74 h)
#Tfull
TfullMethaneDietChagasControl$t.74h<-rep(0,nrow(TfullMethaneDietChagasControl))
#Tind
TindMethaneDietChagasControl$t.74h<-rep(0,nrow(TindMethaneDietChagasControl))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullMethaneDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TfullMethaneDietChagasControl)){
  
  for(j in 1:ncol(TfullMethaneDietChagasControl)){
    
    if(TfullMethaneDietChagasControl[i,j]<0){
      TfullMethaneDietChagasControl[i,j]<-0
    }
    
  }
}
#Tind
apply(TindMethaneDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TindMethaneDietChagasControl)){
  
  for(j in 1:ncol(TindMethaneDietChagasControl)){
    
    if(TindMethaneDietChagasControl[i,j]<0){
      TindMethaneDietChagasControl[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindMethaneDietChagasControl<-vector("list",length=nrow(TfullMethaneDietChagasControl))
for(i in 1:nrow(TfullMethaneDietChagasControl)){
  
  DiffTfullTindMethaneDietChagasControl[[i]]<-abs(TfullMethaneDietChagasControl[i,]-TindMethaneDietChagasControl[i,])
  
}
names(DiffTfullTindMethaneDietChagasControl)<-row.names(TfullMethaneDietChagasControl)

unlist(lapply(DiffTfullTindMethaneDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTfullTindMethaneDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to methane variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindMethaneDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(SindMethaneDietChagasControl)){
  
  for(j in 1:ncol(SindMethaneDietChagasControl)){
    
    if(SindMethaneDietChagasControl[i,j]<0){
      SindMethaneDietChagasControl[i,j]<-0
    }
    
  }
}

PlotSindMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(SindMethaneDietChagasControl)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindMethaneDietChagasControl)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindSindMethaneDietChagasControl<-ggpubr::ggarrange(PlotSindMethaneDietChagasControl,
                                                        PlotTindMethaneDietChagasControl,
                                                        nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindMethaneDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=80,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindMethaneDietChagasControl<-vector("list",length=nrow(TindMethaneDietChagasControl))
for(i in 1:nrow(TindMethaneDietChagasControl)){
  
  DiffTindSindMethaneDietChagasControl[[i]]<-abs(TindMethaneDietChagasControl[i,]-SindMethaneDietChagasControl[i,])
  
}
names(DiffTindSindMethaneDietChagasControl)<-row.names(TindMethaneDietChagasControl)

unlist(lapply(DiffTindSindMethaneDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTindSindMethaneDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind >= Sind indicating that the input parameters contributes to methane variability via their interaction

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullAcetateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Acetate/SfullAcetateDietChagasControl.txt",sep=""),
                                          header=T)
#Total indices
TfullAcetateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Acetate/TfullAcetateDietChagasControl.txt",sep=""),
                                          header=T)
#Independent indices
#First-order indices
SindAcetateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Acetate/SindAcetateDietChagasControl.txt",sep=""),
                                         header=T)
#Total indices
TindAcetateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Acetate/TindAcetateDietChagasControl.txt",sep=""),
                                         header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullAcetateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TfullAcetateDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindAcetateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindAcetateDietChagasControl)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindAcetateDietChagasControl<-ggpubr::ggarrange(PlotTfullAcetateDietChagasControl,
                                                         PlotTindAcetateDietChagasControl,
                                                         nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindAcetateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullAcetateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TfullAcetateDietChagasControl)){
  
  for(j in 1:ncol(TfullAcetateDietChagasControl)){
    
    if(TfullAcetateDietChagasControl[i,j]<0){
      TfullAcetateDietChagasControl[i,j]<-0
    }
    
  }
}
#Tind
apply(TindAcetateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TindAcetateDietChagasControl)){
  
  for(j in 1:ncol(TindAcetateDietChagasControl)){
    
    if(TindAcetateDietChagasControl[i,j]<0){
      TindAcetateDietChagasControl[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindAcetateDietChagasControl<-vector("list",length=nrow(TfullAcetateDietChagasControl))
for(i in 1:nrow(TfullAcetateDietChagasControl)){
  
  DiffTfullTindAcetateDietChagasControl[[i]]<-abs(TfullAcetateDietChagasControl[i,]-TindAcetateDietChagasControl[i,])
  
}
names(DiffTfullTindAcetateDietChagasControl)<-row.names(TfullAcetateDietChagasControl)

unlist(lapply(DiffTfullTindAcetateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTfullTindAcetateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to acetate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindAcetateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(SindAcetateDietChagasControl)){
  
  for(j in 1:ncol(SindAcetateDietChagasControl)){
    
    if(SindAcetateDietChagasControl[i,j]<0){
      SindAcetateDietChagasControl[i,j]<-0
    }
    
  }
}

PlotSindAcetateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(SindAcetateDietChagasControl)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindAcetateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindAcetateDietChagasControl)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindAcetateDietChagasControl<-ggpubr::ggarrange(PlotSindAcetateDietChagasControl,
                                                        PlotTindAcetateDietChagasControl,
                                                        nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindAcetateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindAcetateDietChagasControl<-vector("list",length=nrow(TindAcetateDietChagasControl))
for(i in 1:nrow(TindAcetateDietChagasControl)){
  
  DiffTindSindAcetateDietChagasControl[[i]]<-abs(TindAcetateDietChagasControl[i,]-SindAcetateDietChagasControl[i,])
  
}
names(DiffTindSindAcetateDietChagasControl)<-row.names(TindAcetateDietChagasControl)

unlist(lapply(DiffTindSindAcetateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTindSindAcetateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to acetate variability via their own variability

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Butyrate/SfullButyrateDietChagasControl.txt",sep=""),
                                           header=T)
#Total indices
TfullButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Butyrate/TfullButyrateDietChagasControl.txt",sep=""),
                                           header=T)
#Independent indices
#First-order indices
SindButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Butyrate/SindButyrateDietChagasControl.txt",sep=""),
                                          header=T)
#Total indices
TindButyrateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Butyrate/TindButyrateDietChagasControl.txt",sep=""),
                                          header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullButyrateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                            t(TfullButyrateDietChagasControl)),
                                           aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.02,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindButyrateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TindButyrateDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindButyrateDietChagasControl<-ggpubr::ggarrange(PlotTfullButyrateDietChagasControl,
                                                          PlotTindButyrateDietChagasControl,
                                                          nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindButyrateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullButyrateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TfullButyrateDietChagasControl)){
  
  for(j in 1:ncol(TfullButyrateDietChagasControl)){
    
    if(TfullButyrateDietChagasControl[i,j]<0){
      TfullButyrateDietChagasControl[i,j]<-0
    }
    
  }
}
#Tind
apply(TindButyrateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TindButyrateDietChagasControl)){
  
  for(j in 1:ncol(TindButyrateDietChagasControl)){
    
    if(TindButyrateDietChagasControl[i,j]<0){
      TindButyrateDietChagasControl[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindButyrateDietChagasControl<-vector("list",length=nrow(TfullButyrateDietChagasControl))
for(i in 1:nrow(TfullButyrateDietChagasControl)){
  
  DiffTfullTindButyrateDietChagasControl[[i]]<-abs(TfullButyrateDietChagasControl[i,]-TindButyrateDietChagasControl[i,])
  
}
names(DiffTfullTindButyrateDietChagasControl)<-row.names(TfullButyrateDietChagasControl)

unlist(lapply(DiffTfullTindButyrateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTfullTindButyrateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to butyrate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindButyrateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(SindButyrateDietChagasControl)){
  
  for(j in 1:ncol(SindButyrateDietChagasControl)){
    
    if(SindButyrateDietChagasControl[i,j]<0){
      SindButyrateDietChagasControl[i,j]<-0
    }
    
  }
}

PlotSindButyrateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(SindButyrateDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindButyrateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TindButyrateDietChagasControl)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindButyrateDietChagasControl<-ggpubr::ggarrange(PlotSindButyrateDietChagasControl,
                                                         PlotTindButyrateDietChagasControl,
                                                         nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindButyrateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindButyrateDietChagasControl<-vector("list",length=nrow(TindButyrateDietChagasControl))
for(i in 1:nrow(TindButyrateDietChagasControl)){
  
  DiffTindSindButyrateDietChagasControl[[i]]<-abs(TindButyrateDietChagasControl[i,]-SindButyrateDietChagasControl[i,])
  
}
names(DiffTindSindButyrateDietChagasControl)<-row.names(TindButyrateDietChagasControl)

unlist(lapply(DiffTindSindButyrateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTindSindButyrateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Butyrate variability via their own variability

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullPropionateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Propionate/SfullPropionateDietChagasControl.txt",sep=""),
                                             header=T)
#Total indices
TfullPropionateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Propionate/TfullPropionateDietChagasControl.txt",sep=""),
                                             header=T)
#Independent indices
#First-order indices
SindPropionateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Propionate/SindPropionateDietChagasControl.txt",sep=""),
                                            header=T)
#Total indices
TindPropionateDietChagasControl<-read.table(paste(getwd(),"/Control/Sobolindices_txt/Propionate/TindPropionateDietChagasControl.txt",sep=""),
                                            header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullPropionateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                              t(TfullPropionateDietChagasControl)),
                                             aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindPropionateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                             t(TindPropionateDietChagasControl)),
                                            aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindPropionateDietChagasControl<-ggpubr::ggarrange(PlotTfullPropionateDietChagasControl,
                                                         PlotTindPropionateDietChagasControl,
                                                         nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindPropionateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullPropionateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TfullPropionateDietChagasControl)){
  
  for(j in 1:ncol(TfullPropionateDietChagasControl)){
    
    if(TfullPropionateDietChagasControl[i,j]<0){
      TfullPropionateDietChagasControl[i,j]<-0
    }
    
  }
}
#Tind
apply(TindPropionateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(TindPropionateDietChagasControl)){
  
  for(j in 1:ncol(TindPropionateDietChagasControl)){
    
    if(TindPropionateDietChagasControl[i,j]<0){
      TindPropionateDietChagasControl[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindPropionateDietChagasControl<-vector("list",length=nrow(TfullPropionateDietChagasControl))
for(i in 1:nrow(TfullPropionateDietChagasControl)){
  
  DiffTfullTindPropionateDietChagasControl[[i]]<-abs(TfullPropionateDietChagasControl[i,]-TindPropionateDietChagasControl[i,])
  
}
names(DiffTfullTindPropionateDietChagasControl)<-row.names(TfullPropionateDietChagasControl)

unlist(lapply(DiffTfullTindPropionateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTfullTindPropionateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to Propionate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindPropionateDietChagasControl,1,function(x){which(x<0)})
for(i in 1:nrow(SindPropionateDietChagasControl)){
  
  for(j in 1:ncol(SindPropionateDietChagasControl)){
    
    if(SindPropionateDietChagasControl[i,j]<0){
      SindPropionateDietChagasControl[i,j]<-0
    }
    
  }
}

PlotSindPropionateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                             t(SindPropionateDietChagasControl)),
                                            aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindPropionateDietChagasControl<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                             t(TindPropionateDietChagasControl)),
                                            aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindPropionateDietChagasControl<-ggpubr::ggarrange(PlotSindPropionateDietChagasControl,
                                                           PlotTindPropionateDietChagasControl,
                                                           nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindPropionateDietChagasControl,
                top=text_grob("Control",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindPropionateDietChagasControl<-vector("list",length=nrow(TindPropionateDietChagasControl))
for(i in 1:nrow(TindPropionateDietChagasControl)){
  
  DiffTindSindPropionateDietChagasControl[[i]]<-abs(TindPropionateDietChagasControl[i,]-SindPropionateDietChagasControl[i,])
  
}
names(DiffTindSindPropionateDietChagasControl)<-row.names(TindPropionateDietChagasControl)

unlist(lapply(DiffTindSindPropionateDietChagasControl,function(x){max(x)}))
unlist(lapply(DiffTindSindPropionateDietChagasControl,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasControl,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasControl,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasControl$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Control",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Propionate variability via their own variability

# ---- 2.2. Low treatment (fraction of Asparagopsis taxiformis in the feed = 0.25%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullMethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Methane/SfullMethaneDietChagasLow.txt",sep=""),
                                      header=T)
#Total indices
TfullMethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Methane/TfullMethaneDietChagasLow.txt",sep=""),
                                      header=T)
#Independent indices
#First-order indices
SindMethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Methane/SindMethaneDietChagasLow.txt",sep=""),
                                     header=T)
#Total indices
TindMethaneDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Methane/TindMethaneDietChagasLow.txt",sep=""),
                                     header=T)

#Interpretation
#First step: Comparison of Tfull and Tind

PlotTfullMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TfullMethaneDietChagasLow)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.15,1.1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(TindMethaneDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.15,1.1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTfullTindMethaneDietChagasLow<-ggpubr::ggarrange(PlotTfullMethaneDietChagasLow,
                                                     PlotTindMethaneDietChagasLow,
                                                     nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindMethaneDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullMethaneDietChagasLow,1,function(x){which(x<0)})
apply(TfullMethaneDietChagasLow,1,function(x){which(x>1)})
for(i in 1:nrow(TfullMethaneDietChagasLow)){
  
  for(j in 1:ncol(TfullMethaneDietChagasLow)){
    
    if(TfullMethaneDietChagasLow[i,j]<0){
      TfullMethaneDietChagasLow[i,j]<-0
    }
    
    if(TfullMethaneDietChagasLow[i,j]>1){
      TfullMethaneDietChagasLow[i,j]<-1
    }
    
  }
}
#Tind
apply(TindMethaneDietChagasLow,1,function(x){which(x<0)})
apply(TindMethaneDietChagasLow,1,function(x){which(x>1)})
for(i in 1:nrow(TindMethaneDietChagasLow)){
  
  for(j in 1:ncol(TindMethaneDietChagasLow)){
    
    if(TindMethaneDietChagasLow[i,j]<0){
      TindMethaneDietChagasLow[i,j]<-0
    }
    
    if(TindMethaneDietChagasLow[i,j]>1){
      TindMethaneDietChagasLow[i,j]<-1
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindMethaneDietChagasLow<-vector("list",length=nrow(TfullMethaneDietChagasLow))
for(i in 1:nrow(TfullMethaneDietChagasLow)){
  
  DiffTfullTindMethaneDietChagasLow[[i]]<-abs(TfullMethaneDietChagasLow[i,]-TindMethaneDietChagasLow[i,])
  
}
names(DiffTfullTindMethaneDietChagasLow)<-row.names(TfullMethaneDietChagasLow)

unlist(lapply(DiffTfullTindMethaneDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTfullTindMethaneDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to methane variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindMethaneDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(SindMethaneDietChagasLow)){
  
  for(j in 1:ncol(SindMethaneDietChagasLow)){
    
    if(SindMethaneDietChagasLow[i,j]<0){
      SindMethaneDietChagasLow[i,j]<-0
    }
    
  }
}

PlotSindMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(SindMethaneDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(TindMethaneDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindSindMethaneDietChagasLow<-ggpubr::ggarrange(PlotSindMethaneDietChagasLow,
                                                    PlotTindMethaneDietChagasLow,
                                                    nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindMethaneDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=80,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindMethaneDietChagasLow<-vector("list",length=nrow(TindMethaneDietChagasLow))
for(i in 1:nrow(TindMethaneDietChagasLow)){
  
  DiffTindSindMethaneDietChagasLow[[i]]<-abs(TindMethaneDietChagasLow[i,]-SindMethaneDietChagasLow[i,])
  
}
names(DiffTindSindMethaneDietChagasLow)<-row.names(TindMethaneDietChagasLow)

unlist(lapply(DiffTindSindMethaneDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTindSindMethaneDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind >= Sind indicating that the input parameters contributes to methane variability via their interaction

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullAcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Acetate/SfullAcetateDietChagasLow.txt",sep=""),
                                      header=T)
#Total indices
TfullAcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Acetate/TfullAcetateDietChagasLow.txt",sep=""),
                                      header=T)
#Independent indices
#First-order indices
SindAcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Acetate/SindAcetateDietChagasLow.txt",sep=""),
                                     header=T)
#Total indices
TindAcetateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Acetate/TindAcetateDietChagasLow.txt",sep=""),
                                     header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullAcetateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TfullAcetateDietChagasLow)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.15,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindAcetateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(TindAcetateDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.15,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindAcetateDietChagasLow<-ggpubr::ggarrange(PlotTfullAcetateDietChagasLow,
                                                     PlotTindAcetateDietChagasLow,
                                                     nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindAcetateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullAcetateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TfullAcetateDietChagasLow)){
  
  for(j in 1:ncol(TfullAcetateDietChagasLow)){
    
    if(TfullAcetateDietChagasLow[i,j]<0){
      TfullAcetateDietChagasLow[i,j]<-0
    }
    
  }
}
#Tind
apply(TindAcetateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TindAcetateDietChagasLow)){
  
  for(j in 1:ncol(TindAcetateDietChagasLow)){
    
    if(TindAcetateDietChagasLow[i,j]<0){
      TindAcetateDietChagasLow[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindAcetateDietChagasLow<-vector("list",length=nrow(TfullAcetateDietChagasLow))
for(i in 1:nrow(TfullAcetateDietChagasLow)){
  
  DiffTfullTindAcetateDietChagasLow[[i]]<-abs(TfullAcetateDietChagasLow[i,]-TindAcetateDietChagasLow[i,])
  
}
names(DiffTfullTindAcetateDietChagasLow)<-row.names(TfullAcetateDietChagasLow)

unlist(lapply(DiffTfullTindAcetateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTfullTindAcetateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to acetate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindAcetateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(SindAcetateDietChagasLow)){
  
  for(j in 1:ncol(SindAcetateDietChagasLow)){
    
    if(SindAcetateDietChagasLow[i,j]<0){
      SindAcetateDietChagasLow[i,j]<-0
    }
    
  }
}

PlotSindAcetateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(SindAcetateDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindAcetateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                      t(TindAcetateDietChagasLow)),
                                     aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindAcetateDietChagasLow<-ggpubr::ggarrange(PlotSindAcetateDietChagasLow,
                                                    PlotTindAcetateDietChagasLow,
                                                    nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindAcetateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindAcetateDietChagasLow<-vector("list",length=nrow(TindAcetateDietChagasLow))
for(i in 1:nrow(TindAcetateDietChagasLow)){
  
  DiffTindSindAcetateDietChagasLow[[i]]<-abs(TindAcetateDietChagasLow[i,]-SindAcetateDietChagasLow[i,])
  
}
names(DiffTindSindAcetateDietChagasLow)<-row.names(TindAcetateDietChagasLow)

unlist(lapply(DiffTindSindAcetateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTindSindAcetateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to acetate variability via their own variability

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Butyrate/SfullButyrateDietChagasLow.txt",sep=""),
                                       header=T)
#Total indices
TfullButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Butyrate/TfullButyrateDietChagasLow.txt",sep=""),
                                       header=T)
#Independent indices
#First-order indices
SindButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Butyrate/SindButyrateDietChagasLow.txt",sep=""),
                                      header=T)
#Total indices
TindButyrateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Butyrate/TindButyrateDietChagasLow.txt",sep=""),
                                      header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullButyrateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(TfullButyrateDietChagasLow)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.03,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindButyrateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TindButyrateDietChagasLow)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.03,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindButyrateDietChagasLow<-ggpubr::ggarrange(PlotTfullButyrateDietChagasLow,
                                                      PlotTindButyrateDietChagasLow,
                                                      nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindButyrateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullButyrateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TfullButyrateDietChagasLow)){
  
  for(j in 1:ncol(TfullButyrateDietChagasLow)){
    
    if(TfullButyrateDietChagasLow[i,j]<0){
      TfullButyrateDietChagasLow[i,j]<-0
    }
    
  }
}
#Tind
apply(TindButyrateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TindButyrateDietChagasLow)){
  
  for(j in 1:ncol(TindButyrateDietChagasLow)){
    
    if(TindButyrateDietChagasLow[i,j]<0){
      TindButyrateDietChagasLow[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindButyrateDietChagasLow<-vector("list",length=nrow(TfullButyrateDietChagasLow))
for(i in 1:nrow(TfullButyrateDietChagasLow)){
  
  DiffTfullTindButyrateDietChagasLow[[i]]<-abs(TfullButyrateDietChagasLow[i,]-TindButyrateDietChagasLow[i,])
  
}
names(DiffTfullTindButyrateDietChagasLow)<-row.names(TfullButyrateDietChagasLow)

unlist(lapply(DiffTfullTindButyrateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTfullTindButyrateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to butyrate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindButyrateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(SindButyrateDietChagasLow)){
  
  for(j in 1:ncol(SindButyrateDietChagasLow)){
    
    if(SindButyrateDietChagasLow[i,j]<0){
      SindButyrateDietChagasLow[i,j]<-0
    }
    
  }
}

PlotSindButyrateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(SindButyrateDietChagasLow)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindButyrateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TindButyrateDietChagasLow)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindButyrateDietChagasLow<-ggpubr::ggarrange(PlotSindButyrateDietChagasLow,
                                                     PlotTindButyrateDietChagasLow,
                                                     nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindButyrateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindButyrateDietChagasLow<-vector("list",length=nrow(TindButyrateDietChagasLow))
for(i in 1:nrow(TindButyrateDietChagasLow)){
  
  DiffTindSindButyrateDietChagasLow[[i]]<-abs(TindButyrateDietChagasLow[i,]-SindButyrateDietChagasLow[i,])
  
}
names(DiffTindSindButyrateDietChagasLow)<-row.names(TindButyrateDietChagasLow)

unlist(lapply(DiffTindSindButyrateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTindSindButyrateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Butyrate variability via their own variability

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullPropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Propionate/SfullPropionateDietChagasLow.txt",sep=""),
                                         header=T)
#Total indices
TfullPropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Propionate/TfullPropionateDietChagasLow.txt",sep=""),
                                         header=T)
#Independent indices
#First-order indices
SindPropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Propionate/SindPropionateDietChagasLow.txt",sep=""),
                                        header=T)
#Total indices
TindPropionateDietChagasLow<-read.table(paste(getwd(),"/Low_AT_treatment/Sobolindices_txt/Propionate/TindPropionateDietChagasLow.txt",sep=""),
                                        header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullPropionateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TfullPropionateDietChagasLow)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindPropionateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                         t(TindPropionateDietChagasLow)),
                                        aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindPropionateDietChagasLow<-ggpubr::ggarrange(PlotTfullPropionateDietChagasLow,
                                                        PlotTindPropionateDietChagasLow,
                                                        nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindPropionateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullPropionateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TfullPropionateDietChagasLow)){
  
  for(j in 1:ncol(TfullPropionateDietChagasLow)){
    
    if(TfullPropionateDietChagasLow[i,j]<0){
      TfullPropionateDietChagasLow[i,j]<-0
    }
    
  }
}
#Tind
apply(TindPropionateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(TindPropionateDietChagasLow)){
  
  for(j in 1:ncol(TindPropionateDietChagasLow)){
    
    if(TindPropionateDietChagasLow[i,j]<0){
      TindPropionateDietChagasLow[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindPropionateDietChagasLow<-vector("list",length=nrow(TfullPropionateDietChagasLow))
for(i in 1:nrow(TfullPropionateDietChagasLow)){
  
  DiffTfullTindPropionateDietChagasLow[[i]]<-abs(TfullPropionateDietChagasLow[i,]-TindPropionateDietChagasLow[i,])
  
}
names(DiffTfullTindPropionateDietChagasLow)<-row.names(TfullPropionateDietChagasLow)

unlist(lapply(DiffTfullTindPropionateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTfullTindPropionateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to Propionate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindPropionateDietChagasLow,1,function(x){which(x<0)})
for(i in 1:nrow(SindPropionateDietChagasLow)){
  
  for(j in 1:ncol(SindPropionateDietChagasLow)){
    
    if(SindPropionateDietChagasLow[i,j]<0){
      SindPropionateDietChagasLow[i,j]<-0
    }
    
  }
}

PlotSindPropionateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                         t(SindPropionateDietChagasLow)),
                                        aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindPropionateDietChagasLow<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                         t(TindPropionateDietChagasLow)),
                                        aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindPropionateDietChagasLow<-ggpubr::ggarrange(PlotSindPropionateDietChagasLow,
                                                       PlotTindPropionateDietChagasLow,
                                                       nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindPropionateDietChagasLow,
                top=text_grob("Low",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindPropionateDietChagasLow<-vector("list",length=nrow(TindPropionateDietChagasLow))
for(i in 1:nrow(TindPropionateDietChagasLow)){
  
  DiffTindSindPropionateDietChagasLow[[i]]<-abs(TindPropionateDietChagasLow[i,]-SindPropionateDietChagasLow[i,])
  
}
names(DiffTindSindPropionateDietChagasLow)<-row.names(TindPropionateDietChagasLow)

unlist(lapply(DiffTindSindPropionateDietChagasLow,function(x){max(x)}))
unlist(lapply(DiffTindSindPropionateDietChagasLow,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasLow,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasLow,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasLow$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="Low",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Propionate variability via their own variability

# ---- 2.3. High treatment (fraction of Asparagopsis taxiformis in the feed = 0.50%) ----

# ---- A. Model output: Dynamic of methane output flow of gas phase (mol/h) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullMethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Methane/SfullMethaneDietChagasHigh.txt",sep=""),
                                       header=T)
#Total indices
TfullMethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Methane/TfullMethaneDietChagasHigh.txt",sep=""),
                                       header=T)
#Independent indices
#First-order indices
SindMethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Methane/SindMethaneDietChagasHigh.txt",sep=""),
                                      header=T)
#Total indices
TindMethaneDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Methane/TindMethaneDietChagasHigh.txt",sep=""),
                                      header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(TfullMethaneDietChagasHigh)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.40,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TindMethaneDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.40,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTfullTindMethaneDietChagasHigh<-ggpubr::ggarrange(PlotTfullMethaneDietChagasHigh,
                                                      PlotTindMethaneDietChagasHigh,
                                                      nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindMethaneDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullMethaneDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TfullMethaneDietChagasHigh)){
  
  for(j in 1:ncol(TfullMethaneDietChagasHigh)){
    
    if(TfullMethaneDietChagasHigh[i,j]<0){
      TfullMethaneDietChagasHigh[i,j]<-0
    }
    
  }
}
#Tind
apply(TindMethaneDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TindMethaneDietChagasHigh)){
  
  for(j in 1:ncol(TindMethaneDietChagasHigh)){
    
    if(TindMethaneDietChagasHigh[i,j]<0){
      TindMethaneDietChagasHigh[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindMethaneDietChagasHigh<-vector("list",length=nrow(TfullMethaneDietChagasHigh))
for(i in 1:nrow(TfullMethaneDietChagasHigh)){
  
  DiffTfullTindMethaneDietChagasHigh[[i]]<-abs(TfullMethaneDietChagasHigh[i,]-TindMethaneDietChagasHigh[i,])
  
}
names(DiffTfullTindMethaneDietChagasHigh)<-row.names(TfullMethaneDietChagasHigh)

unlist(lapply(DiffTfullTindMethaneDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTfullTindMethaneDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindMethaneDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindMethaneDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to methane variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindMethaneDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(SindMethaneDietChagasHigh)){
  
  for(j in 1:ncol(SindMethaneDietChagasHigh)){
    
    if(SindMethaneDietChagasHigh[i,j]<0){
      SindMethaneDietChagasHigh[i,j]<-0
    }
    
  }
}

PlotSindMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(SindMethaneDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindMethaneDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TindMethaneDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=4))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~q[CH[4]~",g,out"]),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=55),
        axis.title.y=element_text(size=55),
        plot.title=element_text(size=80,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=50,colour="black"),
        axis.text.y=element_text(size=50,colour="black"))

PlotTindSindMethaneDietChagasHigh<-ggpubr::ggarrange(PlotSindMethaneDietChagasHigh,
                                                     PlotTindMethaneDietChagasHigh,
                                                     nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindMethaneDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=80,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindMethaneDietChagasHigh<-vector("list",length=nrow(TindMethaneDietChagasHigh))
for(i in 1:nrow(TindMethaneDietChagasHigh)){
  
  DiffTindSindMethaneDietChagasHigh[[i]]<-abs(TindMethaneDietChagasHigh[i,]-SindMethaneDietChagasHigh[i,])
  
}
names(DiffTindSindMethaneDietChagasHigh)<-row.names(TindMethaneDietChagasHigh)

unlist(lapply(DiffTindSindMethaneDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTindSindMethaneDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindMethaneDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindMethaneDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~q[CH[4]~",g,out"]~"(mol/h)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind >= Sind indicating that the input parameters contributes to methane variability via their interaction

# ---- B. Model output: Dynamic of acetate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullAcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Acetate/SfullAcetateDietChagasHigh.txt",sep=""),
                                       header=T)
#Total indices
TfullAcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Acetate/TfullAcetateDietChagasHigh.txt",sep=""),
                                       header=T)
#Independent indices
#First-order indices
SindAcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Acetate/SindAcetateDietChagasHigh.txt",sep=""),
                                      header=T)
#Total indices
TindAcetateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Acetate/TindAcetateDietChagasHigh.txt",sep=""),
                                      header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullAcetateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(TfullAcetateDietChagasHigh)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindAcetateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TindAcetateDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.05,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindAcetateDietChagasHigh<-ggpubr::ggarrange(PlotTfullAcetateDietChagasHigh,
                                                      PlotTindAcetateDietChagasHigh,
                                                      nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindAcetateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullAcetateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TfullAcetateDietChagasHigh)){
  
  for(j in 1:ncol(TfullAcetateDietChagasHigh)){
    
    if(TfullAcetateDietChagasHigh[i,j]<0){
      TfullAcetateDietChagasHigh[i,j]<-0
    }
    
  }
}
#Tind
apply(TindAcetateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TindAcetateDietChagasHigh)){
  
  for(j in 1:ncol(TindAcetateDietChagasHigh)){
    
    if(TindAcetateDietChagasHigh[i,j]<0){
      TindAcetateDietChagasHigh[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindAcetateDietChagasHigh<-vector("list",length=nrow(TfullAcetateDietChagasHigh))
for(i in 1:nrow(TfullAcetateDietChagasHigh)){
  
  DiffTfullTindAcetateDietChagasHigh[[i]]<-abs(TfullAcetateDietChagasHigh[i,]-TindAcetateDietChagasHigh[i,])
  
}
names(DiffTfullTindAcetateDietChagasHigh)<-row.names(TfullAcetateDietChagasHigh)

unlist(lapply(DiffTfullTindAcetateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTfullTindAcetateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindAcetateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindAcetateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to acetate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindAcetateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(SindAcetateDietChagasHigh)){
  
  for(j in 1:ncol(SindAcetateDietChagasHigh)){
    
    if(SindAcetateDietChagasHigh[i,j]<0){
      SindAcetateDietChagasHigh[i,j]<-0
    }
    
  }
}

PlotSindAcetateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(SindAcetateDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindAcetateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                       t(TindAcetateDietChagasHigh)),
                                      aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[ac]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindAcetateDietChagasHigh<-ggpubr::ggarrange(PlotSindAcetateDietChagasHigh,
                                                     PlotTindAcetateDietChagasHigh,
                                                     nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindAcetateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindAcetateDietChagasHigh<-vector("list",length=nrow(TindAcetateDietChagasHigh))
for(i in 1:nrow(TindAcetateDietChagasHigh)){
  
  DiffTindSindAcetateDietChagasHigh[[i]]<-abs(TindAcetateDietChagasHigh[i,]-SindAcetateDietChagasHigh[i,])
  
}
names(DiffTindSindAcetateDietChagasHigh)<-row.names(TindAcetateDietChagasHigh)

unlist(lapply(DiffTindSindAcetateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTindSindAcetateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindAcetateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindAcetateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[ac]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to acetate variability via their own variability

# ---- C. Model output: Dynamic of butyrate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Butyrate/SfullButyrateDietChagasHigh.txt",sep=""),
                                        header=T)
#Total indices
TfullButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Butyrate/TfullButyrateDietChagasHigh.txt",sep=""),
                                        header=T)
#Independent indices
#First-order indices
SindButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Butyrate/SindButyrateDietChagasHigh.txt",sep=""),
                                       header=T)
#Total indices
TindButyrateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Butyrate/TindButyrateDietChagasHigh.txt",sep=""),
                                       header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullButyrateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                         t(TfullButyrateDietChagasHigh)),
                                        aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.04,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindButyrateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(TindButyrateDietChagasHigh)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.04,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindButyrateDietChagasHigh<-ggpubr::ggarrange(PlotTfullButyrateDietChagasHigh,
                                                       PlotTindButyrateDietChagasHigh,
                                                       nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindButyrateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullButyrateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TfullButyrateDietChagasHigh)){
  
  for(j in 1:ncol(TfullButyrateDietChagasHigh)){
    
    if(TfullButyrateDietChagasHigh[i,j]<0){
      TfullButyrateDietChagasHigh[i,j]<-0
    }
    
  }
}
#Tind
apply(TindButyrateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TindButyrateDietChagasHigh)){
  
  for(j in 1:ncol(TindButyrateDietChagasHigh)){
    
    if(TindButyrateDietChagasHigh[i,j]<0){
      TindButyrateDietChagasHigh[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindButyrateDietChagasHigh<-vector("list",length=nrow(TfullButyrateDietChagasHigh))
for(i in 1:nrow(TfullButyrateDietChagasHigh)){
  
  DiffTfullTindButyrateDietChagasHigh[[i]]<-abs(TfullButyrateDietChagasHigh[i,]-TindButyrateDietChagasHigh[i,])
  
}
names(DiffTfullTindButyrateDietChagasHigh)<-row.names(TfullButyrateDietChagasHigh)

unlist(lapply(DiffTfullTindButyrateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTfullTindButyrateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindButyrateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindButyrateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to butyrate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindButyrateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(SindButyrateDietChagasHigh)){
  
  for(j in 1:ncol(SindButyrateDietChagasHigh)){
    
    if(SindButyrateDietChagasHigh[i,j]<0){
      SindButyrateDietChagasHigh[i,j]<-0
    }
    
  }
}

PlotSindButyrateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(SindButyrateDietChagasHigh)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindButyrateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                        t(TindButyrateDietChagasHigh)),
                                       aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[bu]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindButyrateDietChagasHigh<-ggpubr::ggarrange(PlotSindButyrateDietChagasHigh,
                                                      PlotTindButyrateDietChagasHigh,
                                                      nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindButyrateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindButyrateDietChagasHigh<-vector("list",length=nrow(TindButyrateDietChagasHigh))
for(i in 1:nrow(TindButyrateDietChagasHigh)){
  
  DiffTindSindButyrateDietChagasHigh[[i]]<-abs(TindButyrateDietChagasHigh[i,]-SindButyrateDietChagasHigh[i,])
  
}
names(DiffTindSindButyrateDietChagasHigh)<-row.names(TindButyrateDietChagasHigh)

unlist(lapply(DiffTindSindButyrateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTindSindButyrateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindButyrateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindButyrateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[bu]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Butyrate variability via their own variability

# ---- D. Model output: Dynamic of propionate concentration (mol/L) ----

#Recovering of Sobol indices (Sfull, Tfull, Sind and Tind)
#Full indices
#First-order indices
SfullPropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Propionate/SfullPropionateDietChagasHigh.txt",sep=""),
                                          header=T)
#Total indices
TfullPropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Propionate/TfullPropionateDietChagasHigh.txt",sep=""),
                                          header=T)
#Independent indices
#First-order indices
SindPropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Propionate/SindPropionateDietChagasHigh.txt",sep=""),
                                         header=T)
#Total indices
TindPropionateDietChagasHigh<-read.table(paste(getwd(),"/High_AT_treatment/Sobolindices_txt/Propionate/TindPropionateDietChagasHigh.txt",sep=""),
                                         header=T)

#Interpretation
#First step: Comparison of Tfull and Tind
PlotTfullPropionateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                           t(TfullPropionateDietChagasHigh)),
                                          aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.02,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^full~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTindPropionateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindPropionateDietChagasHigh)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(-0.02,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))

PlotTfullTindPropionateDietChagasHigh<-ggpubr::ggarrange(PlotTfullPropionateDietChagasHigh,
                                                         PlotTindPropionateDietChagasHigh,
                                                         nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTfullTindPropionateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=40,hjust=0.5))

#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Tfull
apply(TfullPropionateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TfullPropionateDietChagasHigh)){
  
  for(j in 1:ncol(TfullPropionateDietChagasHigh)){
    
    if(TfullPropionateDietChagasHigh[i,j]<0){
      TfullPropionateDietChagasHigh[i,j]<-0
    }
    
  }
}
#Tind
apply(TindPropionateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(TindPropionateDietChagasHigh)){
  
  for(j in 1:ncol(TindPropionateDietChagasHigh)){
    
    if(TindPropionateDietChagasHigh[i,j]<0){
      TindPropionateDietChagasHigh[i,j]<-0
    }
    
  }
}

#Computation of Tfull - Tind over time
DiffTfullTindPropionateDietChagasHigh<-vector("list",length=nrow(TfullPropionateDietChagasHigh))
for(i in 1:nrow(TfullPropionateDietChagasHigh)){
  
  DiffTfullTindPropionateDietChagasHigh[[i]]<-abs(TfullPropionateDietChagasHigh[i,]-TindPropionateDietChagasHigh[i,])
  
}
names(DiffTfullTindPropionateDietChagasHigh)<-row.names(TfullPropionateDietChagasHigh)

unlist(lapply(DiffTfullTindPropionateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTfullTindPropionateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTfullTindPropionateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tfull - Tind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTfullTindPropionateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^full~"-"~T^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=40,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Tfull indicating that the input parameters contributes to Propionate variability via their own variability and/or their interaction

#Second step: Comparison of Tind and Sind
#Set values <0 at 0 and values >1 at 1 (numerical approximation)
#Sind
apply(SindPropionateDietChagasHigh,1,function(x){which(x<0)})
for(i in 1:nrow(SindPropionateDietChagasHigh)){
  
  for(j in 1:ncol(SindPropionateDietChagasHigh)){
    
    if(SindPropionateDietChagasHigh[i,j]<0){
      SindPropionateDietChagasHigh[i,j]<-0
    }
    
  }
}

PlotSindPropionateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(SindPropionateDietChagasHigh)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~S^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindPropionateDietChagasHigh<-ggplot(cbind.data.frame(t=t[RowSimd4][c(1,seq(0,length(t[RowSimd4]),by=60)[-1])],
                                                          t(TindPropionateDietChagasHigh)),
                                         aes(x=t))+
  scale_x_continuous(breaks=seq(72,96,by=2))+ylim(0,1)+
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
  labs(title="",
       x="Time (h)",y=bquote(~T^ind~" -- "~s[pr]~"(mol/L)"),
       color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+guides(color=FALSE)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))

PlotTindSindPropionateDietChagasHigh<-ggpubr::ggarrange(PlotSindPropionateDietChagasHigh,
                                                        PlotTindPropionateDietChagasHigh,
                                                        nrow=1,ncol=2,common.legend=TRUE)
annotate_figure(PlotTindSindPropionateDietChagasHigh,
                top=text_grob("High",color="black",face="bold",size=70,hjust=0.5))

#Computation of Tind - Sind over time
DiffTindSindPropionateDietChagasHigh<-vector("list",length=nrow(TindPropionateDietChagasHigh))
for(i in 1:nrow(TindPropionateDietChagasHigh)){
  
  DiffTindSindPropionateDietChagasHigh[[i]]<-abs(TindPropionateDietChagasHigh[i,]-SindPropionateDietChagasHigh[i,])
  
}
names(DiffTindSindPropionateDietChagasHigh)<-row.names(TindPropionateDietChagasHigh)

unlist(lapply(DiffTindSindPropionateDietChagasHigh,function(x){max(x)}))
unlist(lapply(DiffTindSindPropionateDietChagasHigh,function(x){mean(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasHigh,function(x){median(as.numeric(x))}))
unlist(lapply(DiffTindSindPropionateDietChagasHigh,function(x){quantile(as.numeric(x),0.75)}))

#Tind - Sind over time
ggplot()+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_ndf),color="khyd_ndf"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_nsc),color="khyd_nsc"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_pro),color="khyd_pro"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$khyd_pro),color="khyd_pro"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_su),color="km_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_su),color="km_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_su),color="Ks_su"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_su),color="Ks_su"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_aa),color="km_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_aa),color="km_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_aa),color="Ks_aa"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_aa),color="Ks_aa"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_h2),color="km_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$km_h2),color="km_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_h2),color="Ks_h2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$Ks_h2),color="Ks_h2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$k_br),color="k_br"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$k_br),color="k_br"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p1),color="p1"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p1),color="p1"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p2),color="p2"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p2),color="p2"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p3),color="p3"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p3),color="p3"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p4),color="p4"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p4),color="p4"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p5),color="p5"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p5),color="p5"),
            linetype="solid",linewidth=2)+
  geom_point(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p6),color="p6"),
             size=5,shape=16)+
  geom_line(aes(x=seq(72,96,by=1),y=as.numeric(DiffTindSindPropionateDietChagasHigh$p6),color="p6"),
            linetype="solid",linewidth=2)+
  scale_x_continuous(breaks=seq(72,96,by=2))+
  labs(title="High",x="Time (h)",y=bquote(~"|"~T^ind~"-"~S^ind~"|"~" -- "~s[pr]~"(mol/L)"),color="Input parameters")+
  scale_color_manual(values=colors[-length(colors)],
                     labels=c(bquote(~k[br]),bquote(~k[hyd~","~ndf]),bquote(~k[hyd~","~nsc]),bquote(~k[hyd~","~pro]),
                              bquote(~k[m~","~aa]),bquote(~k[m~","~H[2]]),bquote(~k[m~","~su]),
                              bquote(~K[S~","~aa]),bquote(~K[S~","~H[2]]),bquote(~K[S~","~su]),
                              bquote(~p[1]),bquote(~p[2]),bquote(~p[3]),
                              bquote(~p[4]),bquote(~p[5]),bquote(~p[6])))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        plot.title=element_text(size=30,hjust=0.5),
        legend.title=element_text(size=30),
        legend.text=element_text(size=25),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))
#Tind ~ Sind indicating that the input parameters contributes to Propionate variability via their own variability

stopCluster(cl) #Stop the clusters
