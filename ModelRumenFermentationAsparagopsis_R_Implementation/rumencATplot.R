# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script for running the dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(deSolve)
library(ggplot2)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(pracma)

source(paste(getwd(),"/rumencATload.R",sep=""))
source(paste(getwd(),"/rumencATout.R",sep=""))
source(paste(getwd(),"/rumencAT.R",sep=""))

# ------ Setting input parameters and characteristics of the system studied  ------

#Input parameters considered for the sensitivity analysis
khyd_ndf<-0.024 #Hydrolysis rate constant of cell wall carbohydrates, 1/h 
khyd_nsc<-0.06 #Hydrolysis rate constant of non-structural carbohydrates, 1/h  
khyd_pro<-0.09 #Hydrolysis rate constant of proteins, 1/h 
km_su<-1.0 #Maximum specific utilization rate constant of sugars, mol/(mol h) 
Ks_su<-9e-3 #Monod constant associated with sugar utilization, mol/L   
km_aa<-2.00 #Maximum specific utilization rate constant of amino acids, mol/(mol h)
Ks_aa<-6.4e-3 #Monod constant associated with amino acid utilization, mol/L 
km_h2<-16 #Maximum specific utilization rate constant of hydrogen, mol/(mol h)
Ks_h2<-5.84e-6 #Monod constant associated with hydrogen utilization, mol/L
k_br<-0.095 #Kinetic rate constant of bromoform utilization, 1/h
p1<-7.2551e+04 #Sigmoid function parameter associated with bromoform inhibition factor (Ibr) 
p2<--1.0837e-04 #Sigmoid function parameter associated with bromoform inhibition factor (Ibr)
p3<-0.3655 #Affine function parameter associated with flux allocation parameter towards acetate production (lambda1)
p4<-0.6371 #Affine function parameter associated with flux allocation parameter towards acetate production (lambda1)
p5<-0.3787 #Affine function parameter associated with flux allocation parameter towards propionate production (lambda2)   
p6<-0.1160 #Affine function parameter associated with flux allocation parameter towards propionate production (lambda2) 
Inputparams<-c(khyd_ndf=khyd_ndf,khyd_nsc=khyd_nsc,khyd_pro=khyd_pro,
               km_su=km_su,Ks_su=Ks_su,
               km_aa=km_aa,Ks_aa=Ks_aa,
               km_h2=km_h2,Ks_h2=Ks_h2,
               k_br=k_br,
               p1=p1,p2=p2,
               p3=p3,p4=p4,
               p5=p5,p6=p6)

#Time scale considered
nd<-4 #number of days simulated
ts<-1/60 #step time, hour
t<-seq(0,24*nd,by=ts)#time in hour

#Rumen volume in Rusitec condition
V_l<-0.74 #Volume in liquid phase of the rumen, L (Communication Belanche et al., 2017)
V_g<-0.06 #Volume in gas phase of the rumen, L (Communication Belanche et al., 2017)

# ------ Solving the model for the 3 dietary scenarios ------- 

# ---- Control (fraction of Asparagopsis taxiformis in the feed = 0%) ----

XDietChagasControl<-Cinit_Control #Initial conditions of the system
auxsDietChagasControl<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrControl) 
XmDietChagasControl<-data.frame(ode(y=XDietChagasControl,times=t,func=rumencAT,
                                    parms=auxsDietChagasControl,method="lsode"))

#Biochemical components produced from rumen fermentation (output variables of the model, quantity)
Z_ndfDietChagasControl<-XmDietChagasControl[,2] #Neutral detergent fiber (NDF), g
Z_nscDietChagasControl<-XmDietChagasControl[,3] #Non-structural carbohydrates (NSC), g
Z_proDietChagasControl<-XmDietChagasControl[,4] #Proteins, g 
S_suDietChagasControl<-XmDietChagasControl[,5] #Sugars, mol
S_aaDietChagasControl<-XmDietChagasControl[,6] #Average amino acids, mol

S_acDietChagasControl<-XmDietChagasControl[,7] #Total acetate, mol
S_buDietChagasControl<-XmDietChagasControl[,8] #Total butyrate, mol
S_prDietChagasControl<-XmDietChagasControl[,9] #Total propionate, mol
S_INDietChagasControl<-XmDietChagasControl[,10] #Inorganic nitrogen, mol
S_ICDietChagasControl<-XmDietChagasControl[,11] #Inorganic carbon, mol
S_h2DietChagasControl<-XmDietChagasControl[,12] #Hydrogen in liquid phase, mol
S_ch4DietChagasControl<-XmDietChagasControl[,13] #Methane in the liquid phase, mol
X_suDietChagasControl<-XmDietChagasControl[,14] #Sugars-utilizing microbes, mol
X_aaDietChagasControl<-XmDietChagasControl[,15] #Amino acids-utilizing microbes, mol
X_h2DietChagasControl<-XmDietChagasControl[,16] #Hydrogen-utilizing microbes (methanogens), mol

Sg_co2DietChagasControl<-XmDietChagasControl[,17] #Hydrogen in the gas phase, mol
Sg_h2DietChagasControl<-XmDietChagasControl[,18] #Carbon dioxide in the gas phase, mol
Sg_ch4DietChagasControl<-XmDietChagasControl[,19] #Methane in the gas phase, mol

S_brDietChagasControl<-XmDietChagasControl[,20] #Bromoform, g

#Conversion of the output variables in concentration (quantity/volume)
#Division by rumen volume
z_ndfDietChagasControl<-Z_ndfDietChagasControl/V_l #Neutral detergent fiber (NDF) concentration, g/L
z_nscDietChagasControl<-Z_nscDietChagasControl/V_l #Non-structural carbohydrates (NSC) concentration, g/L
z_proDietChagasControl<-Z_proDietChagasControl/V_l #Proteins concentration, g/L 
s_suDietChagasControl<-S_suDietChagasControl/V_l #Sugars concentration, mol/L
s_aaDietChagasControl<-S_aaDietChagasControl/V_l #Average amino acids concentration, mol/L
s_acDietChagasControl<-S_acDietChagasControl/V_l #Total acetate concentration, mol/L
s_buDietChagasControl<-S_buDietChagasControl/V_l #Total butyrate concentration, mol/L
s_prDietChagasControl<-S_prDietChagasControl/V_l #Total propionate concentration, mol/L
s_INDietChagasControl<-S_INDietChagasControl/V_l #Inorganic nitrogen concentration, mol/L
s_ICDietChagasControl<-S_ICDietChagasControl/V_l #Inorganic carbon concentration, mol/L
s_h2DietChagasControl<-S_h2DietChagasControl/V_l #Hydrogen concentration in liquid phase, mol/L
s_ch4DietChagasControl<-S_ch4DietChagasControl/V_l #Methane concentration in the liquid phase, mol/L
x_suDietChagasControl<-X_suDietChagasControl/V_l #Concentration of sugars-utilizing microbes, mol/L
x_aaDietChagasControl<-X_aaDietChagasControl/V_l #Concentration of amino acids-utilizing microbes, mol/L
x_h2DietChagasControl<-X_h2DietChagasControl/V_l #Concentration of hydrogen-utilizing microbes (methanogens), mol/L
sg_co2DietChagasControl<-Sg_co2DietChagasControl/V_g #Concentration of hydrogen in the gas phase, mol/L
sg_h2DietChagasControl<-Sg_h2DietChagasControl/V_g #Concentration of carbon dioxide in the gas phase, mol/L
sg_ch4DietChagasControl<-Sg_ch4DietChagasControl/V_g #Concentration of methane in the gas phase mol/L
s_brDietChagasControl<-S_brDietChagasControl/V_l #Bromoform concentration, g/L

XmconcentrationDietChagasControl<-cbind.data.frame(time=t,
                                                   z_ndfDietChagasControl,z_nscDietChagasControl,z_proDietChagasControl,
                                                   s_suDietChagasControl,s_aaDietChagasControl,
                                                   s_acDietChagasControl,s_buDietChagasControl,s_prDietChagasControl,
                                                   s_INDietChagasControl,s_ICDietChagasControl,s_h2DietChagasControl,s_ch4DietChagasControl,
                                                   x_suDietChagasControl,x_aaDietChagasControl,x_h2DietChagasControl,
                                                   sg_co2DietChagasControl,sg_h2DietChagasControl,sg_ch4DietChagasControl,
                                                   s_brDietChagasControl)

#Plots of biochemical component concentration produced from rumen fermentation
Plotz_ndfDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=z_ndfDietChagasControl))+
  labs(title="",x="Time (h)",y="NDF (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_nscDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=z_nscDietChagasControl))+
  labs(title="",x="Time (h)",y="NSC (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_proDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=z_proDietChagasControl))+
  labs(title="",x="Time (h)",y="Proteins (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_suDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_suDietChagasControl))+
  labs(title="",x="Time (h)",y="Sugars (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_aaDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_aaDietChagasControl))+
  labs(title="",x="Time (h)",y="Aminoacids (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotz_ndfDietChagasControl,Plotz_nscDietChagasControl,Plotz_proDietChagasControl,
                  Plots_suDietChagasControl,Plots_aaDietChagasControl,
                  nrow=2,ncol=3)

Plots_acDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_acDietChagasControl))+
  labs(title="",x="Time (h)",y="Acetate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_buDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_buDietChagasControl))+
  labs(title="",x="Time (h)",y="Butyrate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_prDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_prDietChagasControl))+
  labs(title="",x="Time (h)",y="Propionate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_INDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_INDietChagasControl))+
  labs(title="",x="Time (h)",y="NH3 (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ICDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_ICDietChagasControl))+
  labs(title="",x="Time (h)",y="CO2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_h2DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_h2DietChagasControl))+
  labs(title="",x="Time (h)",y="H2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ch4DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_ch4DietChagasControl))+
  labs(title="",x="Time (h)",y="CH4 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_suDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=x_suDietChagasControl))+
  labs(title="",x="Time (h)",y="Sugars utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_aaDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=x_aaDietChagasControl))+
  labs(title="",x="Time (h)",y="Amino acids utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_h2DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=x_h2DietChagasControl))+
  labs(title="",x="Time (h)",y="H2 utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_brDietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=s_brDietChagasControl))+
  labs(title="",x="Time (h)",y="Bromoform (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plots_acDietChagasControl,Plots_buDietChagasControl,Plots_prDietChagasControl,
                  Plots_INDietChagasControl,Plots_ICDietChagasControl,Plots_h2DietChagasControl,Plots_ch4DietChagasControl,
                  Plotx_suDietChagasControl,Plotx_aaDietChagasControl,Plotx_h2DietChagasControl,Plots_brDietChagasControl,
                  nrow=3,ncol=4)

Plotsg_co2DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=sg_co2DietChagasControl))+
  labs(title="",x="Time (h)",y="CO2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_h2DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=sg_h2DietChagasControl))+
  labs(title="",x="Time (h)",y="H2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_ch4DietChagasControl<-ggplot(XmconcentrationDietChagasControl,aes(x=time,y=sg_ch4DietChagasControl))+
  labs(title="",x="Time (h)",y="CH4 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotsg_co2DietChagasControl,Plotsg_h2DietChagasControl,Plotsg_ch4DietChagasControl,
                  nrow=1,ncol=3)

# ---- Low treatment (fraction of Asparagopsis taxiformis in the feed = 0.25%) ----

XDietChagasLow<-Cinit_Control #Initial conditions of the system
auxsDietChagasLow<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrLow)
XmDietChagasLow<-data.frame(ode(y=XDietChagasLow,times=t,func=rumencAT,
                                parms=auxsDietChagasLow,method="lsode"))

#Checking that simulated parameters p1 to p8 provide reliable values for I_br, I_H2, lambda_1 and lambda_2
Plotpgas_h2lambda_1DietChagasLow<-ggplot(XmDietChagasLow,aes(x=pgas_h2,y=lambda_1))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[1]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

Plotpgas_h2lambda_2DietChagasLow<-ggplot(XmDietChagasLow,aes(x=pgas_h2,y=lambda_2))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[2]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

Plotpgas_h2I_H2DietChagasLow<-ggplot(XmDietChagasLow,aes(x=pgas_h2,y=I_H2))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=bquote(~I[H[2]]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

Plots_brI_brDietChagasLow<-ggplot(XmDietChagasLow,aes(x=S_br,y=I_br.p1))+
  labs(title="",x=bquote(~S[br]~"(g)"),y=bquote(~I[br]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

ggpubr::ggarrange(Plotpgas_h2lambda_1DietChagasLow,Plotpgas_h2lambda_2DietChagasLow,
                  Plotpgas_h2I_H2DietChagasLow,Plots_brI_brDietChagasLow,
                  nrow=2,ncol=2)

#Biochemical components produced from rumen fermentation (output variables of the model, quantity)
Z_ndfDietChagasLow<-XmDietChagasLow[,2] #Neutral detergent fiber (NDF), g
Z_nscDietChagasLow<-XmDietChagasLow[,3] #Non-structural carbohydrates (NSC), g
Z_proDietChagasLow<-XmDietChagasLow[,4] #Proteins, g 
S_suDietChagasLow<-XmDietChagasLow[,5] #Sugars, mol
S_aaDietChagasLow<-XmDietChagasLow[,6] #Average amino acids, mol

S_acDietChagasLow<-XmDietChagasLow[,7] #Total acetate, mol
S_buDietChagasLow<-XmDietChagasLow[,8] #Total butyrate, mol
S_prDietChagasLow<-XmDietChagasLow[,9] #Total propionate, mol
S_INDietChagasLow<-XmDietChagasLow[,10] #Inorganic nitrogen, mol
S_ICDietChagasLow<-XmDietChagasLow[,11] #Inorganic carbon, mol
S_h2DietChagasLow<-XmDietChagasLow[,12] #Hydrogen in liquid phase, mol
S_ch4DietChagasLow<-XmDietChagasLow[,13] #Methane in the liquid phase, mol
X_suDietChagasLow<-XmDietChagasLow[,14] #Sugars-utilizing microbes, mol
X_aaDietChagasLow<-XmDietChagasLow[,15] #Amino acids-utilizing microbes, mol
X_h2DietChagasLow<-XmDietChagasLow[,16] #Hydrogen-utilizing microbes (methanogens), mol

Sg_co2DietChagasLow<-XmDietChagasLow[,17] #Hydrogen in the gas phase, mol
Sg_h2DietChagasLow<-XmDietChagasLow[,18] #Carbon dioxide in the gas phase, mol
Sg_ch4DietChagasLow<-XmDietChagasLow[,19] #Methane in the gas phase, mol

S_brDietChagasLow<-XmDietChagasLow[,20] #Bromoform, g

#Conversion of the output variables in concentration (quantity/volume)
#Division by rumen volume
z_ndfDietChagasLow<-Z_ndfDietChagasLow/V_l #Neutral detergent fiber (NDF) concentration, g/L
z_nscDietChagasLow<-Z_nscDietChagasLow/V_l #Non-structural carbohydrates (NSC) concentration, g/L
z_proDietChagasLow<-Z_proDietChagasLow/V_l #Proteins concentration, g/L 
s_suDietChagasLow<-S_suDietChagasLow/V_l #Sugars concentration, mol/L
s_aaDietChagasLow<-S_aaDietChagasLow/V_l #Average amino acids concentration, mol/L
s_acDietChagasLow<-S_acDietChagasLow/V_l #Total acetate concentration, mol/L
s_buDietChagasLow<-S_buDietChagasLow/V_l #Total butyrate concentration, mol/L
s_prDietChagasLow<-S_prDietChagasLow/V_l #Total propionate concentration, mol/L
s_INDietChagasLow<-S_INDietChagasLow/V_l #Inorganic nitrogen concentration, mol/L
s_ICDietChagasLow<-S_ICDietChagasLow/V_l #Inorganic carbon concentration, mol/L
s_h2DietChagasLow<-S_h2DietChagasLow/V_l #Hydrogen concentration in liquid phase, mol/L
s_ch4DietChagasLow<-S_ch4DietChagasLow/V_l #Methane concentration in the liquid phase, mol/L
x_suDietChagasLow<-X_suDietChagasLow/V_l #Concentration of sugars-utilizing microbes, mol/L
x_aaDietChagasLow<-X_aaDietChagasLow/V_l #Concentration of amino acids-utilizing microbes, mol/L
x_h2DietChagasLow<-X_h2DietChagasLow/V_l #Concentration of hydrogen-utilizing microbes (methanogens), mol/L
sg_co2DietChagasLow<-Sg_co2DietChagasLow/V_g #Concentration of hydrogen in the gas phase, mol/L
sg_h2DietChagasLow<-Sg_h2DietChagasLow/V_g #Concentration of carbon dioxide in the gas phase, mol/L
sg_ch4DietChagasLow<-Sg_ch4DietChagasLow/V_g #Concentration of methane in the gas phase mol/L
s_brDietChagasLow<-S_brDietChagasLow/V_l #Bromoform concentration, g/L

XmconcentrationDietChagasLow<-cbind.data.frame(time=t,
                                                 z_ndfDietChagasLow,z_nscDietChagasLow,z_proDietChagasLow,
                                                 s_suDietChagasLow,s_aaDietChagasLow,
                                                 s_acDietChagasLow,s_buDietChagasLow,s_prDietChagasLow,
                                                 s_INDietChagasLow,s_ICDietChagasLow,s_h2DietChagasLow,s_ch4DietChagasLow,
                                                 x_suDietChagasLow,x_aaDietChagasLow,x_h2DietChagasLow,
                                                 sg_co2DietChagasLow,sg_h2DietChagasLow,sg_ch4DietChagasLow,
                                                 s_brDietChagasLow)

#Plots of biochemical component concentration produced from rumen fermentation
Plotz_ndfDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=z_ndfDietChagasLow))+
  labs(title="",x="Time (h)",y="NDF (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_nscDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=z_nscDietChagasLow))+
  labs(title="",x="Time (h)",y="NSC (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_proDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=z_proDietChagasLow))+
  labs(title="",x="Time (h)",y="Proteins (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_suDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_suDietChagasLow))+
  labs(title="",x="Time (h)",y="Sugars (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_aaDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_aaDietChagasLow))+
  labs(title="",x="Time (h)",y="Aminoacids (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotz_ndfDietChagasLow,Plotz_nscDietChagasLow,Plotz_proDietChagasLow,
                  Plots_suDietChagasLow,Plots_aaDietChagasLow,
                  nrow=2,ncol=3)

Plots_acDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_acDietChagasLow))+
  labs(title="",x="Time (h)",y="Acetate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_buDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_buDietChagasLow))+
  labs(title="",x="Time (h)",y="Butyrate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_prDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_prDietChagasLow))+
  labs(title="",x="Time (h)",y="Propionate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_INDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_INDietChagasLow))+
  labs(title="",x="Time (h)",y="NH3 (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ICDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_ICDietChagasLow))+
  labs(title="",x="Time (h)",y="CO2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_h2DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_h2DietChagasLow))+
  labs(title="",x="Time (h)",y="H2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ch4DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_ch4DietChagasLow))+
  labs(title="",x="Time (h)",y="CH4 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_suDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=x_suDietChagasLow))+
  labs(title="",x="Time (h)",y="Sugars utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_aaDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=x_aaDietChagasLow))+
  labs(title="",x="Time (h)",y="Amino acids utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_h2DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=x_h2DietChagasLow))+
  labs(title="",x="Time (h)",y="H2 utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_brDietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=s_brDietChagasLow))+
  labs(title="",x="Time (h)",y="Bromoform (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plots_acDietChagasLow,Plots_buDietChagasLow,Plots_prDietChagasLow,
                  Plots_INDietChagasLow,Plots_ICDietChagasLow,Plots_h2DietChagasLow,Plots_ch4DietChagasLow,
                  Plotx_suDietChagasLow,Plotx_aaDietChagasLow,Plotx_h2DietChagasLow,Plots_brDietChagasLow,
                  nrow=3,ncol=4)

Plotsg_co2DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=sg_co2DietChagasLow))+
  labs(title="",x="Time (h)",y="CO2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_h2DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=sg_h2DietChagasLow))+
  labs(title="",x="Time (h)",y="H2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_ch4DietChagasLow<-ggplot(XmconcentrationDietChagasLow,aes(x=time,y=sg_ch4DietChagasLow))+
  labs(title="",x="Time (h)",y="CH4 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotsg_co2DietChagasLow,Plotsg_h2DietChagasLow,Plotsg_ch4DietChagasLow,
                  nrow=1,ncol=3)

# ---- High treatment (fraction of Asparagopsis taxiformis in the feed = 0.50%) ----

XDietChagasHigh<-Cinit_Control #Initial conditions of the system
auxsDietChagasHigh<-c(wNDF=wNDFDietChagas,wNSC=wNSCDietChagas,wPRO=wPRODietChagas,wBr=wBrHigh)
XmDietChagasHigh<-data.frame(ode(y=XDietChagasHigh,times=t,func=rumencAT,
                                 parms=auxsDietChagasHigh,method="lsode"))

#Checking that simulated parameters p1 to p6 provide reliable values for I_br, I_H2, lambda_1 and lambda_2
Plotpgas_h2lambda_1DietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=pgas_h2,y=lambda_1))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[1]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

Plotpgas_h2lambda_2DietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=pgas_h2,y=lambda_2))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[2]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

Plotpgas_h2I_H2DietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=pgas_h2,y=I_H2))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=bquote(~I[H[2]]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

Plots_brI_brDietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=S_br,y=I_br.p1))+
  labs(title="",x=bquote(~S[br]~"(g)"),y=bquote(~I[br]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=30,colour="black"),
        axis.text.y=element_text(size=30,colour="black"))

ggpubr::ggarrange(Plotpgas_h2lambda_1DietChagasHigh,Plotpgas_h2lambda_2DietChagasHigh,
                  Plotpgas_h2I_H2DietChagasHigh,Plots_brI_brDietChagasHigh,
                  nrow=2,ncol=2)

#Biochemical components produced from rumen fermentation (output variables of the model, quantity)
Z_ndfDietChagasHigh<-XmDietChagasHigh[,2] #Neutral detergent fiber (NDF), g
Z_nscDietChagasHigh<-XmDietChagasHigh[,3] #Non-structural carbohydrates (NSC), g
Z_proDietChagasHigh<-XmDietChagasHigh[,4] #Proteins, g 
S_suDietChagasHigh<-XmDietChagasHigh[,5] #Sugars, mol
S_aaDietChagasHigh<-XmDietChagasHigh[,6] #Average amino acids, mol

S_acDietChagasHigh<-XmDietChagasHigh[,7] #Total acetate, mol
S_buDietChagasHigh<-XmDietChagasHigh[,8] #Total butyrate, mol
S_prDietChagasHigh<-XmDietChagasHigh[,9] #Total propionate, mol
S_INDietChagasHigh<-XmDietChagasHigh[,10] #Inorganic nitrogen, mol
S_ICDietChagasHigh<-XmDietChagasHigh[,11] #Inorganic carbon, mol
S_h2DietChagasHigh<-XmDietChagasHigh[,12] #Hydrogen in liquid phase, mol
S_ch4DietChagasHigh<-XmDietChagasHigh[,13] #Methane in the liquid phase, mol
X_suDietChagasHigh<-XmDietChagasHigh[,14] #Sugars-utilizing microbes, mol
X_aaDietChagasHigh<-XmDietChagasHigh[,15] #Amino acids-utilizing microbes, mol
X_h2DietChagasHigh<-XmDietChagasHigh[,16] #Hydrogen-utilizing microbes (methanogens), mol

Sg_co2DietChagasHigh<-XmDietChagasHigh[,17] #Hydrogen in the gas phase, mol
Sg_h2DietChagasHigh<-XmDietChagasHigh[,18] #Carbon dioxide in the gas phase, mol
Sg_ch4DietChagasHigh<-XmDietChagasHigh[,19] #Methane in the gas phase, mol

S_brDietChagasHigh<-XmDietChagasHigh[,20] #Bromoform, g

#Conversion of the output variables in concentration (quantity/volume)
#Division by rumen volume
z_ndfDietChagasHigh<-Z_ndfDietChagasHigh/V_l #Neutral detergent fiber (NDF) concentration, g/L
z_nscDietChagasHigh<-Z_nscDietChagasHigh/V_l #Non-structural carbohydrates (NSC) concentration, g/L
z_proDietChagasHigh<-Z_proDietChagasHigh/V_l #Proteins concentration, g/L 
s_suDietChagasHigh<-S_suDietChagasHigh/V_l #Sugars concentration, mol/L
s_aaDietChagasHigh<-S_aaDietChagasHigh/V_l #Average amino acids concentration, mol/L
s_acDietChagasHigh<-S_acDietChagasHigh/V_l #Total acetate concentration, mol/L
s_buDietChagasHigh<-S_buDietChagasHigh/V_l #Total butyrate concentration, mol/L
s_prDietChagasHigh<-S_prDietChagasHigh/V_l #Total propionate concentration, mol/L
s_INDietChagasHigh<-S_INDietChagasHigh/V_l #Inorganic nitrogen concentration, mol/L
s_ICDietChagasHigh<-S_ICDietChagasHigh/V_l #Inorganic carbon concentration, mol/L
s_h2DietChagasHigh<-S_h2DietChagasHigh/V_l #Hydrogen concentration in liquid phase, mol/L
s_ch4DietChagasHigh<-S_ch4DietChagasHigh/V_l #Methane concentration in the liquid phase, mol/L
x_suDietChagasHigh<-X_suDietChagasHigh/V_l #Concentration of sugars-utilizing microbes, mol/L
x_aaDietChagasHigh<-X_aaDietChagasHigh/V_l #Concentration of amino acids-utilizing microbes, mol/L
x_h2DietChagasHigh<-X_h2DietChagasHigh/V_l #Concentration of hydrogen-utilizing microbes (methanogens), mol/L
sg_co2DietChagasHigh<-Sg_co2DietChagasHigh/V_g #Concentration of hydrogen in the gas phase, mol/L
sg_h2DietChagasHigh<-Sg_h2DietChagasHigh/V_g #Concentration of carbon dioxide in the gas phase, mol/L
sg_ch4DietChagasHigh<-Sg_ch4DietChagasHigh/V_g #Concentration of methane in the gas phase mol/L
s_brDietChagasHigh<-S_brDietChagasHigh/V_l #Bromoform concentration, g/L

XmconcentrationDietChagasHigh<-cbind.data.frame(time=t,
                                                z_ndfDietChagasHigh,z_nscDietChagasHigh,z_proDietChagasHigh,
                                                s_suDietChagasHigh,s_aaDietChagasHigh,
                                                s_acDietChagasHigh,s_buDietChagasHigh,s_prDietChagasHigh,
                                                s_INDietChagasHigh,s_ICDietChagasHigh,s_h2DietChagasHigh,s_ch4DietChagasHigh,
                                                x_suDietChagasHigh,x_aaDietChagasHigh,x_h2DietChagasHigh,
                                                sg_co2DietChagasHigh,sg_h2DietChagasHigh,sg_ch4DietChagasHigh,
                                                s_brDietChagasHigh)

#Checking that simulated parameters p1 to p6 provide reliable values for I_br, I_H2, lambda_1 and lambda_2
Plotpgas_h2lambda_1DietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=pgas_h2,y=lambda_1))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[1]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

Plotpgas_h2lambda_2DietChagasHigh<-ggplot(XmDietChagasHigh,aes(x=pgas_h2,y=lambda_2))+
  labs(title="",x=bquote(~p[H[2]]~"(bar)"),y=expression(lambda[2]))+
  geom_line(size=1)+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

Plots_brI_brDietChagasHigh<-ggplot()+
  geom_line(aes(x=XmconcentrationDietChagasHigh$s_brDietChagasHigh*1000,y=XmDietChagasHigh$I_br.p1),size=1)+
  labs(title="",x=bquote(~s[br]~"(mg/L)"),y=bquote(~I[br]))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))

ggpubr::ggarrange(Plots_brI_brDietChagasHigh,Plotpgas_h2lambda_1DietChagasHigh,Plotpgas_h2lambda_2DietChagasHigh,
                  nrow=1,ncol=3)

#Plots of biochemical component concentration produced from rumen fermentation
Plotz_ndfDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=z_ndfDietChagasHigh))+
  labs(title="",x="Time (h)",y="NDF (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_nscDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=z_nscDietChagasHigh))+
  labs(title="",x="Time (h)",y="NSC (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotz_proDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=z_proDietChagasHigh))+
  labs(title="",x="Time (h)",y="Proteins (g/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_suDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_suDietChagasHigh))+
  labs(title="",x="Time (h)",y="Sugars (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_aaDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_aaDietChagasHigh))+
  labs(title="",x="Time (h)",y="Aminoacids (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotz_ndfDietChagasHigh,Plotz_nscDietChagasHigh,Plotz_proDietChagasHigh,
                  Plots_suDietChagasHigh,Plots_aaDietChagasHigh,
                  nrow=2,ncol=3)

Plots_acDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_acDietChagasHigh))+
  labs(title="",x="Time (h)",y="Acetate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_buDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_buDietChagasHigh))+
  labs(title="",x="Time (h)",y="Butyrate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_prDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_prDietChagasHigh))+
  labs(title="",x="Time (h)",y="Propionate (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_INDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_INDietChagasHigh))+
  labs(title="",x="Time (h)",y="NH3 (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ICDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_ICDietChagasHigh))+
  labs(title="",x="Time (h)",y="CO2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_h2DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_h2DietChagasHigh))+
  labs(title="",x="Time (h)",y="H2 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_ch4DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_ch4DietChagasHigh))+
  labs(title="",x="Time (h)",y="CH4 liquid (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_suDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=x_suDietChagasHigh))+
  labs(title="",x="Time (h)",y="Sugars utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_aaDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=x_aaDietChagasHigh))+
  labs(title="",x="Time (h)",y="Amino acids utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotx_h2DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=x_h2DietChagasHigh))+
  labs(title="",x="Time (h)",y="H2 utilisers (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plots_brDietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=s_brDietChagasHigh))+
  labs(title="",x="Time (h)",y="Bromoform (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plots_acDietChagasHigh,Plots_buDietChagasHigh,Plots_prDietChagasHigh,
                  Plots_INDietChagasHigh,Plots_ICDietChagasHigh,Plots_h2DietChagasHigh,Plots_ch4DietChagasHigh,
                  Plotx_suDietChagasHigh,Plotx_aaDietChagasHigh,Plotx_h2DietChagasHigh,Plots_brDietChagasHigh,
                  nrow=3,ncol=4)

Plotsg_co2DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=sg_co2DietChagasHigh))+
  labs(title="",x="Time (h)",y="CO2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_h2DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=sg_h2DietChagasHigh))+
  labs(title="",x="Time (h)",y="H2 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))

Plotsg_ch4DietChagasHigh<-ggplot(XmconcentrationDietChagasHigh,aes(x=time,y=sg_ch4DietChagasHigh))+
  labs(title="",x="Time (h)",y="CH4 gas (mol/L)")+
  geom_line(size=1)+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))


ggpubr::ggarrange(Plotsg_co2DietChagasHigh,Plotsg_h2DietChagasHigh,Plotsg_ch4DietChagasHigh,
                  nrow=1,ncol=3)

# ---- Comparison of the 3 dietary scenarios (Control, Low Asparagopsis taxiformis treatment, High Asparagopsis taxiformis treatment) ----

ColorsCompTreatment<-c("Control"="black","Low"="red","High"="green")
#Dynamic of methane output flow of gas phase (mol/h)
Plotq_ch4g_outDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+
  geom_line(aes(x=t,y=XmDietChagasControl$q_ch4g_out,color="Control"),size=2)+
  geom_line(aes(x=t,y=XmDietChagasLow$q_ch4g_out,color="Low"),size=2)+
  geom_line(aes(x=t,y=XmDietChagasHigh$q_ch4g_out,color="High"),size=2)+
  labs(title="",
       x="Time (h)",y=bquote(~q[CH[4]~",g,out"]~"(mol/h)"),
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=40),
        legend.text=element_text(size=35),
        axis.text.x=element_text(size=35,colour="black"),
        axis.text.y=element_text(size=35,colour="black"))
Plotq_ch4g_outDietChagasControlLowHighComp

#Computation of CH4 reduction percentage between the control and 2 Asparagopsis taxiformis treatments
#Conversion of CH4 from mol/h to g/h (1 mol = 16 g of CH4)
q_ch4g_outDietChagasControlgh<-XmDietChagasControl$q_ch4g_out*16
q_ch4g_outDietChagasLowgh<-XmDietChagasLow$q_ch4g_out*16
q_ch4g_outDietChagasHighgh<-XmDietChagasHigh$q_ch4g_out*16

Plotq_ch4g_outDietChagasControlLowHighCompgh<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=4))+
  geom_line(aes(x=t,y=q_ch4g_outDietChagasControlgh,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=q_ch4g_outDietChagasLowgh,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=q_ch4g_outDietChagasHighgh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y=bquote(~q[CH[4]~",g,out"]~"(g/h)"),
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=20,colour="black"),
        axis.text.y=element_text(size=20,colour="black"))
Plotq_ch4g_outDietChagasControlLowHighCompgh

#Computation of CH4 in g/day (computation of the area under the cure with the trapezoidal method)
q_ch4g_outDietChagasControlgd<-trapz(t,q_ch4g_outDietChagasControlgh) 
q_ch4g_outDietChagasLowgd<-trapz(t,q_ch4g_outDietChagasLowgh) 
q_ch4g_outDietChagasHighgd<-trapz(t,q_ch4g_outDietChagasHighgh) 

#Percentage of reduction between Control and Low/High Asparagopsis taxiformis treatment
PerRedDietChagasControlLow<-((q_ch4g_outDietChagasControlgd-q_ch4g_outDietChagasLowgd)/q_ch4g_outDietChagasControlgd)*100 #Low treatment
PerRedDietChagasControlHigh<-((q_ch4g_outDietChagasControlgd-q_ch4g_outDietChagasHighgd)/q_ch4g_outDietChagasControlgd)*100 #High treatment

#Computation of vfa fractions, %C2/%C3, %C3/%C2 and (%C2 + %C4)/%C3
#Control
x_acDietChagasControl<-s_acDietChagasControl/(s_acDietChagasControl+s_buDietChagasControl+s_prDietChagasControl)
x_buDietChagasControl<-s_buDietChagasControl/(s_acDietChagasControl+s_buDietChagasControl+s_prDietChagasControl)
x_prDietChagasControl<-s_prDietChagasControl/(s_acDietChagasControl+s_buDietChagasControl+s_prDietChagasControl)
x_ac_pr_ratioDietChagasControl<-x_acDietChagasControl/x_prDietChagasControl
x_pr_ac_ratioDietChagasControl<-x_prDietChagasControl/x_acDietChagasControl
x_ac_bu_pr_ratioDietChagasControl<-(x_acDietChagasControl+x_buDietChagasControl)/x_prDietChagasControl
#Low treatment
x_acDietChagasLow<-s_acDietChagasLow/(s_acDietChagasLow+s_buDietChagasLow+s_prDietChagasLow)
x_buDietChagasLow<-s_buDietChagasLow/(s_acDietChagasLow+s_buDietChagasLow+s_prDietChagasLow)
x_prDietChagasLow<-s_prDietChagasLow/(s_acDietChagasLow+s_buDietChagasLow+s_prDietChagasLow)
x_ac_pr_ratioDietChagasLow<-x_acDietChagasLow/x_prDietChagasLow
x_pr_ac_ratioDietChagasLow<-x_prDietChagasLow/x_acDietChagasLow
x_ac_bu_pr_ratioDietChagasLow<-(x_acDietChagasLow+x_buDietChagasLow)/x_prDietChagasLow
#High Treatment
x_acDietChagasHigh<-s_acDietChagasHigh/(s_acDietChagasHigh+s_buDietChagasHigh+s_prDietChagasHigh)
x_buDietChagasHigh<-s_buDietChagasHigh/(s_acDietChagasHigh+s_buDietChagasHigh+s_prDietChagasHigh)
x_prDietChagasHigh<-s_prDietChagasHigh/(s_acDietChagasHigh+s_buDietChagasHigh+s_prDietChagasHigh)
x_ac_pr_ratioDietChagasHigh<-x_acDietChagasHigh/x_prDietChagasHigh
x_pr_ac_ratioDietChagasHigh<-x_prDietChagasHigh/x_acDietChagasHigh
x_ac_bu_pr_ratioDietChagasHigh<-(x_acDietChagasHigh+x_buDietChagasHigh)/x_prDietChagasHigh

#Dynamic of acetate fraction (%) 
Plotx_acDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=8))+ylim(0.35,0.60)+
  geom_line(aes(x=t,y=x_acDietChagasControl,color="Control"),size=2)+
  geom_line(aes(x=t,y=x_acDietChagasLow,color="Low"),size=2)+
  geom_line(aes(x=t,y=x_acDietChagasHigh,color="High"),size=2)+
  labs(title="",
       x="Time (h)",y="Acetate fraction (%)",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=32),
        axis.title.y=element_text(size=32),
        plot.title=element_text(size=20),
        legend.title=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))
Plotx_acDietChagasControlLowHighComp

#Percentage of reduction between Control and Low/High Asparagopsis taxiformis treatment
AcetatePerRedDietChagasControlLow<-((x_acDietChagasControl[length(x_acDietChagasControl)]-x_acDietChagasLow[length(x_acDietChagasLow)])/x_acDietChagasControl[length(x_acDietChagasControl)])*100 #Low treatment
AcetatePerRedDietChagasControlHigh<-((x_acDietChagasControl[length(x_acDietChagasControl)]-x_acDietChagasHigh[length(x_acDietChagasHigh)])/x_acDietChagasControl[length(x_acDietChagasControl)])*100 #High treatment

#Dynamic of butyrate fraction (%) 
Plotx_buDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=8))+ylim(0.15,0.30)+
  geom_line(aes(x=t,y=x_buDietChagasControl,color="Control"),size=2)+
  geom_line(aes(x=t,y=x_buDietChagasLow,color="Low"),size=2)+
  geom_line(aes(x=t,y=x_buDietChagasHigh,color="High"),size=2)+
  labs(title="",
       x="Time (h)",y="Butyrate fraction (%)",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=32),
        axis.title.y=element_text(size=32),
        plot.title=element_text(size=20),
        legend.title=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))
Plotx_buDietChagasControlLowHighComp

#Percentage of increase between Control and Low/High Asparagopsis taxiformis treatment
ButyratePerIncDietChagasControlLow<-((x_buDietChagasLow[length(x_buDietChagasLow)]-x_buDietChagasControl[length(x_buDietChagasControl)])/x_buDietChagasControl[length(x_buDietChagasControl)])*100 #Low treatment
ButyratePerIncDietChagasControlHigh<-((x_buDietChagasHigh[length(x_buDietChagasHigh)]-x_buDietChagasControl[length(x_buDietChagasControl)])/x_buDietChagasControl[length(x_buDietChagasControl)])*100 #High treatment

#Dynamic of propionate fraction (%)
Plotx_prDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=8))+ylim(0.24,0.32)+
  geom_line(aes(x=t,y=x_prDietChagasControl,color="Control"),size=2)+
  geom_line(aes(x=t,y=x_prDietChagasLow,color="Low"),size=2)+
  geom_line(aes(x=t,y=x_prDietChagasHigh,color="High"),size=2)+
  labs(title="",
       x="Time (h)",y="Propionate fraction (%)",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=32),
        axis.title.y=element_text(size=32),
        plot.title=element_text(size=20),
        legend.title=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))
Plotx_prDietChagasControlLowHighComp

#Percentage of increase between Control and Low/High Asparagopsis taxiformis treatment
PropionatePerIncDietChagasControlLow<-((x_prDietChagasLow[length(x_prDietChagasLow)]-x_prDietChagasControl[length(x_prDietChagasControl)])/x_prDietChagasControl[length(x_prDietChagasControl)])*100 #Low treatment
propionatePerIncDietChagasControlHigh<-((x_prDietChagasHigh[length(x_prDietChagasHigh)]-x_prDietChagasControl[length(x_prDietChagasControl)])/x_prDietChagasControl[length(x_prDietChagasControl)])*100 #High treatment

#Dynamic of C2/C3 ratio
Plotx_ac_pr_ratioDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+
  geom_line(aes(x=t,y=x_ac_pr_ratioDietChagasControl,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=x_ac_pr_ratioDietChagasLow,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=x_ac_pr_ratioDietChagasHigh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y="Acetate/Propionate ratio",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))
Plotx_ac_pr_ratioDietChagasControlLowHighComp

#Dynamic of C3/C2 ratio
Plotx_pr_ac_ratioDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=8))+
  geom_line(aes(x=t,y=x_pr_ac_ratioDietChagasControl,color="Control"),size=2)+
  geom_line(aes(x=t,y=x_pr_ac_ratioDietChagasLow,color="Low"),size=2)+
  geom_line(aes(x=t,y=x_pr_ac_ratioDietChagasHigh,color="High"),size=2)+
  labs(title="",
       x="Time (h)",y="Propionate:Acetate ratio",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=32),
        axis.title.y=element_text(size=32),
        plot.title=element_text(size=20),
        legend.title=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.x=element_text(size=28,colour="black"),
        axis.text.y=element_text(size=28,colour="black"))
Plotx_pr_ac_ratioDietChagasControlLowHighComp

#Percentage of increase between Control and Low/High Asparagopsis taxiformis treatment
PropionateAcetateratioPerRedDietChagasControlLow<-((x_pr_ac_ratioDietChagasLow[length(x_pr_ac_ratioDietChagasLow)]-x_pr_ac_ratioDietChagasControl[length(x_pr_ac_ratioDietChagasControl)])/x_pr_ac_ratioDietChagasControl[length(x_pr_ac_ratioDietChagasControl)])*100 #Low treatment
PropionateAcetateratioPerRedDietChagasControlHigh<-((x_pr_ac_ratioDietChagasHigh[length(x_pr_ac_ratioDietChagasHigh)]-x_pr_ac_ratioDietChagasControl[length(x_pr_ac_ratioDietChagasControl)])/x_pr_ac_ratioDietChagasControl[length(x_pr_ac_ratioDietChagasControl)])*100 #High treatment

#Percentage of increase between Control and Low/High Asparagopsis taxiformis treatment
ButyratePerRedDietChagasControlLow<-((x_buDietChagasLow-x_buDietChagasControl)/x_buDietChagasControl)*100 #Low treatment
ButyratePerRedDietChagasControlHigh<-((x_buDietChagasHigh-x_buDietChagasControl)/x_buDietChagasControl)*100 #High treatment

#Dynamic of (C2+C4)/C3 ratio
Plotx_ac_bu_pr_ratioDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+
  geom_line(aes(x=t,y=x_ac_bu_pr_ratioDietChagasControl,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=x_ac_bu_pr_ratioDietChagasLow,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=x_ac_bu_pr_ratioDietChagasHigh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y="(Acetate+Butyrate)/Propionate ratio",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))
Plotx_ac_bu_pr_ratioDietChagasControlLowHighComp

#Dynamic of bromoform concentration (g/L)
Plots_brDietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+ylim(0,4e-04)+
  geom_line(aes(x=t,y=s_brDietChagasControl,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=s_brDietChagasLow,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=s_brDietChagasHigh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y="Bromoform (g/L)",
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))
Plots_brDietChagasControlLowHighComp

#Dynamic of hydrogen concentration in gas phase (mol/L)
Plotsg_h2DietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+
  geom_line(aes(x=t,y=sg_h2DietChagasControl,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=sg_h2DietChagasLow,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=sg_h2DietChagasHigh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y=bquote(~H[2]~"(mol/L)"),
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))
Plotsg_h2DietChagasControlLowHighComp

#Dynamic of methane concentration in gas phase (mol/L)
Plotsg_ch4DietChagasControlLowHighComp<-ggplot()+
  scale_x_continuous(breaks=seq(0,24*nd,by=6))+
  geom_line(aes(x=t,y=sg_ch4DietChagasControl,color="Control"),size=1.5)+
  geom_line(aes(x=t,y=sg_ch4DietChagasLow,color="Low"),size=1.5)+
  geom_line(aes(x=t,y=sg_ch4DietChagasHigh,color="High"),size=1.5)+
  labs(title="",
       x="Time (h)",y=bquote(~CH[4]~"(mol/L)"),
       color="")+
  scale_color_manual(values=ColorsCompTreatment,
                     labels=c("Control","High treatment","Low treatment"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x=element_text(size=15,colour="black"),
        axis.text.y=element_text(size=15,colour="black"))
Plotsg_ch4DietChagasControlLowHighComp

#Comparison of treatments for sensitivity analysis output variables
ggpubr::ggarrange(Plotx_acDietChagasControlLowHighComp,Plotx_buDietChagasControlLowHighComp,
                  Plotx_prDietChagasControlLowHighComp,Plotx_ac_pr_ratioDietChagasControlLowHighComp,
                  Plotx_ac_bu_pr_ratioDietChagasControlLowHighComp,Plots_brDietChagasControlLowHighComp,
                  Plotsg_h2DietChagasControlLowHighComp,Plotsg_ch4DietChagasControlLowHighComp,
                  nrow=4,ncol=2,common.legend=TRUE)

ggpubr::ggarrange(Plotx_acDietChagasControlLowHighComp,Plotx_buDietChagasControlLowHighComp,
                  Plotx_prDietChagasControlLowHighComp,Plotx_pr_ac_ratioDietChagasControlLowHighComp,
                  nrow=2,ncol=2,common.legend=TRUE)
