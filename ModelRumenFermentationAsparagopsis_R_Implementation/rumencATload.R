# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script loading the initial conditions of the rumen system and simulating the dynamic intake (Muñoz-Tamayo et al., 2019) and dietary scenarios (Chagas et al., 2019)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(ggplot2)

# ---- Initial conditions of the system (Set by running the model with the control during 50 days) ----

#Model ran with values from Rusitec condition (Belanche et al., 2017)
Cinit_Control<-read.table(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/InitialconditionsControl.txt",sep=""),header=T)
NamesCinit_Control<-colnames(Cinit_Control)
Cinit_Control<-as.numeric(Cinit_Control)
names(Cinit_Control)<-NamesCinit_Control

# ---- Simulation of the dynamic intake and setting of the diet composition scenarios ----

# -- Dynamic intake (Muñoz-Tamayo et al., 2019) --

#Intake scenario (Belanche et al., 2017 and expert knowledge) -
DMtotal<-11.25 #total dry matter ingested in one day, g/d
ni<-2 #number of intake in one day
DMtotalRep<-c(0.70,0.30) #Distribution of DMtotal between meals, %
k<-0.015 #the intake rate

#Plot of the dry matter intake
source(paste(paste(unlist(strsplit(getwd(),"/"))[-length(unlist(strsplit(getwd(),"/")))],collapse="/"),"/ModelRumenFermentationAsparagopsis_R_Implementation/DynamicIntake.R",sep="")) #function simulating the dynamic intake 
DMIdyn<-DynamicIntake(DMtotal,DMtotalRep,ni,k)
tDMI<-DMIdyn[,1] #time of measurements of DMI, h 
yDMI<-DMIdyn[,2] #measurements of DMI, g/h 

plotDMI<-ggplot()+
  geom_line(aes(x=tDMI,y=yDMI),size=1)+
  labs(title="",x="Time (h)",y="DMI (g/h)")+
  theme_classic(base_size=30)+
  theme(axis.title.x=element_text(size=40),
        axis.title.y=element_text(size=40),
        plot.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.x=element_text(size=40,colour="black"),
        axis.text.y=element_text(size=40,colour="black"))
print(plotDMI)

# -- Diet composition scenarios (Chagas et al., 2019) --

#Fraction of bromoform in the feed
x_br_AT<-6.84/1000 #Content of bromoform in Asparagopsis taxiformis, g/g
#Control
w_ATControl<-0 #fraction of Asparagopsis taxiformis in the feed, % 
wBrControl<-w_ATControl*x_br_AT #fraction of Bromoform in the feed, %
#Low Asparagopsis taxiformis treatment
w_ATLow<-0.25 #fraction of Asparagopsis taxiformis in the feed, % 
wBrLow<-w_ATLow*x_br_AT #fraction of Bromoform in the feed, % 
#High Asparagopsis taxiformis treatment
w_ATHigh<-0.5 #fraction of Asparagopsis taxiformis in the feed, % 
wBrHigh<-w_ATHigh*x_br_AT #fraction of Bromoform in the feed, % 

#Feed inputs composition
#Fractions of diet chemical composition
wNDFDietChagas<-38.7 #fraction of NDF in the feed, %
wPRODietChagas<-16 #fraction of Protein in the feed, %
wNSCDietChagas<-39.7 #fraction of NSC in the feed, %
wEEDietChagasControl<-100-(wNDFDietChagas+wPRODietChagas+wNSCDietChagas+w_ATControl) #fraction of esthers in the feed for the control, %
wEEDietChagasLow<-100-(wNDFDietChagas+wPRODietChagas+wNSCDietChagas+w_ATLow) #fraction of esthers in the feed for the low Asparagopsis taxiformis treatment, %
wEEDietChagasHigh<-100-(wNDFDietChagas+wPRODietChagas+wNSCDietChagas+w_ATHigh) #fraction of esthers in the feed for the high Asparagopsis taxiformis treatment, %