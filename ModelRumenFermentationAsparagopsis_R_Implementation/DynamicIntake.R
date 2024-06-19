# ---------------------------------------------------------------------------------------------------

# Implementation of the dynamic model of the dry matter intake (DMI, g/h) (Muñoz-Tamayo et al., 2019)

# ---------------------------------------------------------------------------------------------------

DynamicIntake<-function(DMtotal,DMtotalRep,ni,k) {

  #################################################################################
  # List of inputs to this function:
  # DMtotal:total dry matter ingested in one day, g
  # DMtotalRep:Distribution of DMtotal between meals, %
  # ni:number of intake in one day 
  # k:the intake rate
  #################################################################################
  
  #Defining the DM (g) dynamically
  time<-seq(0,24,by=1/60) #time in hour
  tm<-time%%(time[length(time)]/ni) #applies modulo to work with ni symetric responses
  
  #Computation of DM (g)
  DM<-c()
  for (i in 1:ni){
    
  DM_max<-DMtotalRep[i]*DMtotal #Maximum DMI for one meal (according to the distributions of the meals given by DMtotalRep)
  DM<-c(DM,DM_max*exp(-k*60*tm[((which(tm==0)[i]):(which(tm==0)[i+1]-1))])) #equation from Muñoz-Tamayo et al., 2019
  
  }
  DM[(time[length(time)]*60)+1]<-DMtotalRep[1]*DMtotal
  
  #Computation of DMI (g/h ; DMI = d/dt DM)
  DMI<-c(0) #initial condition of DMI
  for (j in 2:length(DM)){
    
  DMI[j]<-(DM[j-1]-DM[j])/(time[2]-time[1])
  DMI[j]<-max(0,DMI[j]) 
  
  }
  
  return(cbind(time,DMI))

}