# -------------------------------------------------------------

# Function selecting the measured output variables of the model 

# -------------------------------------------------------------

rumencATout<-function(X,V_l,V_g) {

  #######################################################################
  # List of inputs to this function:
  # X:Matrix of the output variables simulated with the mechanistic model
  # V_l:Volume in the liquid phase of the rumen, L
  # V_g:Volume in the gas phase of the rumen, L
  #######################################################################
  
  time<-X[,"time"] #Time in hour
  
  #Biochemical components produced from the rumen fermentation (output variables of the model, quantity)
  #Polymer components
  Z_ndf<-X[,"Z_ndf"] #Neutral detergent fiber (NDF), g
  Z_nsc<-X[,"Z_nsc"] #Non-structural carbohydrates (NSC), g
  Z_pro<-X[,"Z_pro"] #Proteins, g 
  #Soluble components
  S_su<-X[,"S_su"] #Sugars, mol
  S_aa<-X[,"S_aa"] #Average amino acids, mol
  S_ac<-X[,"S_ac"] #Total acetate, mol
  S_bu<-X[,"S_bu"] #Total butyrate, mol
  S_pr<-X[,"S_pr"] #Total propionate, mol
  S_IN<-X[,"S_IN"] #Inorganic nitrogen, mol
  S_IC<-X[,"S_IC"] #Inorganic carbon, mol
  S_h2<-X[,"S_h2"] #Hydrogen in liquid phase, mol
  S_ch4<-X[,"S_ch4"] #Methane in the liquid phase, mol
  S_br<-X[,"S_br"] #Bromoform, g
  #Microbial functional groups
  X_su<-X[,"X_su"] #Sugars-utilizing microbes, mol
  X_aa<-X[,"X_aa"] #Amino acids-utilizing microbes, mol
  X_h2<-X[,"X_h2"] #Hydrogen-utilizing microbes (methanogens), mol
  #Gas phase
  Sg_co2<-X[,"Sg_co2"] #Hydrogen in the gas phase, mol
  Sg_h2<-X[,"Sg_h2"] #Carbon dioxide in the gas phase, mol
  Sg_ch4<-X[,"Sg_ch4"] #Methane in the gas phase, mol
  q_ch4g_out<-X[,"q_ch4g_out"] #Methane output flow of gas phase, mol/h
  
  
  #Conversion of the output variables in concentration (quantity/volume)
  #Polymer components
  z_ndf<-Z_ndf/V_l #Neutral detergent fiber (NDF) concentration, g/L
  z_nsc<-Z_nsc/V_l #Non-structural carbohydrates (NSC) concentration, g/L
  z_pro<-Z_pro/V_l #Proteins concentration, g/L 
  #Soluble components
  s_su<-S_su/V_l #Sugars concentration, mol/L
  s_aa<-S_aa/V_l #Average amino acids concentration, mol/L
  s_ac<-S_ac/V_l #Total acetate concentration, mol/L
  s_bu<-S_bu/V_l #Total butyrate concentration, mol/L
  s_pr<-S_pr/V_l #Total propionate concentration, mol/L
  s_IN<-S_IN/V_l #Inorganic nitrogen concentration, mol/L
  s_IC<-S_IC/V_l #Inorganic carbon concentration, mol/L
  s_h2<-S_h2/V_l #Hydrogen concentration in liquid phase, mol/L
  s_ch4<-S_ch4/V_l #Methane concentration in the liquid phase, mol/L
  s_br<-S_br/V_l #Bromoform concentration, g/L
  #Microbial functional groups
  x_su<-X_su/V_l #Concentration of sugars-utilizing microbes, mol/L
  x_aa<-X_aa/V_l #Concentration of amino acids-utilizing microbes, mol/L
  x_h2<-X_h2/V_l #Concentration of hydrogen-utilizing microbes (methanogens), mol/L
  #Gas phase
  sg_co2<-Sg_co2/V_g #Concentration of hydrogen in the gas phase, mol/L
  sg_h2<-Sg_h2/V_g #Concentration of carbon dioxide in the gas phase, mol/L
  sg_ch4<-Sg_ch4/V_g #Concentration of methane in the gas phase mol/L
  
  #Selection of the studied output variables (acetate, butyrate, propionate and methane)
  Ym<-cbind.data.frame(time,s_ac,s_bu,s_pr,sg_ch4,q_ch4g_out)
  
  return(Ym)
  
}
  
  