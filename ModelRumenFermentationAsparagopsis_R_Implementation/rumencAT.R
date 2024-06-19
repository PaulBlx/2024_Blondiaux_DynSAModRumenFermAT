# ------------------------------------------------------------------------------------------------------------------------------------------------------

# Dynamic model of rumen fermentation under in vitro continuous condition accounting for the effect of Asparagopsis taxiformis (Muñoz-Tamayo et al 2021)

# ------------------------------------------------------------------------------------------------------------------------------------------------------

rumencAT<- function(t,X,auxs){
  
  #################################################################################
  # List of inputs to this function:
  # t: vector of time
  # X: Initial values of the State variables of the model
  # auxs:parameters used in the model 
  #################################################################################
  
  with(as.list(c(X,auxs)),{ 
    
      # ---- State variables of the model ----
    
      #Defining the intake dynamically with the equation of Muñoz-Tamayo et al., 2019
      tm<-t%%24
      DMI<-interp1(tDMI,yDMI,tm)
      
      #Defining the rumen volume 
      V_l<-0.74 #Volume in the liquid phase of the rumen, L (Communication Belanche et al., 2017)
      V_g<-0.06 #Volume in the gas phase of the rumen, L (Communication Belanche et al., 2017)
      
      z_ndf<-Z_ndf/V_l #Neutral detergent fiber (NDF) concentration, g/L
      z_nsc<-Z_nsc/V_l #Non-structural carbohydrates (NSC) concentration, g/L
      z_pro<-Z_pro/V_l #Proteins concentration, g/L 
      s_su<-S_su/V_l #Sugars concentration, mol/L
      s_aa<-S_aa/V_l #Average amino acids concentration, mol/L
      s_ac<-S_ac/V_l #Total acetate concentration, mol/L
      s_bu<-S_bu/V_l #Total butyrate concentration, mol/L
      s_pr<-S_pr/V_l #Total propionate concentration, mol/L
      s_IN<-S_IN/V_l #Inorganic nitrogen concentration, mol/L
      s_IC<-S_IC/V_l #Inorganic carbon concentration, mol/L
      s_h2<-S_h2/V_l #Hydrogen concentration in liquid phase, mol/L
      s_ch4<-S_ch4/V_l #Methane concentration in the liquid phase, mol/L
      x_su<-X_su/V_l #Concentration of sugars-utilizing microbes, mol/L
      x_aa<-X_aa/V_l #Concentration of amino acids-utilizing microbes, mol/L
      x_h2<-X_h2/V_l #Concentration of hydrogen-utilizing microbes (methanogens), mol/L
      sg_co2<-Sg_co2/V_g #Concentration of hydrogen in the gas phase, mol/L
      sg_h2<-Sg_h2/V_g #Concentration of carbon dioxide in the gas phase, mol/L
      sg_ch4<-Sg_ch4/V_g #Concentration of methane in the gas phase, mol/L
      s_br<-S_br/V_l #Bromoform concentration, g/L
      
      # ----- Parameters of the model -----

      #Microbial parameters 
      k_d<-8.3333e-004 #Death cell rate constant for all microbial groups, 1/h 
      khyd_ndf<-Inputparams["khyd_ndf"] #Hydrolysis rate constant of cell wall carbohydrates, 1/h 
      khyd_nsc<-Inputparams["khyd_nsc"] #Hydrolysis rate constant of non-structural carbohydrates, 1/h  
      khyd_pro<-Inputparams["khyd_pro"] #Hydrolysis rate constant of proteins, 1/h 
      km_su<-Inputparams["km_su"] #Maximum specific utilization rate constant of sugars, mol/(mol h) 
      Ks_su <-Inputparams["Ks_su"] #Monod constant associated with the utilization of sugars, mol/L   
      Ysu<-0.16 #Microbial biomass yield factor of sugars utilizers, mol/mol 
      km_aa<-Inputparams["km_aa"] #Maximum specific utilization rate constant of amino acids, mol/(mol h)
      Ks_aa<-Inputparams["Ks_aa"] #Monod constant associated with the utilization of amino acids mol/L 
      Yaa<-0.31 #Microbial biomass yield factor of amino acids utilizers, mol/mol 
      km_h2<-Inputparams["km_h2"] #Maximum specific utilization rate constant of hydrogen, mol/(mol h)
      Ks_h2<-Inputparams["Ks_h2"] #Monod constant associated with the utilization of hydrogen, mol/L
      Yh2<-0.006 #Microbial biomass yield factor of hydrogen utilizers, mol/mol 
      kLa<-8.3333 #Liquid-gas transfer constant, 1/h
      fac_aa<-0.67 #Stoichiometric coefficient of acetate production from amino acids fermentation 
      fpr_aa<-0.062 #Stoichiometric coefficient of propionate production from amino acids fermentation 
      fbu_aa<-0.24 #Stoichiometric coefficient of butyrate production from amino acids fermentation  
      fh2_aa<-0.82 #Stoichiometric coefficient of hydrogen production from amino acids fermentation 
      fIC_aa<-0.88 #Stoichiometric coefficient of inorganic carbon production from amino acids fermentation 
      fch_x<-0.20 #Fraction of carbohydrates of the biomass, g/g
      fpro_x<-0.55 #Fraction of proteins of the biomass, g/g 
      k_br<-Inputparams["k_br"] #Kinetic rate constant of bromoform utilization, 1/h

      #Physicochemical parameters 
      Trumen<-312.15 #Temperature, K
      R<-8.314*1e-2 #Ideal gas constant, bar*L/(mol*K)
      pgas_h2o<-0.08274 #Partial pressure water vapour, bar
      
      deltaH0_KH_co2<--19410 #J
      deltaH0_KH_ch4<--14240 #J
      deltaH0_KH_h2<--4180 #J
      
      #Henrys constant  M/bar, at T = 25C (298.15K)
      KH_co2_s<-0.035
      KH_ch4_s<-0.0014
      KH_h2_s<- 7.8e-4 
      
      KH_co2<-KH_co2_s*exp(-(19410/(R*100))*(1/298.15-1/Trumen))
      KH_ch4<-KH_ch4_s*exp(-(14240/(R*100))*(1/298.15-1/Trumen))
      KH_h2<-KH_h2_s*exp(-(4180/(R*100))*(1/298.15-1/Trumen))
      
      #Equilibrium constants 
      deltaH0_Ka_w<-55900   
      deltaH0_Ka_co2<-7646
      deltaH0_Ka_nh4<-51965
      
      #Acid-base constants mol/L
      K_w<-exp(deltaH0_Ka_w/(R*100)*(1/298.15-1/Trumen))*1e-14
      K_a_ac<-10^(-4.76)
      K_a_bu<-10^(-4.82)
      K_a_pr<-10^(-4.88)
      K_a_vfa<-K_a_ac
      K_a_co2<-10^(-6.35)*exp(deltaH0_Ka_co2/(R*100)*(1/298.15-1/Trumen))
      K_a_nh4<- 10^(-9.25)*exp(deltaH0_Ka_nh4/(R*100)*(1/298.15-1/Trumen))
      
      #Molecular weights (g/mol)
      w_h2<-2
      w_nh3<-17   
      w_vfa<-60.05 
      w_ch4<-16.0425
      w_co2<-44.0095 
      
      w_mb<-113 
      
      N_mb<-1
      w_su<-180.16   
      w_aa<-134    
      
      N_aa<-2
      
      w_ac<-60.05 
      w_bu<-88.1051 
      w_pr<-74.08
      
      pH<-6.6         #Average pH 
      ionH<-10^(-pH) 
      
      s_co2<-s_IC-K_a_co2*s_IC/(K_a_co2+ionH) #Concentration of carbon dioxide in the liquid phase, mol/L
      
      #liquid-gas transfer
      pgas_h2<-sg_h2*R*Trumen 
      pgas_ch4<-sg_ch4*R*Trumen
      pgas_co2<-sg_co2*R*Trumen
      
      Ptot<-pgas_h2+pgas_ch4+pgas_co2+pgas_h2o
      Ptot<-max(1.01325,Ptot)
      
      #Liquid-gas transfer phenomena
      pgas_co2<-R*Trumen*s_co2 #Partial pressure of CO2, bars
      pgas_ch4<-R*Trumen*sg_ch4 #Partial pressure of CH4, bars
      pgas_h2<-R*Trumen*sg_h2 #Partial pressure of H2, bars
      pgas_h2<-max(0,pgas_h2) 
      pgas_co2<-max(0,pgas_co2) 
      pgas_ch4<-max(0,pgas_ch4) 
      
      rhoT_h2<-kLa*(s_h2  - KH_h2*pgas_h2) #Liquid-gas transfer rate of H2, mol/(L h)
      rhoT_co2<-kLa*(s_co2 - KH_co2*pgas_co2) #Liquid-gas transfer rate of CO2, mol/(L h)
      rhoT_ch4<-kLa*(s_ch4 - KH_ch4*pgas_ch4) #Liquid-gas transfer rate of CH4, mol/(L h)
      
      #Microbial fermentation 
      #First order kinetics of hydrolysis
      rho_ndf<-khyd_ndf*z_ndf #Hydrolysis rate of NDF crabohydrates, g/(L h)
      rho_nsc<-khyd_nsc*z_nsc #Hydrolysis rate of NSC crabohydrates, g/(L h)
      rho_pro<-khyd_pro*z_pro #Hydrolysis rate of proteins, g/(L h)
      
      #Stoichiometry of the rumen fermentation 
      #Glucose utilization  
      # C6H12O6 + 2H2O  -> 2CH3COOH +2CO2 + 4H2                   (R1)
      # 3C6H12O6        -> 2CH3COOH  + 4CH3CH2COOH + 2CO2 + 2H2O  (R2)
      # C6H12O6         -> CH3CH2CH2COOH + 2CO2 + 2H2             (R3)
      # 5C6H12O6 + 6NH3 -> 6C5H7O2N + 18H2O                       (R4)
      fsu<-1-(5/6)*Ysu #Fraction of glucose utilized in catabolism
      p3<-Inputparams["p3"]
      p4<-Inputparams["p4"]
      lambda_1<-p3-p4*pgas_h2 #Molar fraction of the sugars utilized via reaction 1 
      
      p5<-Inputparams["p5"]   
      p6<-Inputparams["p6"]
      lambda_2<-p5+p6*pgas_h2 #Molar fraction of the sugars utilized via reaction 2

      lambda_1<-max(0,lambda_1)
      lambda_2<-max(0,lambda_2)
      lambda_3<-1- (lambda_1 + lambda_2) #Molar fraction of the sugars utilized via reaction 3
      lambda_3<-max(0,lambda_3)
      
      Yac_su<-fsu*(2*lambda_1 + (2/3)*lambda_2) #Yield factor of acetate during sugars utilization, mol/mol                        
      Ypr_su<-fsu*((4/3)*lambda_2) #Yield factor of propionate during sugars utilization, mol/mol
      Ybu_su<-fsu*(1*lambda_3) #Yield factor of butyrate during sugars utilization, mol/mol
      Yh2_su<-fsu*(4*lambda_1 + 2*lambda_3) #Yield factor of hydrogen during sugars utilization, mol/mol
      YIC_su<-fsu*(2*lambda_1 + (2/3)*lambda_2 + 2*lambda_3) #Yield factor of inorganic carbon during sugars utilization, mol/mol
      YIN_su<--Ysu #Yield factor of inorganic nitrogen during sugars utilization, mol/mol
      
      #Amino acids utilization
      Yac_aa<-(1-Yaa)*fac_aa #Yield factor of acetate during amino acids utilization, mol/mol
      Ypr_aa<-(1-Yaa)*fpr_aa #Yield factor of propionate during amino acids utilization, mol/mol
      Ybu_aa<-(1-Yaa)*fbu_aa #Yield factor of butyrate during amino acids utilization, mol/mol
      Yh2_aa<-(1-Yaa)*fh2_aa #Yield factor of hydrogen during amino acids utilization, mol/mol
      YIC_aa<-(1-Yaa)*fIC_aa #Yield factor of inorganic carbon during amino acids utilization, mol/mol
      YIN_aa<-N_aa -Yaa*N_mb #Yield factor of inorganic nitrogen during amino acids utilization, mol/mol
      
      #H2 utilization 
      #4H2 + CO2         -> CH4 + 2H2O      (R5)
      #10H2 + 5CO2 + NH3 -> C5H7O2N + 8H2O  (R6)
      fh2<-1-10*Yh2 #Fraction of H2 utilized in catabolism
      Ych4_h2<-fh2*(1/4) #Yield factor of methane during hydrogen utilization, mol/mol
      YIC_h2<--((1/4)*fh2 + 5*Yh2) #Yield factor of inorganic carbon during hydrogen utilization, mol/mol
      YIN_h2<--Yh2 #Yield factor of inorganic nitrogen during hydrogen utilization, mol/mol
      
      #Microbial kinetic rates
      #Nitrogen limitation 
      K_S_IN<-2.0e-4 #Nitrogen limitation constant, mol/L      
      I_IN_lim<-1/(1 + K_S_IN/s_IN) #Nitrogen limitation factor
      
      I_H2<-1.0 #Hydrogen control factor for sugars utilization 
      
      p1<-Inputparams["p1"] 
      p2<-Inputparams["p2"]
      I_br<-1-(1/(1+exp(-p1*(s_br+p2)))) #Inhibition factor of the methanogens growth rate by the action of bromoform 
      
      rho_su<-km_su*s_su*x_su*I_H2*I_IN_lim/(Ks_su + s_su) #Utilization rate of sugars with H2 regulation, mol/(L h)
      rho_aa<-km_aa*s_aa*x_aa/(Ks_aa + s_aa) #Utilization rate of amino acids, mol/(L h)
      rho_h2<-I_br*km_h2*s_h2*x_h2*I_IN_lim/(Ks_h2 + s_h2) #Utilization rate of hydrogen with bromoform inhibition, mol/(L h) 
      rho_xsu<-k_d*x_su #Cell death rate of sugars utilizers, mol/(L h)
      rho_xaa<-k_d*x_aa #Cell death rate of amino acids utilizers, mol/(L h)
      rho_xh2<-k_d*x_h2 #Cell death rate of hydrogen utilizers, mol/(L h)
      
      #Input fluxes (the feed intake)
      F_ndf_l_in<-wNDF*DMI/(100*V_l) #g/L h, input of NDF
      F_nsc_l_in<-wNSC*DMI/(100*V_l) #g/L h, input of NSC
      F_pro_l_in<-wPRO*DMI/(100*V_l) #g/L h, input of Proteins
      F_s_br_in<-wBr*DMI/(100*V_l) #g/L h, input of Bromoform
      
      #Output fluxes in the liquid phase
      D<-0.035 #1/h. Dilution rate for small particles from Bayat et al 2011
      F_ndf_out<-D*z_ndf #g/L h
      F_nsc_out<-D*z_nsc #g/L h
      F_pro_out<-D*z_pro #g/L h
      F_su_out<-D*s_su #mol/L h
      F_aa_out<-D*s_aa #mol/L h
      F_ac_out<-D*s_ac #mol/L h
      F_bu_out<-D*s_bu #mol/L h
      F_pr_out<-D*s_pr #mol/L h
      F_IN_out<-D*s_IN #mol/L h
      F_IC_out<-D*s_IC #mol/L h
      F_h2_out<-D*s_h2 #mol/L h
      F_ch4_out<-D*s_ch4 #mol/L h
      F_xsu_out<-D*x_su #mol/L h
      F_xaa_out<-D*x_aa #mol/L h
      F_xh2_out<-D*x_h2 #mol/L h
      F_s_br_out<-D*s_br #g/L h
      
      #Output fluxes in the gas phase
      q_g<-R*Trumen*V_l*(rhoT_h2 + rhoT_co2 + rhoT_ch4 )/(Ptot-pgas_h2o) #output flow of gas phase, L/h
      q_g<-max(q_g,0)
      
      F_co2g_out<-q_g*sg_co2/V_g #mol/L h
      F_h2g_out <-q_g*sg_h2/V_g  #mol/L h
      F_ch4g_out<-q_g*sg_ch4/V_g #mol/L h 
      
      # ---- State equations -----

      Z_ndf<-(F_ndf_l_in - rho_ndf -F_ndf_out)*V_l #Neutral detergent fiber (NDF), g  
      Z_nsc<-(F_nsc_l_in - rho_nsc + (fch_x*w_mb)*(rho_xsu + rho_xaa + rho_xh2) - F_nsc_out)*V_l #Non-structural carbohydrates (NSC), g 
      Z_pro<-(F_pro_l_in-rho_pro + (fpro_x*w_mb)*(rho_xsu + rho_xaa + rho_xh2) - F_pro_out)*V_l #Proteins, g
      
      S_su<-(rho_ndf/w_su + rho_nsc/w_su  - rho_su - F_su_out)*V_l #Sugars, mol 
      S_aa<-(rho_pro/w_aa - rho_aa - F_aa_out)*V_l #Average amino acids, mol 
      
      S_ac<-(Yac_su*rho_su + Yac_aa*rho_aa - F_ac_out)*V_l #Total acetate, mol 
      S_bu<-(Ybu_su*rho_su + Ybu_aa*rho_aa - F_bu_out)*V_l #Total butyrate, mol 
      S_pr<-(Ypr_su*rho_su + Ypr_aa*rho_aa - F_pr_out)*V_l #Total propionate, mol 
      
      S_IN<-(YIN_su*rho_su + YIN_aa*rho_aa + YIN_h2*rho_h2 - F_IN_out)*V_l #Inorganic nitrogen, mol 
      S_IC<-(YIC_su*rho_su + YIC_aa*rho_aa + YIC_h2*rho_h2 - rhoT_co2 - F_IC_out)*V_l #Inorganic carbon, mol 
      S_h2<-(Yh2_su*rho_su + Yh2_aa*rho_aa - rho_h2 - rhoT_h2 - F_h2_out)*V_l #Hydrogen in the liquid phase, mol 
      S_ch4<-(Ych4_h2*rho_h2 - rhoT_ch4 - F_ch4_out)*V_l #Methane in the liquid phase, mol 
      
      X_su<-(Ysu*rho_su - rho_xsu - F_xsu_out)*V_l #Sugars-utilizing microbes, mol 
      X_aa<-(Yaa*rho_aa - rho_xaa - F_xaa_out)*V_l #Amino acids-utilizing microbes, mol 
      X_h2<-(Yh2*rho_h2 - rho_xh2 - F_xh2_out)*V_l #Hydrogen-utilizing microbes (methanogens), mol 
      
      Sg_co2<-(V_l*rhoT_co2/V_g - F_co2g_out)*V_g #Carbon dioxide in the gas phase, mol 
      Sg_h2<-(V_l*rhoT_h2/V_g  - F_h2g_out)*V_g #Hydrogen in the gas phase, mol 
      Sg_ch4<-(V_l*rhoT_ch4/V_g - F_ch4g_out)*V_g #Methane in the gas phase, mol 
      S_br<-(F_s_br_in - k_br*s_br - F_s_br_out)*V_l #Bromoform, g 

      return (list(c(Z_ndf,Z_nsc,Z_pro,S_su,S_aa,S_ac,S_bu,S_pr,S_IN,S_IC,S_h2,S_ch4,X_su,X_aa,X_h2,Sg_co2,Sg_h2,Sg_ch4,S_br),
                   pgas_h2=pgas_h2,lambda_1=lambda_1,lambda_2=lambda_2,I_H2=I_H2,I_br=I_br,q_ch4g_out=q_g*sg_ch4))   
  })
}