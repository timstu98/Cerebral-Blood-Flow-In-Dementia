class p:
    #outside params, _o stands for original 
    H = 0.42 #no units,ratio #Hematocrit assumed to be constant
    cHb = 0.2 #mL_O2/mL #Taken from table 2 from Wiley Payne paper
    paO2_bar_t = 15 #mmHG #Taken from table 2 from Wiley Payne paper
    K_o = 5e-8 #µL/(mm*s*mmHg) #payne paper and boas et al
    # alpha_t = (2.6e-5)**-1 #mL_O2/(mL*mmHg) from payne paper, solutbility coeff of oxygen in brain tissue
    alpha_t = 2.6e-5 #mL_O2/(mL*mmHg) from payne paper, solutbility coeff of oxygen in brain tissue

    #convert param to my SU
    M_constant = 8.20e-4 # cm3_O2/(cm3*s) , still unsure about the exact conversion so will just input in this section
    K = K_o * (1e6) # µm2/(s*mmHg)

    #model input params
    hill_constant = 2 # mmHg
    optimised_value_tau = 2.292929292929293*60
    phi_min_baseline = 0.153
    n=1