### import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
from tqdm import tqdm
from pathlib import Path
import datetime
folder = Path.cwd().parent
save = True

### import required files
baselines = pd.read_csv(folder / '_storage/baselines' / 'baselines_met.csv', index_col=0)

### set parameters
#outside params, _o stands for original
H = 0.42 #no units,ratio #Hematocrit assumed to be constant
cHb = 0.2 #mL_O2/mL #Taken from table 2 from Wiley Payne paper
paO2_bar_t = 15 #mmHG #Taken from table 2 from Wiley Payne paper
K_o = 5e-8 #µL/(mm*s*mmHg) #payne paper and boas et al
alpha_t = 2.6e-5 #mL_O2/(mL*mmHg) from payne paper, solutbility coeff of oxygen in brain tissue

#convert param to my SU
M_constant = 8.20e-4 # cm3_O2/(cm3*s) , still unsure about the exact conversion so will just input in this section
K = K_o * (1e6) # µm2/(s*mmHg)

#model input params
delay_pressure_drop = 4
time_for_drop = 1
optimised_value_tau = 2.292929292929293
phi_min_baseline = 0.153
n=1
ratio_drop = 0.5
hill_constant = 2 # mmHg
kp_constant = 1/300 # s-1
kn_constant = 1/600 # s-1

#loop info
# PREVIOUS: no = 400001 max_time = 300
no = 40001
max_time = 30

### set the arrays
t = pd.Series(np.linspace(0,max_time,no))
dt = max_time/(no-1)
drop_value = 60 - 34.18*ratio_drop

print('Values initialised.')

# PRessure stuff
pressure_in = pd.Series(np.zeros(len(t)))
for i in range(len(t)):
    if t[i] <= delay_pressure_drop:
        pressure_in[i] = 60
    elif t[i] <= delay_pressure_drop + time_for_drop:
        pressure_in[i] = 60 - ((t[i] - delay_pressure_drop)/time_for_drop)*(60-drop_value)
    else:
        pressure_in[i] = drop_value
pressure_out = pd.Series(np.zeros(len(t))) + 60 - 34.18
pressure_difference = pressure_in - pressure_out

print('Pressure arrays finished.')

def total_R(baselines,phi,alpha):
    """Calculates the resistance of the network.
    
    """
    
    C_ = alpha*(baselines.at[6,'Resistance for Q']/2)/phi**4  + (1-alpha)*(baselines.at[6,'Resistance for Q']/2) #### IF THIS CORRECT?
    C_6 = (C_ + baselines.at[5,'Resistance for Q'] + baselines.at[7,'Resistance for Q'])/2
    C_65 = (C_6 + baselines.at[4,'Resistance for Q'] + baselines.at[8,'Resistance for Q'])/2
    C_654 = (C_65 + baselines.at[3,'Resistance for Q'] + baselines.at[9,'Resistance for Q'])/2
    C_6543 = (C_654 + baselines.at[2,'Resistance for Q'] + baselines.at[10,'Resistance for Q'])/2
    C_65432 = (C_6543 + baselines.at[1,'Resistance for Q'] + baselines.at[11,'Resistance for Q'])/2
    C_654321 = C_65432 + baselines.at[0,'Resistance for Q'] + baselines.at[12,'Resistance for Q']
    R_total = C_654321
    return R_total

if i % 0.5 < 0.01:
    

def fn_dphidt(i,phi,Q_norm,phi_min,delay_number): #fn on phi delayed but not t? how do i get this involved?
    R_tot = total_R(baselines,phi[i],alpha)
    Q_tot = pressure_difference[i]/R_tot
    Q_norm = Q_tot / baselines['Q in single(µm3/s)'][0]
    dphidt = (1/optimised_value_tau)*( -phi[i-delay_number] + Q_norm*(1-phi_min[i]) + phi_min[i] )
    return dphidt

def euler(i,phi,Q_norm,phi_min,delay_number):
    k1 = fn_dphidt(i,phi,Q_norm,phi_min,delay_number)
    k2 = fn_dphidt()
    k3 = 0
    k4 = 0
    return k1

# def RK4():
#     k1 = dydt
#     k2 = 0
#     k3 = 0
#     k4 = 0
#     return dydt

# def RK4(y,t,dt,dydt[i]):
#     k1 = dydt(i)
#     # l1 = 
#     # m1 = 

def run(alpha,delay,baselines):
    delay_number = round(delay/dt)
    no_vessels = 13
    
    R_tot = pd.Series(np.zeros(len(t)))
    Q_tot = pd.Series(np.zeros(len(t)))
    Q_norm = pd.Series(np.ones(len(t)))
    
    phi = pd.Series(np.ones(no))
    phi_min = pd.Series(np.ones(no))
    dphidt = pd.Series(np.zeros(no))

    
    
    for i in tqdm(range(len(t))):
        if i >=delay_number: 
            phi[i+1] = phi[i] + euler(i,phi,Q_norm,phi_min,delay_number)*dt
            
            R_tot[i] = total_R(baselines,phi[i],alpha)
            Q_tot[i] = pressure_difference[i]/R_tot[i]
            Q_norm[i] = Q_tot[i] / baselines['Q in single(µm3/s)'][0]
            phi_min[i+1] = 1-(alpha**0.25)*(1-phi_min_baseline)*(1-Q_norm[i])**n
        else:
            ###############################################################################
            #URGENTLY NEED TO FIX THIS BIT UP PROPERLY#
            ###############################################################################
            Q_tot[i] = 3.713608e+06
#             Q_tot[i] = pressure_difference[0]/R_tot[0] # Havent tried putting this into model but dont see why it wouldnt work
            
    combined = pd.DataFrame()
    combined['t'] = t
    combined['pressure_difference'] = pressure_difference
    combined['phi'] = phi
    combined['phi_min'] = phi_min
    combined['dphidt'] = dphidt
    combined['R_tot'] = R_tot
    combined['Q_tot'] = Q_tot
    combined['Q_norm'] = Q_norm
            
    return combined

alpha = 0.2
delay = 2.99
out = run(alpha,delay,baselines)
out.to_csv('look_at.csv')

plt.plot(t,out['Q_norm'])
plt.show()

display(combined)