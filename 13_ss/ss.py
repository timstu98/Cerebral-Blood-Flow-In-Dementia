import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path

folder = Path.cwd().parent

baselines = pd.read_csv(folder / '_storage/baselines' / 'baselines_met.csv', index_col=0)

H = 0.42 #no units,ratio #Hematocrit assumed to be constant
# M_o = 30* 10e-9 #mol_ O2/(mL/s) #Taken from table 2 from Wiley Payne paper
cHb = 0.2 #mL_O2/m #Taken from table 2 from Wiley Payne paper
K_o = 5e-8 #µL/(mm*s*mmHg) #payne paper and boas et al

#convert param to my SU
M_constant = 8.20e-4 # cm3_O2/(cm3*s) , still unsure about the exact conversion so will just input in this section
K = K_o * (10**6) # µm2/(s*mmHg)

hill_constant = 2 # mmHg

############

its = 1000
I = [i for i in range(its)]

def create_df():
    name = pd.DataFrame(np.NaN, index=np.arange(13), columns=I)
    return name
    
M = create_df()
M[0] = [M_constant for i in range(13)]
Sat_in = create_df()
Sat_in.iloc[0,:] = 0.94
Sat_out = create_df()
p_diff = create_df() # pb-pt
Sat_ave = create_df()
pb = create_df()
pt = create_df()

flow = pd.DataFrame(np.zeros(its))



for i in range(its):
    if i!=0:
        M[i] = M_constant*pt[i-1]/(pt[i-1]+hill_constant)
    p_diff[i] = M[i]*ss['wall thickness(µm)']*ss['Vt(µm3)']/(2*np.pi*K*0.5*ss['Diameter(µm)']*ss['Length(µm)'])
    for j in range(len(ss)):
        if j != 0:
            Sat_in.iloc[j,i] = Sat_out.iloc[j-1,i]

        Sat_out.iloc[j,i] = Sat_in.iloc[j,i] - ( (2*np.pi*K*0.5*ss.at[j,'Diameter(µm)']*ss.at[j,'Length(µm)']) / (ss.at[j,'wall thickness(µm)']*cHb*H) ) * p_diff.iloc[j,i] / ss.at[j,'Q(µm3/s)']

    Sat_ave[i] = ( Sat_in[i] + Sat_out[i] ) / 2

    for j in range(len(Sat_ave)):
            root_temp = None
            sols = 0
            root_temp = np.roots([1,0,150,23400*Sat_ave.iloc[j,i]/(Sat_ave.iloc[j,i]-1)]) #gives mmHg , *133.322
            sols = check_imag_roots_real(root_temp) 
            #sols = sols*133.322 ###Remove this step if want in mmHg
            pb.iloc[j,i] = sols

    pt[i] = pb[i] - p_diff[i]
    
    

# display(ss)