### import modules
print('\n' * 1,'############ NEW RUN ############','\n' * 1)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
from tqdm import tqdm
from pathlib import Path
import datetime
from pydub import AudioSegment
from pydub.playback import play
from parameters import p

### Save_info
save = False
print('Save:', save)
save_changes = True
save_vessels = True
save_combined = True
folder = Path.cwd().parent

### import required files
baselines = pd.read_csv(folder / '_storage/baselines' / 'baselines_met.csv', index_col=0)
baseline = baselines.drop(['Viscocity(mmHg*s)','Resistance for U', 'U in single(µm/s)'],axis=1).rename(columns={'Saturation ave': 'Sav'})
constant = baseline.copy().drop(['Q in single(µm3/s)','Pressure Drop(mmHg)','Saturation in','Saturation out','Resistance for Q','dS','Sav','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)'], axis=1)
print('Baselines read in.')

### parameters
#most already loaded in under p class
alpha = 0.2
delay = 1

max_time = 80
no = 1001

delay_pressure_drop = 2
time_for_drop = 2
ratio_drop = 0.5

kp_constant = 1/300 # s-1
kn_constant = 1/600 # s-1

R = 1
k = 0.1
tau_c = 600 #tau_c = 1/600 #s-1

### Create Arrays

t = pd.Series(np.linspace(0,max_time,no))
dt = max_time/(no-1)
print('Time array done.')
drop_value = 60 - 34.18*ratio_drop
pressure_in = pd.Series(np.zeros(len(t)))
for i in range(len(t)):
    # if i % (no/4) < 0.01 or i - no/20 < 0.01:
    #     print(100*i/no,'percent done')
    if t[i] <= delay_pressure_drop:
        pressure_in[i] = 60
    elif t[i] <= delay_pressure_drop + time_for_drop:
        pressure_in[i] = 60 - ((t[i] - delay_pressure_drop)/time_for_drop)*(60-drop_value)
    else:
        pressure_in[i] = drop_value
pressure_out = pd.Series(np.zeros(len(t))) + 60 - 34.18
pressure_difference = pressure_in - pressure_out
print('Pressure array done.')

################################## Functions ##################################
#Used in main script
def network_initial_conditions(network):
    network_row0 = network.iloc[0,:].copy()

    network_row0['phi'] = 1
    network_row0['phi_min'] = 1
    network_row0['Ap'] = 1
    network_row0['Dp'] = 0
    network_row0['An'] = 1
    network_row0['Dn'] = 0
    network_row0['R_tot'] = total_R(baseline,network_row0['phi'],alpha)
    network_row0['Q_tot'] = network_row0['pressure_difference'] / network_row0['R_tot']
    network_row0['Q_norm'] = network_row0['Q_tot'] / baseline_total_flow
    network_row0['pt_averaged'] = baseline_pt_averaged
    network_row0['kp_p'] = (kp_constant/2)*(1-np.tanh( (network_row0['pt_averaged']-pt_averaged_50) / pt_averaged_5))
    network_row0['kp_n'] = (kp_constant/2)*(1+np.tanh( (network_row0['pt_averaged']-pt_averaged_50) / pt_averaged_5 ))
    network_row0['c'] = c_baseline 
    network_row0['kn_p'] = (kn_constant/2)*(1+np.tanh( (network_row0['c']-c_50) / c_5))
    network_row0['kn_n'] = (kn_constant/2)*(1-np.tanh( (network_row0['c']-c_50) / c_5 ))

    network.iloc[0,:] = network_row0.copy()
    return network
def total_R(baselines,phi,alpha):
    C_ = alpha*(baselines.at[6,'Resistance for Q']/2)/phi**4  + (1-alpha)*(baselines.at[6,'Resistance for Q']/2)
    C_6 = (C_ + baselines.at[5,'Resistance for Q'] + baselines.at[7,'Resistance for Q'])/2
    C_65 = (C_6 + baselines.at[4,'Resistance for Q'] + baselines.at[8,'Resistance for Q'])/2
    C_654 = (C_65 + baselines.at[3,'Resistance for Q'] + baselines.at[9,'Resistance for Q'])/2
    C_6543 = (C_654 + baselines.at[2,'Resistance for Q'] + baselines.at[10,'Resistance for Q'])/2
    C_65432 = (C_6543 + baselines.at[1,'Resistance for Q'] + baselines.at[11,'Resistance for Q'])/2
    C_654321 = C_65432 + baselines.at[0,'Resistance for Q'] + baselines.at[12,'Resistance for Q']
    R_total = C_654321
    return R_total

#Used in master and after
def master(alpha,delay):
    #Set up dicts to save to:
    network, all_vessels, all_changes, all_ks = {}
    print(network)



    for i in tqdm(range(no-1)):
        do_nothing=0


    return network, all_vessels, all_changes, all_ks
################################## Model Run ##################################

# Initialise Vessels and Networks.
first_iteration_vessels = baseline.copy().loc[:,['Name','Number','Diameter(µm)','Length(µm)','Saturation in','Saturation out','Sav','Resistance for Q','Q in single(µm3/s)','Vt(µm3)','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)']]

network_column_names = ['t', 'pressure_difference', 'phi','phi_min','R_tot','Q_tot','Q_norm','pt_averaged','kp_p','kp_n','Ap','Dp','kn_p','kn_n','An','Dn','c']
network = pd.DataFrame(columns = network_column_names)
network['t'] = t
network['pressure_difference'] = pressure_difference

# Save useful values
delay_number = round(delay/dt)
no_vessels = len(first_iteration_vessels)
baseline_total_flow = first_iteration_vessels['Q in single(µm3/s)'][0]

#Further model params
baseline_pt_averaged = ( (first_iteration_vessels['tissue partials(mmHg)']*first_iteration_vessels['Vt(µm3)']).sum() ) / (first_iteration_vessels['Vt(µm3)'].sum())
pt_averaged_50 = baseline_pt_averaged*0.5
pt_averaged_5 = baseline_pt_averaged*0.1  

c_baseline = R/(1+k)
c_50 = c_baseline*10
c_5 = c_baseline

# Initial Conditions
network = network_initial_conditions(network)

#Run model
print('Running model.')
network, all_vessels, all_changes, all_ks = master(alpha,delay)

################################## SAVING ##################################
def combine_network_changes(network,RK4):
    network['dphidt'] = np.nan
    network['dApdt'] = np.nan
    network['dDpdt'] = np.nan
    network['dAndt'] = np.nan
    network['dDndt'] = np.nan
    network['dcdt'] = np.nan
    for i in tqdm(range(len(network)-5)):
        network.at[i,'dphidt'] = RK4[i][4]
        network.at[i,'dApdt'] = RK4[i][0]
        network.at[i,'dDpdt'] = RK4[i][1]
        network.at[i,'dAndt'] = RK4[i][2]
        network.at[i,'dDndt'] = RK4[i][3]
        network.at[i,'dcdt'] = RK4[i][5]
    network.to_csv(main_path / 'combined_network.csv')

def combine_vessels_changes(vessels_alltime,RK4):
    for i in tqdm(range(len(vessels_alltime)-5)):
        vessels_alltime[i]['dptdt'] = RK4[i][6]['dptdt']
        vessels_alltime[i]['dSoutdt'] = RK4[i][7]['dSoutdt']
        vessels_alltime[i]['dSindt'] = RK4[i][8]['dSindt']
    save = pd.concat(vessels_alltime)
    save.to_csv(main_path / 'combined_vessels.csv')

if save == True:
    print('Saving ...')
    txt_file = [
                'Time of simulation: ' + str(max_time) + ' seconds\n' ,
                'Number of iterations: ' + str(no) + '\n',
                'therefore time step is: ' + str(dt) + ' seconds\n' ,
                'Delay for pressure drop: ' + str(delay_pressure_drop) + ' seconds\n' ,
                'Time over which pressure drops: ' + str(time_for_drop) + ' seconds\n',
                'Ratio pressure drop: ' + str(ratio_drop) + '\n',
                'kp_constant: ' + str(kp_constant) + '\n',
                'kn_constant: ' + str(kn_constant) + '\n',
                'Date and time: ' + str(datetime.datetime.now()) + '\n',
                'save_changes: ' + str(save_changes) + '\n'
                ]

    count = 1
    while count < 1000:
        if len(str(count)) == 3:
            subfoulder = str(count)
        elif len(str(count)) == 2:
            subfoulder = str(0) + str(count)
        else:
            subfoulder = str(0) + str(0) + str(count)
        attempt = folder / '_storage' / 'main' / subfoulder
        if attempt.exists():
            count = count+1
        else:
            main_path = attempt
            main_path.mkdir()
            break

    print('save number: ',count)
    file= open(main_path / "info.txt","w") 
    file.writelines(txt_file)
    file.close()  
    print('text file saved')

    if save_combined == True:
        combine_network_changes(N,RK4_change_alltime)
        print('combined network saved')
        combine_vessels_changes(V,RK4_change_alltime)
        print('combined vessels saved')
        plt.plot(N['t'],N['Q_norm'])
        plt.ylabel('Normalised Flow')
        plt.xlabel('Time / seconds')
        plt.grid(which='both')
        plt.savefig(main_path / 'Q_vs_t')
        print('ALL DONE.')
        pickle_out3 = open(main_path / 'ks.pickle' , 'wb' )
        pickle.dump(ks_alltime, pickle_out3)
        pickle_out3.close()
        print('ks saved')
    else:
        N.to_csv(main_path / 'network.csv')
        print('Network saved')
        if save_vessels == True:
            pickle_out4 = open(main_path / 'vessels.pickle' , 'wb' )
            pickle.dump(V, pickle_out4)
            pickle_out4.close()
            print('Vessels saved')
        if save_changes == True:
            pickle_out3 = open(main_path / 'changes.pickle' , 'wb' )
            pickle.dump(RK4_change_alltime, pickle_out3)
            pickle_out3.close()
            print('changes saved')
        plt.plot(N['t'],N['Q_norm'])
        plt.ylabel('Normalised Flow')
        plt.xlabel('Time / seconds')
        plt.grid(which='both')
        plt.savefig(main_path / 'Q_vs_t')
        print('ALL DONE.')
sound = AudioSegment.from_mp3('minecraft.mp3')
play(sound)
print('Finished.')
    

# ks_names = pd.array(['0: dApdt','1: dDpdt','2: dAndt','3: dDndt','4: dphidt','5: dcdt','6:v dptdt','7v: dSoudt','8v: dSindt'])
# values_2_change_names = pd.array(['0: Ap','1: Dp','2: An','3: Dn','4: phi','5: c','6:v tissue partials(mmHg)','7v: Saturation out','8v: Saturation in'])