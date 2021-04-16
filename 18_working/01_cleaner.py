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
save = True
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
alpha = 0.4
delay = 2

max_time = 600
no = 225001

delay_pressure_drop = 2
time_for_drop = 2
ratio_drop = 0.3

kp_constant = 1/300 # s-1
kn_constant = 1/600 # s-1

R = 1
k = 0.1
tau_c = 1/600 #s-1

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

    network_row0['t'] = t[0]
    network_row0['pressure_difference'] = pressure_difference[0]
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
def check_imag_roots_real(two_roots):
    """Takes two imaginary roots from the Servinghaus equation and returns the real one.
    
    Will print if there is an error.
    Returns real root."""
    
    counter = 0
    sols = 0
    for i in range(len(two_roots)):
        if two_roots[i].imag <= 10e-13 and two_roots[i].imag >= -10e-13:
            sols = two_roots[i].real
            counter += 1
    if counter != 1:
        print('Something seems wrong with the partial presure to Saturation cubic solver as not one solution is returned.')

    return sols

def interpolate(array,i):
    if i % 1 > 0.01:
        value = ( array[int(i+0.5)] + array[int(i-0.5)] ) / 2
    else:
        value = array[i]
    return value

def update_v_n(the_change, multiplier, vessels, network_row):

    # This is flawed as it doesnt update all the other values for these new values and save properly. everything else will be fucked up for a while
    network_row_new = network_row.copy()

    # print(the_change[0])
    network_row_new['Ap'] = network_row['Ap'] + multiplier* the_change[0]*dt
    network_row_new['Dp'] = network_row['Dp'] + multiplier* the_change[1]*dt
    network_row_new['An'] = network_row['An'] + multiplier* the_change[2]*dt
    network_row_new['Dn'] = network_row['Dn'] + multiplier* the_change[3]*dt
    network_row_new['phi'] = network_row['phi'] +multiplier* the_change[4]*dt
    network_row_new['c'] = network_row['c'] + multiplier* the_change[5]*dt

    vessels['tissue partials(mmHg)'] = vessels['tissue partials(mmHg)'] + multiplier* the_change[6]['dptdt']*dt
    vessels['Saturation out'] = vessels['Saturation out'] + multiplier* the_change[7]['dSoutdt']*dt
    vessels['Saturation in'] = vessels['Saturation in'] + multiplier* the_change[8]['dSindt']*dt

    return vessels, network_row_new

def combine_ks(k1s,k2s,k3s,k4s):

    out = {}
    out[6] = empty_column('dptdt')
    out[7] = empty_column('dSoutdt')
    out[8] = empty_column('dSindt')

    for n in range(len(k1s)):
        if 0<=n<=5:
            out[n] = ( k1s[n] + k2s[n]*2 + k3s[n]*2 + k4s[n])/6
        else:
            out[n].iloc[:,1] = ( k1s[n].iloc[:,1] + k2s[n].iloc[:,1]*2 + k3s[n].iloc[:,1]*2 + k4s[n].iloc[:,1])/6

    return out

def empty_column(string_name):
    df = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'] })
    df[string_name] = 0 #np.nan
    return df

def RK4(i,network_row_main_updated,vessels_main_updated):

    network_row_in = network_row_main_updated.copy()
    vessels_in = vessels_main_updated.copy()

    # no need to update
    k1s , network_row_all_updated, vessels_all_updated = single_RK_run( i , network_row_in , vessels_in)

    vessels , network_row = update_v_n(k1s,0.5,vessels_in,network_row_in)
    k2s , _, _ = single_RK_run(i+0.5,network_row,vessels)

    vessels,network_row = update_v_n(k2s,0.5,vessels_in,network_row_in)
    k3s, _, _ = single_RK_run(i+0.5,network_row,vessels)

    vessels,network_row = update_v_n(k3s,1,vessels_in,network_row_in)
    k4s, _, _ = single_RK_run(i+1,network_row,vessels)

    all_ks = {}
    all_ks[0] = k1s
    all_ks[1] = k2s
    all_ks[2] = k3s
    all_ks[3] = k4s

    change = combine_ks(k1s,k2s,k3s,k4s)

    return change , network_row_all_updated , vessels_all_updated, all_ks

def single_RK_run(i,network_row,vessels):

    ks = {} #change empty column back
    ks[0] = 0
    ks[1] = 0
    ks[2] = 0
    ks[3] = 0
    ks[4] = 0 
    ks[5] = 0
    ks[6] = empty_column('dptdt')
    ks[7] = empty_column('dSoutdt')
    ks[8] = empty_column('dSindt')

    network_row['R_tot'] = total_R(baselines,network_row['phi'],alpha)
    network_row['Q_tot'] = interpolate(network['pressure_difference'],i) / network_row['R_tot']
    network_row['Q_norm'] = network_row['Q_tot'] / baseline_total_flow
    network_row['phi_min'] = 1-(alpha**0.25)*(1-p.phi_min_baseline)*(1-network_row['Q_norm'])**p.n

    vessels['pt_by_Vt'] = vessels['tissue partials(mmHg)'] * vessels['Vt(µm3)']
    network_row['pt_volume_averaged'] = vessels['pt_by_Vt'].sum()/13 
    
    network_row['kp_p'] = (kp_constant/2)*(1-np.tanh( (network_row['pt_averaged']-pt_averaged_50) / pt_averaged_5))
    network_row['kp_n'] = (kp_constant/2)*(1+np.tanh( (network_row['pt_averaged']-pt_averaged_50) / pt_averaged_5 ))
    ks[0]= -network_row['kp_p']*network_row['Ap'] + network_row['kp_n']*(1-network_row['Ap']-network_row['Dp']) #dApdt
    ks[1] = network_row['kp_p']*(1-network_row['Ap']-network_row['Dp']) #dDpdt
    
    network_row['kn_p'] = (kn_constant/2)*(1+np.tanh( (network_row['c']-c_50) / c_5))
    network_row['kn_n'] = (kn_constant/2)*(1-np.tanh( (network_row['c']-c_50) / c_5 ))
    ks[2] = -network_row['kn_p']*network_row['An'] + network_row['kn_n']*(1-network_row['An']-network_row['Dn']) #dAndt
    ks[3] = network_row['kn_p']*(1-network_row['An']-network_row['Dn']) #dDndt
    
    ks[5] = (1/tau_c)*(R- k*network_row['c'] - network_row['Q_norm']*network_row['c'] )

    if i-delay_number <0: # As cant index before 0 as not initialised for that.
        ks[4] =  (1/p.optimised_value_tau)*( -                    1                       + network_row['Ap']*(network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
    else:
        ks[4] =  (1/p.optimised_value_tau)*( - interpolate(network['phi'],i-delay_number) + network_row['Ap']*(network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
     
    #Update other values for new flow
    vessels.at[6,'Diameter(µm)'] = constant.at[6,'Diameter(µm)']*network_row['phi']    
    vessels.at[6,'X-Area wall(???)'] = np.pi * ( 0.16*(vessels.at[6,'Diameter(µm)'])**2 + 1.4*(vessels.at[6,'Diameter(µm)']) + 14 )
    vessels.at[6,'wall thickness(µm)'] = ( -(vessels.at[6,'Diameter(µm)']) + np.sqrt((vessels.at[6,'Diameter(µm)'])**2 + (4*(vessels.at[6,'X-Area wall(???)'])/np.pi))) / 2
    vessels['Sav'] = (vessels['Saturation in']+vessels['Saturation out'])/2
    vessels['M'] = (network_row['An'])*p.M_constant*vessels['tissue partials(mmHg)']/(vessels['tissue partials(mmHg)']+p.hill_constant)

    for j in range(no_vessels):
        root_temp = None
        sols = 0
        root_temp = np.roots([1,0,150,23400*vessels.at[j,'Sav']/(vessels.at[j,'Sav']-1)]) #gives mmHg
        sols = check_imag_roots_real(root_temp) 
        vessels.at[j,'partial pressure blood(mmHg)'] = sols    

    vessels['Q in single(µm3/s)'] = first_iteration_vessels['Q in single(µm3/s)']*network_row['Q_norm']
    
    for j in range(len(vessels)):
        """For each vessel
        Update Q in single based on updated Q_tot
        Calculate dSoutdt and make this dSindt for next vessel."""
        if j ==0:
            ks[8].at[j,'dSindt'] = 0
        else:
            ks[8].at[j,'dSindt'] = ks[7].at[j-1,'dSoutdt']

        ks[7].at[j,'dSoutdt'] = 2/(np.pi*((vessels.at[j,'Diameter(µm)']**2)/4)*vessels.at[j,'Length(µm)']) \
                                * ( - (2 * np.pi * p.K * (vessels.at[j,'Diameter(µm)']/2) * constant.at[j,'Length(µm)'])/(vessels.at[j,'wall thickness(µm)']*p.cHb*p.H) \
                                * (vessels.at[j,'partial pressure blood(mmHg)']-vessels.at[j,'tissue partials(mmHg)'])\
                                - vessels.at[j,'Q in single(µm3/s)']*(vessels.at[j,'Saturation out']-vessels.at[j,'Saturation in']) )\
                                - ks[8].at[j,'dSindt']
    
    ks[6]['dptdt'] = (1/(p.alpha_t*constant['Vt(µm3)']))\
                    *(((2 * np.pi * p.K * (vessels['Diameter(µm)']/2) * constant['Length(µm)'])/(vessels['wall thickness(µm)']))\
                    *(vessels['partial pressure blood(mmHg)']-vessels['tissue partials(mmHg)'])\
                    - vessels['M'] * constant['Vt(µm3)'])

    return ks,network_row,vessels

def master(alpha,delay):
    #Set up dicts to save to:
    all_vessels, all_changes, all_ks = {},{},{}

    network_row_main_updated = network.iloc[0,:]
    vessels_main_updated = first_iteration_vessels

    for i in tqdm(range(no-1)):

        network_row_main_updated['t'] = t[i]
        network_row_main_updated['pressure_difference'] = pressure_difference[i]

        changes , network_row_all_updated, vessels_all_updated, ks = RK4(i, network_row_main_updated, vessels_main_updated)

        network.iloc[i,:] = network_row_all_updated
        all_vessels[i] = vessels_all_updated
        all_changes[i] = changes
        all_ks[i] = ks

        vessels_main_updated, network_row_main_updated = update_v_n(changes, 1, vessels_all_updated, network_row_all_updated)

    return network, all_vessels, all_changes, all_ks
################################## Model Run ##################################

# Initialise Vessels and Networks.
first_iteration_vessels = baseline.copy().loc[:,['Name','Number','Diameter(µm)','Length(µm)','Saturation in','Saturation out','Sav','Resistance for Q','Q in single(µm3/s)','Vt(µm3)','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)']]

network_column_names = ['t', 'pressure_difference', 'phi','phi_min','R_tot','Q_tot','Q_norm','pt_averaged','kp_p','kp_n','Ap','Dp','kn_p','kn_n','An','Dn','c']
network = pd.DataFrame(np.nan,index=np.arange(0,no,1).tolist(), columns = network_column_names)
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
# c_50 = c_baseline*10 ######## CHANGING ##########
# c_50 = c_baseline*2
c_50 = 1.1*c_baseline
c_5 = 0.06*c_baseline
# c_5 = c_baseline

# Initial Conditions
network = network_initial_conditions(network)

#Run model
print('Running model.')
# print(first_iteration_vessels)
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
    save = pd.concat(vessels_alltime,sort=False)
    save.to_csv(main_path / 'combined_vessels.csv')

if save == True:
    print('Saving ...')
    txt_file = [
                'Time of simulation: ' + str(max_time) + ' seconds\n' ,
                'Number of iterations: ' + str(no) + '\n',
                'therefore time step is: ' + str(dt) + ' seconds\n' ,
                'Delay for pressure drop: ' + str(delay_pressure_drop) + ' seconds\n' ,
                'Time over which pressure drops: ' + str(time_for_drop) + ' seconds\n',
                'Alpha: ' + str(alpha) + '\n',
                'delay: ' + str(delay) + '\n',
                'Ratio pressure drop: ' + str(ratio_drop) + '\n',
                'kp_constant: ' + str(kp_constant) + '\n',
                'kn_constant: ' + str(kn_constant) + '\n',
                'Date and time: ' + str(datetime.datetime.now()) + '\n',
                'R: ' + str(R) + '\n',
                'k: ' + str(k) + '\n',
                'tau_c: ' + str(tau_c) + '\n',
                'save_changes: ' + str(save_changes) + '\n'
                ]

    count = 1
    while count < 10000:
        if len(str(count)) == 4:
            subfoulder = str(count)
        elif len(str(count)) == 3:
            subfoulder = str(0) + str(count)
        elif len(str(count)) == 2:
            subfoulder = str(0) + str(0) + str(count)
        else:
            subfoulder = str(0) + str(0) + str(0) + str(count)
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

    combine_network_changes(network,all_changes)
    print('combined network saved')
    combine_vessels_changes(all_vessels,all_changes)
    print('combined vessels saved')
    fig = plt.figure()
    plt.plot(network['t'],network['Q_norm'])
    plt.ylabel('Normalised Flow')
    plt.xlabel('Time / seconds')
    plt.ylim(bottom=0)
    plt.grid(which='both')
    plt.savefig(main_path / 'Q_vs_t')
    fig = plt.figure()
    plt.plot(network['t'],network['c'])
    plt.ylabel('concentration')
    plt.xlabel('Time / seconds')
    plt.ylim(bottom=0)
    plt.grid(which='both')
    plt.savefig(main_path / 'c_vs_t')
    fig = plt.figure()
    plt.plot(network['t'],network['An'])
    plt.ylabel('An')
    plt.xlabel('Time / seconds')
    plt.ylim(bottom=0)
    plt.grid(which='both')
    plt.savefig(main_path / 'An_vs_t')
    print('Graphs saved')
    pickle_out3 = open(main_path / 'ks.pickle' , 'wb' )
    pickle.dump(all_ks, pickle_out3)
    pickle_out3.close()
    print('ks saved')
    print('ALL DONE.')

sound = AudioSegment.from_mp3('minecraft.mp3')
play(sound)
print('Finished.')
    

# ks_names = pd.array(['0: dApdt','1: dDpdt','2: dAndt','3: dDndt','4: dphidt','5: dcdt','6:v dptdt','7v: dSoudt','8v: dSindt'])
# values_2_change_names = pd.array(['0: Ap','1: Dp','2: An','3: Dn','4: phi','5: c','6:v tissue partials(mmHg)','7v: Saturation out','8v: Saturation in'])