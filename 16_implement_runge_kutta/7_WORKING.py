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
from pydub import AudioSegment
from pydub.playback import play
save = True
save_changes = True
save_vessels = True

print('\n' * 1)
print('############ NEW RUN ############')
print('\n' * 1)

### import required files
baselines = pd.read_csv(folder / '_storage/baselines' / 'baselines_met.csv', index_col=0)

pickle_in = open(folder / '_storage/baselines' / 'baselines_pt_averaged.pickle' , 'rb')
baseline_pt_averaged = pickle.load(pickle_in)
pickle_in.close()

print('baselines read in.')

### set parameters
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
delay_pressure_drop = 2
time_for_drop = 2
optimised_value_tau = 2.292929292929293
phi_min_baseline = 0.153
n=1
ratio_drop = 0.5
hill_constant = 2 # mmHg
kp_constant = 1/300 # s-1
kn_constant = 1/600 # s-1

#loop info
# max_time = 10
# no = 800000
max_time = 600
no = 400001

### set the arrays
t = pd.Series(np.linspace(0,max_time,no))
dt = max_time/(no-1)
drop_value = 60 - 34.18*ratio_drop

print('Values initialised.')

# PRessure stuff
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
    network_row['Q_norm'] = network_row['Q_tot'] / baselines_total_flow
    network_row['phi_min'] = 1-(alpha**0.25)*(1-phi_min_baseline)*(1-network_row['Q_norm'])**n


    vessels['pt_by_Vt'] = vessels['tissue partials(mmHg)'] * vessels['Vt(µm3)']
    network_row['pt_volume_averaged'] = vessels['pt_by_Vt'].sum()/13 
    
    #New bits:
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
        ks[4] =  (1/optimised_value_tau)*( -                    1                       + network_row['Ap']*(network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
        # ks[4] =  (1/optimised_value_tau)*( -                    1                       +                      (network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
    else:
        ks[4] =  (1/optimised_value_tau)*( - interpolate(network['phi'],i-delay_number) + network_row['Ap']*(network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
        # ks[4] =  (1/optimised_value_tau)*( - interpolate(network['phi'],i-delay_number) +                     (network_row['Q_norm']*(1-network_row['phi_min']) + network_row['phi_min']) )
    
    # # Euler from previous model. Necessary for: ...
    # if i == 0:                                  
    #     vessels = first_iteration.copy()
    # else:
    #     z = 0 # wanted to do something atleast
        
    #Update other values for new flow
    vessels.at[6,'Diameter(µm)'] = constant.at[6,'Diameter(µm)']*network_row['phi']    
    vessels.at[6,'X-Area wall(???)'] = np.pi * ( 0.16*(vessels.at[6,'Diameter(µm)'])**2 + 1.4*(vessels.at[6,'Diameter(µm)']) + 14 )
    vessels.at[6,'wall thickness(µm)'] = ( -(vessels.at[6,'Diameter(µm)']) + np.sqrt((vessels.at[6,'Diameter(µm)'])**2 + (4*(vessels.at[6,'X-Area wall(???)'])/np.pi))) / 2
    vessels['Sav'] = (vessels['Saturation in']+vessels['Saturation out'])/2
    vessels['M'] = (network_row['An'])*M_constant*vessels['tissue partials(mmHg)']/(vessels['tissue partials(mmHg)']+hill_constant)
    # print(vessels)
    for j in range(no_vessels):
        root_temp = None
        sols = 0
        root_temp = np.roots([1,0,150,23400*vessels.at[j,'Sav']/(vessels.at[j,'Sav']-1)]) #gives mmHg
        sols = check_imag_roots_real(root_temp) 
        vessels.at[j,'partial pressure blood(mmHg)'] = sols    

    vessels['Q in single(µm3/s)'] = first_iteration['Q in single(µm3/s)']*network_row['Q_norm']    # have just finished this.
    
    for j in range(len(vessels)):
        """For each vessel
        Update Q in single based on updated Q_tot
        Calculate dSoutdt and make this dSindt for next vessel."""
        #Could definitely get rid of this by not updating. I ahve all the values I realistically need with one column. Just for completeness. Could just have dSoutdt for each vessel and know this is thenext ones dSinDt. and know for A1 is is 0.
        if j ==0:
            ks[8].at[j,'dSindt'] = 0
        else:
            ks[8].at[j,'dSindt'] = ks[7].at[j-1,'dSoutdt']

        ks[7].at[j,'dSoutdt'] = 2/(np.pi*((vessels.at[j,'Diameter(µm)']**2)/4)*vessels.at[j,'Length(µm)']) \
                                * ( - (2 * np.pi * K * (vessels.at[j,'Diameter(µm)']/2) * constant.at[j,'Length(µm)'])/(vessels.at[j,'wall thickness(µm)']*cHb*H) \
                                * (vessels.at[j,'partial pressure blood(mmHg)']-vessels.at[j,'tissue partials(mmHg)'])\
                                - vessels.at[j,'Q in single(µm3/s)']*(vessels.at[j,'Saturation out']-vessels.at[j,'Saturation in']) )\
                                - ks[8].at[j,'dSindt']
        
    ## The main show: calculating dpdt and other useful representations.
    # print((1/(alpha_t*constant['Vt(µm3)'])))
    # print(((2 * np.pi * K * (vessels['Diameter(µm)']/2) * constant['Length(µm)'])/(vessels['wall thickness(µm)'])))
    # print((vessels['partial pressure blood(mmHg)']-vessels['tissue partials(mmHg)']))
    # print(vessels['M'] * constant['Vt(µm3)'])
    
    ks[6]['dptdt'] = (1/(alpha_t*constant['Vt(µm3)']))\
                    *(((2 * np.pi * K * (vessels['Diameter(µm)']/2) * constant['Length(µm)'])/(vessels['wall thickness(µm)']))\
                    *(vessels['partial pressure blood(mmHg)']-vessels['tissue partials(mmHg)'])\
                    - vessels['M'] * constant['Vt(µm3)'])

    if i ==3:
        # print(network.at[i,'t'])
        # print(empty_network_row(i))
        # print(empty_network_row['t'])
        # print(network)
        # print(vessels) 
        ignore = 0

    return ks,network_row,vessels

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

    # network_row_new

    return vessels, network_row_new

def empty_column(string_name):
    df = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'] })
    df[string_name] = 0 #np.nan
    return df

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

def empty_network_row(m):
    network_empty_row = pd.DataFrame()
    # network_empty_row.at[0,'t'] = network.at[m,'t']
    # network_empty_row.at[0,'pressure_difference'] = network.at[m,'pressure_difference']
    network_empty_row['phi'] = np.NaN
    network_empty_row['phi_min'] = np.NaN
    network_empty_row['R_tot'] = np.NaN
    network_empty_row['Q_tot'] = np.NaN
    network_empty_row['Q_norm'] = np.NaN
    network_empty_row['pt_averaged'] = np.NaN

    network_empty_row['kp_p'] = np.NaN
    network_empty_row['kp_n'] = np.NaN
    network_empty_row['Ap'] = np.NaN
    network_empty_row['Dp'] = np.NaN

    network_empty_row['kn_p'] = np.NaN
    network_empty_row['kn_n'] = np.NaN
    network_empty_row['An'] = np.NaN
    network_empty_row['Dn'] = np.NaN

    return network_empty_row

def RK4(i,network_row_main_updated,vessels_main_updated):
    #Can now drop the changing values from the vessels and also the network. Should create a one column vessels dataFrame just to record changes for the values. Could combine this for both. Could save as vessels_changes. eg two columns for both dptdt and dSoutdt. 
    #Need to work out how i want to consider dSoutdt and dSindt

    # vessels_in = vessels.copy()  # could probably just load in first_iteration or a zeros version (may do this) as we recalculate everything. Not strictly true as we use some of the constants values. Could also use first_iteration as that helps with the flow stuff
    
    network_row_in = network_row_main_updated.copy()
    vessels_in = vessels_main_updated.copy()

    # no need to update
    k1s , network_row_all_updated, vessels_all_updated = single_RK_run( i , network_row_in , vessels_in)
    # if i == 1 or i ==2 or i ==3:
    #     print(network_row_TIMESTEP)

    network_row = network_row_all_updated.copy()

    vessels , network_row = update_v_n(k1s,0.5,vessels_in,network_row_in)
    k2s , ignore_n, ignore_v= single_RK_run(i+0.5,network_row,vessels)

    vessels,network_row = update_v_n(k2s,0.5,vessels_in,network_row_in)
    k3s,ignore_n,ignore_v = single_RK_run(i+0.5,network_row,vessels)

    vessels,network_row = update_v_n(k3s,1,vessels_in,network_row_in)
    k4s,ignore_n,ignore_v = single_RK_run(i+1,network_row,vessels)

    all_ks = {}
    all_ks[0] = k1s
    all_ks[1] = k2s
    all_ks[2] = k3s
    all_ks[3] = k4s

    change = combine_ks(k1s,k2s,k3s,k4s)

    return change , network_row_all_updated , vessels_all_updated,all_ks

def network_set_up_values(vessels):
    # print(vessels)
    network.at[0,'R_tot'] = total_R(baselines,network.at[0,'phi'],alpha)
    network.at[0,'Q_tot'] = network.at[0,'pressure_difference'] / network.at[0,'R_tot']
    network.at[0,'Q_norm'] = network.at[0,'Q_tot'] / baselines_total_flow

    vessels['pt_by_Vt'] = vessels['tissue partials(mmHg)'] * vessels['Vt(µm3)']
    network.at[0,'pt_volume_averaged'] = vessels['pt_by_Vt'].sum()/13

    network.at[0,'kp_p'] = (kp_constant/2)*(1-np.tanh( (network.at[0,'pt_averaged']-pt_averaged_50) / pt_averaged_5))
    network.at[0,'kp_n'] = (kp_constant/2)*(1+np.tanh( (network.at[0,'pt_averaged']-pt_averaged_50) / pt_averaged_5 ))

    network.at[0,'kn_p'] = (kn_constant/2)*(1+np.tanh( (network.at[0,'c']-c_50) / c_5))
    network.at[0,'kn_n'] = (kn_constant/2)*(1-np.tanh( (network.at[0,'c']-c_50) / c_5 ))

def master(alpha,delay,baselines):
    RK4_change_alltime = {}
    vessels_alltime = {}
    ks_alltime = {}

    vessels = first_iteration.copy()
    network_set_up_values(vessels)

    network_row_main_updated = network.iloc[0,:]
    vessels_main_updated = vessels

    for i in tqdm(range(len(t)-1)): #Could maybe itterrows this column now as simply saving?
        #Coul;d create a list that would be easier to read then pass that through and unpack the other side.

        network_row_main_updated['t'] = network.at[i,'t']
        network_row_main_updated['pressure_difference'] = network.at[i,'pressure_difference']

        delta_values , network_row_all_updated, vessels_all_updated,ks = RK4(i, network_row_main_updated, vessels_main_updated)
        network.iloc[i,:] = network_row_all_updated
        vessels_alltime[i] = vessels_all_updated
        RK4_change_alltime[i] = delta_values
        ks_alltime[i] = ks

        vessels_main_updated, network_row_main_updated = update_v_n(delta_values, 1, vessels_all_updated, network_row_all_updated)
    
    ### To let you know whether it managed to run all the way through.
    if i == no-1:
        print('finished')
    else:
        print('broken')
    
    return network,vessels_alltime,RK4_change_alltime,ks_alltime

##############################################################################################################################################################################################################################################################################################

# Initialise Data, set up DFs
baseline = baselines.copy().drop(['Viscocity(mmHg*s)','Resistance for U', 'U in single(µm/s)'],axis=1).rename(columns={'Saturation ave': 'Sav'})
first_iteration = baseline.copy().loc[:,['Name','Number','Diameter(µm)','Length(µm)','Saturation in','Saturation out','Sav','Resistance for Q','Q in single(µm3/s)','Vt(µm3)','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)']] #.assign(dptdt = np.zeros(len(baseline)))
constant = baseline.copy().drop(['Q in single(µm3/s)','Pressure Drop(mmHg)','Saturation in','Saturation out','Resistance for Q','dS','Sav','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)'], axis=1)
pt_averaged = pd.Series(np.zeros(len(t)))
terms_data_alltime = {}
# print(baseline)
# print(first_iteration)
# print(constant)
# empty_network_row = pd.DataFrame({columns=['']})

alpha = 0.4
delay = 2.99

delay_number = round(delay/dt)
no_vessels = 13
baselines_total_flow = baselines['Q in single(µm3/s)'][0]

network = pd.DataFrame()
network['t'] = t
network['pressure_difference'] = pressure_difference
network['phi'] = pd.Series(np.zeros(no))
network['phi_min'] = pd.Series(np.zeros(no))
network['R_tot'] = pd.Series(np.zeros(no))
network['Q_tot'] = pd.Series(np.zeros(no))
network['Q_norm'] = pd.Series(np.zeros(no))
network['pt_averaged'] = pd.Series(np.zeros(no))

network['kp_p'] = pd.Series(np.zeros(no))
network['kp_n'] = pd.Series(np.zeros(no))
network['Ap'] = pd.Series(np.zeros(no))
network['Dp'] = pd.Series(np.zeros(no))

network['kn_p'] = pd.Series(np.zeros(no))
network['kn_n'] = pd.Series(np.zeros(no))
network['An'] = pd.Series(np.zeros(no))
network['Dn'] = pd.Series(np.zeros(no))

### Initialise values
network.at[0,'phi'] = 1
network.at[0,'phi_min'] = 1
# network.at[0,'dphidt'] = 0    # not necassarily
network.at[0,'Ap'] = 1
network.at[0,'Dp'] = 0
network.at[0,'pt_averaged'] = baseline_pt_averaged # kind of a bit of a different way to the rest. So I can embed the new parts where they should be in for loop.

pt_averaged_50 = baseline_pt_averaged*0.5
pt_averaged_5 = baseline_pt_averaged*0.1  

R = 1
k = 0.1
# tau_c = 600 
tau_c = 1/600 #s-1


c_baseline = R/(1+k)
network.at[0,'c'] = c_baseline
c_50 = c_baseline*10
c_5 = c_baseline
network.at[0,'An'] = 1
network.at[0,'Dn'] = 0

ks_names = pd.array(['0: dApdt','1: dDpdt','2: dAndt','3: dDndt','4: dphidt','5: dcdt','6:v dptdt','7v: dSoudt','8v: dSindt'])
values_2_change_names = pd.array(['0: Ap','1: Dp','2: An','3: Dn','4: phi','5: c','6:v tissue partials(mmHg)','7v: Saturation out','8v: Saturation in'])

print('Other arrays initialised.')


N,V,RK4_change_alltime,ks_alltime = master(alpha,delay,baselines)
print('Run complete.')

##############################################################################################################################################################################################################################################################################################

# pickle_out = open('network_first_long_run.pickle' , 'wb' )
# pickle.dump(N, pickle_out)
# pickle_out.close()

# pickle_out2 = open('RK4.pickle' , 'wb' )
# pickle.dump(RK4_change_alltime,pickle_out2)
# pickle_out2.close()
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

print('Save:', save)

save_combined = True

if save == True:
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
    