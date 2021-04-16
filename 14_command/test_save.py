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

pickle_in = open(folder / '_storage/baselines' / 'baselines_pt_averaged.pickle' , 'rb')
baseline_pt_averaged = pickle.load(pickle_in)
pickle_in.close()

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
delay_pressure_drop = 2
time_for_drop = 2
optimised_value_tau = 2.292929292929293
phi_min_baseline = 0.153
n=1
ratio_drop = 0.2
hill_constant = 2 # mmHg
kp_constant = 1/300 # s-1
kn_constant = 1/600 # s-1

#loop info
max_time = 600
no = 18000001

### set the arrays
t = pd.Series(np.linspace(0,max_time,no))
dt = max_time/(no-1)
drop_value = 60 - 34.18*ratio_drop

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

def master(alpha,delay,baselines):
    """Entire model.
    
    """  
    ###Set up
    
    delay_number = round(delay/dt) #finds nearest index to the 
    no_vessels = 13
    
    network = pd.DataFrame()
    network['t'] = t
    network['pressure_difference'] = pressure_difference
    network['phi'] = pd.Series(np.zeros(no))
    network['phi_min'] = pd.Series(np.zeros(no))
    network['dphidt'] = pd.Series(np.zeros(no))
    network['R_tot'] = pd.Series(np.zeros(no))
    network['Q_tot'] = pd.Series(np.zeros(no))
    network['Q_norm'] = pd.Series(np.zeros(no))
    network['pt_averaged'] = pd.Series(np.zeros(no))
    
    network['pt_averaged_50'] = pd.Series(np.zeros(no))
    network['pt_averaged_5'] = pd.Series(np.zeros(no))
    network['kp_p'] = pd.Series(np.zeros(no))
    network['kp_n'] = pd.Series(np.zeros(no))
    network['dApdt'] = pd.Series(np.zeros(no))
    network['Ap'] = pd.Series(np.zeros(no))
    network['dDpdt'] = pd.Series(np.zeros(no))
    network['Dp'] = pd.Series(np.zeros(no))
    
    ### Initialise values
    
    network.at[0,'phi'] = 1
    network.at[0,'phi_min'] = 1
    network.at[0,'dphidt'] = 0
    network.at[0,'Ap'] = 1
    network.at[0,'Dp'] = 0
    network.at[0,'pt_averaged'] = baseline_pt_averaged # kind of a bit of a different way to the rest. So I can embed the new parts where they should be in for loop.
    
    pt_averaged_50 = baseline_pt_averaged*0.5
    pt_averaged_5 = baseline_pt_averaged*0.1
    
    R = 1
    k = 0.1
    tau_c = 1/600 #s-1
    c_baseline = R/(1+k)
    network.at[0,'c'] = c_baseline
    c_50 = c_baseline*10
    c_5 = c_baseline
    network.at[0,'An'] = 1
    network.at[0,'Dn'] = 0
    
    
    for i in tqdm(range(len(t))):
        """For each time step:
        
        """  
        ### Update network paramaters
        
        network.at[i,'R_tot'] = total_R(baselines,network.at[i,'phi'],alpha)
        network.at[i,'Q_tot'] = network.at[i,'pressure_difference']/network.at[i,'R_tot']
        network.at[i,'Q_norm'] = network.at[i,'Q_tot'] / baselines['Q in single(µm3/s)'][0]
        
        #New bits:
        network.at[i,'kp_p'] = (kp_constant/2)*(1-np.tanh( (network.at[i,'pt_averaged']-pt_averaged_50) / pt_averaged_5))
        network.at[i,'kp_n'] = (kp_constant/2)*(1+np.tanh( (network.at[i,'pt_averaged']-pt_averaged_50) / pt_averaged_5 ))
        network.at[i,'dApdt'] = -network.at[i,'kp_p']*network.at[i,'Ap'] + network.at[i,'kp_n']*(1-network.at[i,'Ap']-network.at[i,'Dp'])
        network.at[i,'dDpdt'] = network.at[i,'kp_p']*(1-network.at[i,'Ap']-network.at[i,'Dp'])
        
        network.at[i,'kn_p'] = (kn_constant/2)*(1+np.tanh( (network.at[i,'c']-c_50) / c_5))
        network.at[i,'kn_n'] = (kn_constant/2)*(1-np.tanh( (network.at[i,'c']-c_50) / c_5 ))
        network.at[i,'dAndt'] = -network.at[i,'kn_p']*network.at[i,'An'] + network.at[i,'kn_n']*(1-network.at[i,'An']-network.at[i,'Dn'])
        network.at[i,'dDndt'] = network.at[i,'kn_p']*(1-network.at[i,'An']-network.at[i,'Dn'])
        
        network.at[i,'dcdt'] = (1/tau_c)*(R-network.at[i,'c']*(k+network.at[i,'Q_norm']))
 

        if i-delay_number <0: # As cant index before 0 as not initialised for that.
            network.at[i,'dphidt'] =  (1/optimised_value_tau)*( - 1 + network.at[i,'Q_norm']*(1-network.at[i,'phi_min']) + network.at[i,'phi_min'] )
        else:
            network.at[i,'dphidt'] =  (1/optimised_value_tau)*( -network.at[i-delay_number,'phi'] + network.at[i,'Ap']*(network.at[i,'Q_norm']*(1-network.at[i,'phi_min']) + network.at[i,'phi_min']) )
            
        #Euler    
        network.at[i+1,'phi'] = network.at[i,'dphidt']*dt + network.at[i,'phi']
        network.at[i+1,'phi_min'] = 1-(alpha**0.25)*(1-phi_min_baseline)*(1-network.at[i,'Q_norm'])**n  
        
        network.at[i+1,'Ap'] = network.at[i,'dApdt']*dt + network.at[i,'Ap']   
        network.at[i+1,'Dp'] = network.at[i,'dDpdt']*dt + network.at[i,'Dp']  
        
        network.at[i+1,'An'] = network.at[i,'dAndt']*dt + network.at[i,'An']   
        network.at[i+1,'Dn'] = network.at[i,'dDndt']*dt + network.at[i,'Dn']  
        network.at[i+1,'c'] = network.at[i,'dcdt']*dt + network.at[i,'c']  
        #### Update Vessel paramaters
        
        # Euler from previous model. Necessary for: ...
        if i == 0:                                  
            vessels = first_iteration.copy()
        else:
            # vessels = vessel_data_alltime[i-1].copy() ##### REMOVED THIS ON THE 3_REDUDEC_SAVE
            #Following works as I have loaded in previous values. Then updateding for the new ones Forwards Euler.
            vessels['tissue partials(mmHg)'] = vessels['tissue partials(mmHg)'] + vessels['dptdt']*dt
            vessels['Saturation out'] = vessels['Saturation out'] + vessels['dSoutdt']*dt
            vessels['Saturation in'] = vessels['Saturation in'] + vessels['dSindt']*dt 
#             vessels['Sav'] = vessels['Sav'] + vessels['dSavdt']*dt
        #Update other values for new flow
        vessels.at[6,'Diameter(µm)'] = baseline.at[6,'Diameter(µm)']*network.at[i,'phi']    
        vessels.at[6,'X-Area wall(???)'] = np.pi * ( 0.16*(vessels.at[6,'Diameter(µm)'])**2 + 1.4*(vessels.at[6,'Diameter(µm)']) + 14 )
        vessels.at[6,'wall thickness(µm)'] = ( -(vessels.at[6,'Diameter(µm)']) + np.sqrt((vessels.at[6,'Diameter(µm)'])**2 + (4*(vessels.at[6,'X-Area wall(???)'])/np.pi))) / 2
        vessels['Sav'] = (vessels['Saturation in']+vessels['Saturation out'])/2
        vessels['M'] = (network.at[i,'An'])*M_constant*vessels['tissue partials(mmHg)']/(vessels['tissue partials(mmHg)']+hill_constant)
        for j in range(no_vessels):
            root_temp = None
            sols = 0
            root_temp = np.roots([1,0,150,23400*vessels.at[j,'Sav']/(vessels.at[j,'Sav']-1)]) #gives mmHg
            sols = check_imag_roots_real(root_temp) 
            vessels.at[j,'partial pressure blood(mmHg)'] = sols    
        #Could just do the flow update as ratio of previous flow to new flow then times previous values by this ratio as all flows scale linearly.
        vessels['Q in single(µm3/s)'] = vessels['Q in single(µm3/s)']*network.at[i,'Q_norm']  # DONT QUITWE UNDERSTAND WHY THIS WORKS. CAN REMEMBER THAT IT WORKED THOUGH JUST NEED TO UNDERstAND
        
        for j in range(len(vessels)):
            """For each vessel
            Update Q in single based on updated Q_tot
            Calculate dSoutdt and make this dSindt for next vessel."""
            #Could definitely get rid of this by not updating. I ahve all the values I realistically need with one column. Just for completeness. Could just have dSoutdt for each vessel and know this is thenext ones dSinDt. and know for A1 is is 0.
            if j ==0:
                vessels.at[j,'dSindt'] = 0
            else:
                vessels.at[j,'dSindt'] = vessels.at[j-1,'dSoutdt']

            vessels.at[j,'dSoutdt'] = 2/(np.pi*((vessels.at[j,'Diameter(µm)']**2)/4)*vessels.at[j,'Length(µm)']) \
                                    * ( - (2 * np.pi * K * (vessels.at[j,'Diameter(µm)']/2) * constant.at[j,'Length(µm)'])/(vessels.at[j,'wall thickness(µm)']*cHb*H) \
                                    * (vessels.at[j,'partial pressure blood(mmHg)']-vessels.at[j,'tissue partials(mmHg)'])\
                                    - vessels.at[j,'Q in single(µm3/s)']*(vessels.at[j,'Saturation out']-vessels.at[j,'Saturation in']) )\
                                    - vessels.at[j,'dSindt']
            
        ### The main show: calculating dpdt and other useful representations.
        
        vessels['dptdt'] = (1/(alpha_t*constant['Vt(µm3)']))\
                        *(((2 * np.pi * K * (vessels['Diameter(µm)']/2) * constant['Length(µm)'])/(vessels['wall thickness(µm)']))\
                        *(vessels['partial pressure blood(mmHg)']-vessels['tissue partials(mmHg)'])\
                        - vessels['M'] * constant['Vt(µm3)'])
        vessels['pt_by_Vt'] = vessels['tissue partials(mmHg)'] * vessels['Vt(µm3)']

        network.at[i+1,'pt_averaged'] = vessels['pt_by_Vt'].sum()/13 ###### BADDDDDDd

        ### Save vessels to a dictionary

        # if i == 500:
        #     print(vessels)
        #     print(network)
        
        # vessel_data_alltime[i] = vessels
    
    ### To let you know whether it managed to run all the way through.
    
    if i == no-1:
        print('finished')
    else:
        print('broken')
    
    return network#, vessel_data_alltime

# Initialise Data, set up DFs
baseline = baselines.drop(['Viscocity(mmHg*s)','Resistance for U', 'U in single(µm/s)'],axis=1).rename(columns={'Saturation ave': 'Sav'})
first_iteration = baseline.loc[:,['Name','Number','Diameter(µm)','Length(µm)','Saturation in','Saturation out','Sav','Resistance for Q','Q in single(µm3/s)','Vt(µm3)','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)']].assign(dptdt = np.zeros(len(baseline)))
constant = baseline.drop(['Q in single(µm3/s)','Pressure Drop(mmHg)','Saturation in','Saturation out','Resistance for Q','dS','Sav','partial pressure blood(mmHg)','X-Area wall(???)','wall thickness(µm)','tissue partials(mmHg)'], axis=1)
# vessel_data_alltime = {}
pt_averaged = pd.Series(np.zeros(len(t)))
terms_data_alltime = {}

alpha = 0.2
delay = 2.99

# N,V_alltime = master(alpha,delay,baselines)
N = master(alpha,delay,baselines)
# plt.plot(N['t'],N['Q_norm'])
# plt.show()
# plt.plot(N['t'],N['c'])
# plt.show()
print('min:',min(N['Q_norm']))
print('final:', N.at[len(N)-2,'Q_norm'])

txt_file = [
    'Time of simulation: ' + str(max_time) + ' seconds\n' ,
    'Number of iterations: ' + str(no) + '\n',
    'therefore time step is: ' + str(dt) + ' seconds\n' ,
    'Delay for pressure drop: ' + str(delay_pressure_drop) + ' seconds\n' ,
    'Time over which pressure drops: ' + str(time_for_drop) + ' seconds\n',
    'Ratio pressure drop: ' + str(ratio_drop) + '\n',
    'kp_constant: ' + str(kp_constant) + '\n',
    'kn_constant: ' + str(kn_constant) + '\n',
    'Date and time: ' + str(datetime.datetime.now())
           ]
print(txt_file)
if save == True:
    count = 1
    while count < 1000:
    
        if len(str(count)) == 3:
            subfoulder = str(i)
        elif len(str(i)) == 2:
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

    print('hi')
    file= open(main_path / "info.txt","w") 
    file.writelines(txt_file)
    file.close()  
    print('text file saved')
    
    # pickle_out = open(main_path / 'vessels_first_long_run.pickle' , 'wb' )
    # pickle.dump(V_alltime, pickle_out)
    # pickle_out.close()
    # print('Vessels saved')
    
    pickle_out2 = open(main_path / 'network_first_long_run.pickle' , 'wb' )
    pickle.dump(N, pickle_out2)
    pickle_out2.close()
    print('Network saved')
    
    plt.plot(N['t'],N['Q_norm'])
    plt.ylabel('Normalised Flow')
    plt.xlabel('Time / seconds')
    plt.grid(which='both')
    plt.savefig(main_path / 'Q_vs_t')
    