import pandas as pd
import numpy as np

import pickle

# ks_names = pd.array(['0: dApdt','1: dDpdt','2: dAndt','3: dDndt','4: dphidt','5: dcdt','6:v dptdt','7v: dSoudt'])
# ks = pd.array(np.zeros(len(ks_names)))
# k1s = ks.copy()
# k2s = ks.copy()
# k3s = ks.copy()
# k4s = ks.copy()

# print(ks_names)
# print(ks)
# print(k1s)
# print(len(ks_names))
# print(len(ks))

# t = np.linspace(0,100)
# pressure_difference = np.ones(len(t))
# no = len(t)

# network = pd.DataFrame()
# network['t'] = t
# network['pressure_difference'] = pressure_difference
# network['phi'] = pd.Series(np.zeros(no))
# network['phi_min'] = pd.Series(np.zeros(no))
# network['dphidt'] = pd.Series(np.zeros(no))
# network['R_tot'] = pd.Series(np.zeros(no))
# network['Q_tot'] = pd.Series(np.zeros(no))
# network['Q_norm'] = pd.Series(np.zeros(no))
# network['pt_averaged'] = pd.Series(np.zeros(no))

# network['pt_averaged_50'] = pd.Series(np.zeros(no))
# network['pt_averaged_5'] = pd.Series(np.zeros(no))
# network['kp_p'] = pd.Series(np.zeros(no))
# network['kp_n'] = pd.Series(np.zeros(no))
# network['dApdt'] = pd.Series(np.zeros(no))
# network['Ap'] = pd.Series(np.zeros(no))
# network['dDpdt'] = pd.Series(np.zeros(no))
# network['Dp'] = pd.Series(np.zeros(no))

# network_row = network.iloc[4,:].copy()
# network_row['Q_tot'] = 999
# print(network_row)

# ks_names = pd.array(['0: dApdt','1: dDpdt','2: dAndt','3: dDndt','4: dphidt','5: dcdt','6:v dptdt','7v: dSoudt'])
# ks = {}
# k1s = ks.copy()
# k2s = ks.copy()
# k3s = ks.copy()
# k4s = ks.copy()

# ks[7] = pd.DataFrame({'column_1':[0,1,4],'column_2':[9,5,6]})

# print(ks)

# ks[7] = 0.5* pd.DataFrame({'column_1':[0,1,4],'column_2':[9,5,6]})

# print(ks)

# egg = {}

# egg[3] = 6
# print(egg)

# def empty_column(string_name):
#     df = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'] })
#     df[string_name] = np.NaN
#     return df

# # print(empty_column('pt'))

# ks = {}
# ks[6] = empty_column('dptdt')
# ks[7] = empty_column('dSoutdt')
# ks[8] = empty_column('dSindt')

# print(ks[7])

# k1s = {}
# k1s[6] = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'],
#                         'pt': np.zeros(13)
                            
                            
#                             })

# k1s[6] = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'] })

# print(k1s)

# pickle_in = open("k1s.pickle","rb") 
# k1s = pickle.load(pickle_in)
# pickle_in.close() 

# pickle_in = open("k2s.pickle","rb") 
# k2s = pickle.load(pickle_in)
# pickle_in.close() 

# pickle_in = open("k3s.pickle","rb") 
# k3s = pickle.load(pickle_in)
# pickle_in.close() 

# pickle_in = open("k4s.pickle","rb") 
# k4s = pickle.load(pickle_in)
# pickle_in.close() 

# print(k1s)
# k1s = {0,2,3}

# for key,val in k1s.items():
#     2*val

# for i in range(len(k1s)):
#     if 0<=i<=5:
#         print(k1s[i])
#     else:
#         print('new')
#         print(k1s[i])
    

# for n in range(len(k1s)):
#     if type(k1s[n]) == "class 'pandas.core.frame.DataFrame'":
#         print(n)


# print(type(k1s[6]))

# def empty_column(string_name):
#     df = pd.DataFrame({'Name': ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'C', 'V6', 'V5', 'V4', 'V3', 'V2', 'V1'] })
#     df[string_name] = np.NaN
#     return df

# def combine_ks(k1s,k2s,k3s,k4s):

#     out = {}
#     out[6] = empty_column('dptdt')
#     out[7] = empty_column('dSoutdt')
#     out[8] = empty_column('dSindt')

#     for n in range(len(k1s)):
#         if 0<=n<=5:
#             out[n] = ( k1s[n] + k2s[n]*2 + k3s[n]*2 + k4s[n])/6
#         else:
#             out[n].iloc[:,1] = ( k1s[n].iloc[:,1] + k2s[n].iloc[:,1]*2 + k3s[n].iloc[:,1]*2 + k4s[n].iloc[:,1])/6

#     return out

# print('k1s',k1s,'k2s:',k2s,'k3s:',k3s,'k4s',k4s)

# out = combine_ks(k1s,k2s,k3s,k4s)
# print(out)

network_empty_row = pd.DataFrame()
network_empty_row.at[0,'t'] = network.at[m,'t']
network_empty_row.at[0,['pressure_difference'] = network.at[m,'pressure_difference']
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



        