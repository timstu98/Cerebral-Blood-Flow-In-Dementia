import pandas as pd 
import numpy as np 
from pathlib import Path
import glob
import re
from tqdm import tqdm


def dec_or_frac(string):
    if not re.search('([0-9]*\/[0-9]*)',string,re.IGNORECASE):
        value = float(string)
    else:
        value = eval((re.search('([0-9]*\/[0-9]*)',string,re.IGNORECASE).group(1)))
    return value
            
def final_value(folder_path_name,columns):
    db_red = pd.read_csv(folder_path_name + '/combined_network.csv')[columns]
    final_values = np.zeros(len(columns))
    for j in range(len(columns)):
        final_values[j] = db_red[columns[j]].iloc[-2]
    return final_values

s_f = Path.cwd().parent / '_storage' / 'main' 
s_f_len = len(str(s_f)) + 1

search = str(s_f / "*/")

main = []

# print(search)
all_folders = glob.glob(search)
# print(all_folders)
for i in tqdm(range(len(all_folders))):
    folder_path = all_folders[i]
    folder_name = folder_path[s_f_len:]
    remove = len(folder_path) + 1
    # print(folder_name)

    txt_s_path = folder_path+ "/*.txt"
    text_file_path = glob.glob(txt_s_path)[0]
    file= open(text_file_path,"rt") 
    tf = file.read()
    file.close()  
    
    time = float(re.search('Time of simulation: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    alpha = float(re.search('Alpha: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    delay = float(re.search('delay: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    p_drop = float(re.search('Ratio pressure drop: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    if re.search('tau_c: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE) == None:
        tau_c = 'Not printed'
    else:
        tau_c = dec_or_frac(re.search('tau_c: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    kp_constant = dec_or_frac(re.search('kp_constant: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    if re.search('pt_averaged_50: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE) == None:
        pt_50 = 'Not printed'
    else:
        pt_50 = float(re.search('pt_averaged_50: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    if re.search('pt_averaged_5: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE) == None:
        pt_5 = 'Not printed'
    else:
        pt_5 = float(re.search('pt_averaged_5: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    kn_constant = dec_or_frac(re.search('kn_constant: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    if re.search('c_50: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE) == None:
        c_50 = 'Not printed'
    else:
        c_50 = float(re.search('c_50: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    if re.search('c_5: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE) == None:
        c_5 = 'Not printed'
    else:
        c_5 = float(re.search('c_5: (-?\\b\\d+(?:\\.\\d+)?(?:/\\d+(?:\\.\\d+)?)?\\b)', tf, re.IGNORECASE).group(1))
    its = int(re.search('Number of iterations: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    delay_p_drop = float(re.search('Delay for pressure drop: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    time_p_drop = float(re.search('Time over which pressure drops: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    R = float(re.search('R: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))
    k = float(re.search('k: ([\d.]+)\s+(\S+)', tf, re.IGNORECASE).group(1))

    # if tau_c < 1e-5:
    #     tau_c = 1e-40
    # else:
    #     tau_c = tau_c
    # if kp_constant <1e-5:
    #     kp_constant = 1e-40
    # else:
    #     kp_constant = kp_constant
    # if kn_constant <1e-5:
    #     kn_constant = 1e-40
    # else:
    #     kn_constant = kn_constant

    network_s_path = folder_path + '/combined_network.csv'
    if not glob.glob(network_s_path):
        network_file = False
    else:
        network_file = True
    vessels_s_path = folder_path + '/combined_vessels.csv'
    if not glob.glob(vessels_s_path):
        vessels_file = False
    else:
        vessels_file = True

    final_values = final_value(folder_path,['Q_norm','c'])
    c_end = final_values[1]
    Q_end = final_values[0]

    main.append(    
        {
                'folder name':folder_name,
                'com vessels': vessels_file,
                'com network': network_file,
                'time': time,
                'alpha': alpha,
                'delay': delay,
                'pressure_drop':p_drop,
                'tau_c_f':tau_c,
                'tau_c_s': '1/'+str(round(1/tau_c)),
                'kp_constant_f':kp_constant,
                'kp_constant_s': '1/'+str(round(1/kp_constant)),
                'pt_50':pt_50,
                'pt_5':pt_5,
                'kn_constant_f':kn_constant,
                'km_constant_s': '1/'+str(round(1/kn_constant)),
                'c_50':c_50,
                'c_5':c_5,
                'Number of iterations':its,
                'delay pressure drop':delay_p_drop,
                'time for pressure drop':time_p_drop,
                'R':R,
                'k':k,
                'final Q': Q_end,
                'final c': c_end
        }
        )
    
out = pd.DataFrame(main)
out.to_csv('db.csv',index=False)