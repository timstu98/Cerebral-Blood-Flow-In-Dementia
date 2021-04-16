import fn
from pathlib import Path
import os
import matplotlib.pyplot as plt
import sys

results_folder = Path.cwd().parent.parent / 'LaTex files' / 'images' / 'results'

subdirectories = [x for x in results_folder.iterdir() if x.is_dir()]
# print(subdirectories)
# sys.exit()


print('start')

latex_folder = 'start'

# subdirectory = 

# data.save_or_show = 'show'

for i in range(len(subdirectories)): #range(len(subdirectories)) range(1)
    # subdirectory = '/Users/Debs/OneDrive - Nexus365/4YP/LaTex files/images/results/inc_peri'
    subdirectory = subdirectories[i]
    print('------------------------------------------------------------------------------------')
    print(i+1,'/',len(subdirectories))
    print(str(subdirectory))
    latex_folder = subdirectory

    vessel_path = str(subdirectory) + '/combined_vessels.csv'

    if not os.path.exists(vessel_path):
        vessels_off = True
    else:
        vessels_off = False

    # print(vessels_off)

    data = fn.import_and_set_up(subdirectory,vessels_off,already_have_path=True)

    data.path = str(subdirectory)
    print(data.path)

    #TEXT COLOR IDEA
    fn.plt_single(data,'Q_norm', '$Q_[norm]$')
    fn.plt_single(data,'phi','$\phi$')
    fn.plt_single(data,'phi_min','$\phi_[min]$')
    # fn.plt_single(data,'norm R_tot','$R_[norm]$')
    # print('singles')
    fn.plt_triple(data,'Ap','Vp','Dp','p','Proportion of pericytes')
    fn.plt_triple(data,'An','Vn','Dn','n','Proportion of neurons')
    # # print('doubles')
    # fn.plt_base(data,'c','$c$',0.909090909)
    # fn.plt_base(data,'pt_volume_averaged','$p_[t]$',44.113829801035656)

    fn.plt_norm(data,'c','$c_[norm]$')
    fn.plt_norm(data,'pt_volume_averaged','$\hat{p}_[t,norm]$')
    fn.plt_norm(data,'R_tot','$R_[norm]$')

    if not os.path.exists(vessel_path):
        tim = 'hi'
    else:
        fn.plt_vessels(data,'partial pressure blood(mmHg)','partial pressure blood',normal=1)

    plt.close('all')

print('DONE')