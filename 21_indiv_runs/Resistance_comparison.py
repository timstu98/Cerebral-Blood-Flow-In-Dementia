import fn
from pathlib import Path
import os
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('pgf')

import pandas as pd
import numpy as np

def R_plots(data):


    C_v_baseline = 4.513708e-05
    C_s_baseline = C_v_baseline/2**6

    Tot = data.network['R_tot']
    Tot_baseline = data.network.at[0,'R_tot']

    Without_C = Tot_baseline - C_s_baseline

    C_s = (data.alpha*C_v_baseline/data.network['phi']**4 + (1-data.alpha)*(C_v_baseline)  )/2**6


    C_s_normal = C_s/Tot_baseline
    Tot_normal = Tot/Tot_baseline

    print(C_s_baseline/Tot_baseline)



    fig,ax = plt.subplots(figsize=data.size)
    plt.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')
    plt.plot(data.t,C_s_normal,linewidth=1.2,linestyle="-",color=data.colours[0],label='$R_[C\:seg,norm]$')
    plt.plot(data.t,Tot_normal,linewidth=1.2,linestyle="-",color=data.colours[1],label='$R_[Tot,norm]$')

    

    fn.settings(data,ax)
    ax.set_ylabel('$R_[norm]$')
    ax.set_yscale('log')
    ax.set_xlim(0)
    # ax.set_ylim(0)
    ax.legend()
    fig.tight_layout()
    plt.show()
    # fig.savefig(str(latex_folder) + '/' + 'R_comp' + '.pgf')


    C_s_ratio = C_s/Tot
    # C_s_ratio2baseline


    fig1,ax1a = plt.subplots(figsize=data.size)

    fn.settings(data,ax1a)
    ax1b = ax1a.twinx()
    ax1a.tick_params(axis='y', labelcolor=data.colours[0])
    ax1b.tick_params(axis='y', labelcolor=data.colours[1])

    ax1a.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')

    ax1a.plot(data.t,C_s_ratio,linewidth=1.2,linestyle="-",color=data.colours[0],label='$R_[C\:seg]/R_[Tot]$')
    ax1b.plot(data.t,data.network['phi'],linewidth=1.2,linestyle="-",color=data.colours[1],label='$\phi$')

    # ax1.plot(data.t,C_s_ratio,linewidth=1,linestyle="-",color='red',label='ratio')
    # ax1.plot(data.t,data.network['phi'],linewidth=1,linestyle="-",color='blue',label='phi')

    ax1a.set_ylabel('Ratio of $R_[C\:seg]/R_[Tot]$',color=data.colours[0])
    ax1b.set_ylabel('$\phi$',color=data.colours[1])
    ax1a.set_yscale('log')
    ax1b.set_yscale('log')
    ax1a.set_xlim(0)

    ax1a.set_ylim(0.065,1.2)
    ax1b.set_ylim(0.065,1.2)

    lines1a, labels1a = ax1a.get_legend_handles_labels()
    lines1b, labels1b = ax1b.get_legend_handles_labels()
    ax1b.legend(lines1a + lines1b, labels1a + labels1b,loc=4)
    
    # ax1a.legend()
    # ax1b.legend
    fig1.tight_layout()
    plt.show()
    # plt.savefig(str(latex_folder) + '/' + 'R_vs_phi' + '.pgf')

    # fig2,ax2 = plt.subplots(figsize=data.size)
    # ax2.plot(data.network['phi'],Tot_normal,linewidth=1.2,linestyle="-",color=data.colours[0],label='temp')
    # # fn.settings(data,ax2)
    # ax2.set_yscale('log')
    # ax2.set_xscale('log')
    # plt.show()

    




    # save_name = column_string
    # fig.tight_layout()
    # plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')



print('start')
subdirectory = '/Users/Debs/OneDrive - Nexus365/4YP/LaTex files/images/results/default'
latex_folder = '/Users/Debs/OneDrive - Nexus365/4YP/LaTex files/images/results/default'
vessel_path = str(subdirectory) + '/combined_vessels.csv'

if not os.path.exists(vessel_path):
    vessels_off = True
else:
    vessels_off = False

data = fn.import_and_set_up(subdirectory,vessels_off,already_have_path=True)

data.path = str(subdirectory)
data.size = (9,6)
data.alpha = 0.3


print('plotting')

R_plots(data)

