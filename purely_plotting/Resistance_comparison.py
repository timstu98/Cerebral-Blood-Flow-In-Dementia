import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from pathlib import Path
import glob
import re
from tqdm import tqdm
import os 
import shutil 
import sys
from pylab import *
# import tikzplotlib
import matplotlib
matplotlib.use('pgf')

def import_and_set_up(folder_name,vessels_off,already_have_path=False):
    class data:
        size = (3.2*1.05,2.625*1.05)
        # size = (2*3.2,2*2.625)
    if not already_have_path:
        latex_folder = Path.cwd().parent.parent / 'LaTex files' / 'images' / folder_name
        if not os.path.exists(str(latex_folder)):
            latex_folder.mkdir() 
        else:
            # sys.exit('Folder ALREADY THERE')
            print('foler there')
        print('cont')

        data_folder = Path.cwd().parent / '_storage' / 'main' / folder_name
        s_f_len = len(str(data_folder)) + 1



        txt_s_path = str(data_folder) + "/*.txt"
        text_file_path = glob.glob(txt_s_path)[0]
        print('text file:', str(text_file_path)[s_f_len:])
        shutil.copyfile(str(text_file_path),str(latex_folder)+'/info.txt')

        shutil.copyfile(str(data_folder)+'/combined_network.csv',str(latex_folder)+'/combined_network.csv')
        print('copied')

        if not vessels_off:
            shutil.copyfile(str(data_folder)+'/combined_vessels.csv',str(latex_folder)+'/combined_vessels.csv')
    else:
        latex_folder = folder_name

    data.network = pd.read_csv(str(latex_folder)+'/combined_network.csv')
    if not vessels_off:
        data.vessels = pd.read_csv(str(latex_folder)+'/combined_vessels.csv')
    data.network = pd.read_csv(str(latex_folder)+'/combined_network.csv')
    data.network['Vn'] = 1-data.network['An']-data.network['Dn']
    data.network['Vp'] = 1-data.network['Ap']-data.network['Dp']
    data.t = data.network['t']
    data.p_diff = data.network['pressure_difference']
    data.network['norm R_tot'] = data.network['R_tot']/data.network.at[0,'R_tot'] 

    data.colours = {0:'#e41a1c',1:'#377eb8',2:'#4daf4a', 3:'#984ea3', 4:'#ff7f00'}

    return data

def settings(data,ax):
    plt.xlabel('Time [s]')
    # plt.xlim(0, t.max()*1.1)
    # pl    t.ylim(0, network[column_string].max()*1.1)

    ax.set_xlim(0)
    
    max_t = data.network.at[len(data.network)-1,'t']

    if max_t < 151 and max_t > 149:
        plt.xticks([0,30,60,90,120,150])
    elif max_t < 301 and max_t > 299:
        plt.xticks([0,60,120,180,240,300])
    elif max_t < 601 and max_t > 599:
        plt.xticks([0,150,300,450,600])
    elif max_t < 1201 and max_t > 1199:
        plt.xticks([0,200,400,600,800,1000,1200])
    elif max_t < 2401 and max_t > 2399:
        plt.xticks([0,400,800,1200,1600,2000,2400])

    # plt.grid(axis='y',which='minor',linewidth=0.5)
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.tick_params(left=True,right=True,bottom=True,top=True)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)

def crazy_plot():
    fig,ax = plt.subplots(figsize=size)
    baselineC_R = 0.00226883720930232
    R_C = baselineC_R/(network['phi']**4)
    R_C_b = R_C/baselineC_R
    baseline_R_tot = network.at[0,'R_tot']
    R_withoutC = baseline_R_tot-baselineC_R
    plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0],label=None, color='0.8')
    plt.plot(t,R_C_b,linewidth=1,linestyle="-",color='black')

    # plt.ylabel(column_label)

    settings(ax)
    ax.set_xlim(0)
    ax.set_ylim(0)

    # plt.axvline(x=0,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')
    plt.axhline(y=R_withoutC/baseline_R_tot,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')

    # save_name = column_string
    # plt.show()
    # fig.tight_layout()
    # plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')


    # baselineC_R = 0.00226883720930232
    # R_C = baselineC_R/(network['phi']**4)