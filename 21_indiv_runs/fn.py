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

def plt_single(data,column_string,column_label):
    fig,ax = plt.subplots(figsize=data.size)
    ax.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')
    ax.plot(data.t,data.network[column_string],linewidth=1.2,linestyle="-",color='black',label='Model')

    plt.ylabel(column_label)

    settings(data,ax)
    # ax.set_xlim(0,max(data.network['t']))
    # plt.hlines(0.8033433301148208,0,150,linewidth=1,linestyle="--",color=data.colours[0],label='Steady state from\niterative solver')
    # ax.legend()

    ax.set_ylim(0)


    save_name = column_string
    # plt.show()
    fig.tight_layout()
    if not data.path:
        plt.savefig(str(latex_path) + '/' +save_name + '.pgf')
        # tizplotlib.clean_figure()
        # tikzplotlib.save(str(latex_path) + '/' +save_name + '.tex')
    else:
        plt.savefig(str(data.path) + '/' +save_name + '.pgf')
        # tikzplotlib.clean_figure()
        # tikzplotlib.save(str(data.path) + '/' +save_name + '.tex')

def plt_norm(data,column_string,column_label):
    fig,ax = plt.subplots(figsize=data.size)
    plt.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')
    plt.plot(data.t,data.network[column_string]/data.network.at[0,column_string],linewidth=1.2,linestyle="-",color='black')

    plt.ylabel(column_label)

    settings(data,ax)
    ax.set_ylim(0)

    save_name = column_string + '_norm'
    # plt.show()
    fig.tight_layout()
    if not data.path:
        plt.savefig(str(latex_path) + '/' +save_name + '.pgf')
    else:
        plt.savefig(str(data.path) + '/' +save_name + '.pgf')


def plt_base(data,column_string,column_label,baseline):
    fig,ax = plt.subplots(figsize=data.size)
    plt.plot(data.t,data.network[column_string],linewidth=1.2,linestyle="-",color='black')

    plt.ylabel(column_label)
    settings(data,ax)
    ax.set_ylim(0)

    plt.axhline(y=baseline,linewidth=0.5,linestyle='--',dashes=(5, 10), color='grey')

    save_name = column_string
    # plt.show()
    fig.tight_layout()
    if not data.path:
        plt.savefig(str(latex_path) + '/' +save_name + '.pgf')
    else:
        plt.savefig(str(data.path) + '/' +save_name + '.pgf')

def plt_triple(data,column_string_1,column_string_2,column_string_3,n_or_p,column_label):
    fig,ax = plt.subplots(figsize=data.size)
    plt.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')

    if n_or_p == 'p':
        plt.plot(data.t,data.network[column_string_1],linewidth=1.2,linestyle="-",color=data.colours[0],label='$A_[p]$')
        plt.plot(data.t,data.network[column_string_2],linewidth=1.2,linestyle="-",color=data.colours[1],label='$V_[p]$')
        plt.plot(data.t,data.network[column_string_3],linewidth=1.2,linestyle="-",color=data.colours[2],label='$D_[p]$')
    else:
        plt.plot(data.t,data.network[column_string_1],linewidth=1.2,linestyle="-",color=data.colours[0],label='$A_[n]$')
        plt.plot(data.t,data.network[column_string_2],linewidth=1.2,linestyle="-",color=data.colours[1],label='$V_[n]$')
        plt.plot(data.t,data.network[column_string_3],linewidth=1.2,linestyle="-",color=data.colours[2],label='$D_[n]$')

    plt.ylabel(column_label)
    
    settings(data,ax)
    # ax.set_ylim(0)

    ax.legend()

    save_name = column_string_1
    # plt.show()
    fig.tight_layout()
    if not data.path:
        plt.savefig(str(latex_path) + '/' +save_name + '.pgf')
        # tizplotlib.clean_figure()
        # tikzplotlib.save(str(latex_path) + '/' +save_name + '.tex')
    else:
        plt.savefig(str(data.path) + '/' +save_name + '.pgf')
        # tikzplotlib.clean_figure()
        # tikzplotlib.save(str(data.path) + '/' +save_name + '.tex')


def plt_vessels(data,column_string,column_label,normal=False):

    column_string = 'partial pressure blood(mmHg)'
    column_label = '$p_[t,norm]$'

    if not normal:
        column_label = '$p_[t]$'
    else:
        column_label = '$p_[t,norm]$'

    vessels = data.vessels



    A1_data = vessels[vessels['Name'] == 'A1'][column_string]
    A2_data = vessels[vessels['Name'] == 'A2'][column_string]
    A3_data = vessels[vessels['Name'] == 'A3'][column_string]
    A4_data = vessels[vessels['Name'] == 'A4'][column_string]
    A5_data = vessels[vessels['Name'] == 'A5'][column_string]
    A6_data = vessels[vessels['Name'] == 'A6'][column_string]
    C_data = vessels[vessels['Name'] == 'C'][column_string]
    V6_data = vessels[vessels['Name'] == 'V1'][column_string]
    V5_data = vessels[vessels['Name'] == 'V2'][column_string]
    V4_data = vessels[vessels['Name'] == 'V3'][column_string]
    V3_data = vessels[vessels['Name'] == 'V4'][column_string]
    V2_data = vessels[vessels['Name'] == 'V5'][column_string]
    V1_data = vessels[vessels['Name'] == 'V6'][column_string]



    fig,ax = plt.subplots(figsize=data.size)

    if not normal:
        plt.plot(t[0:len(A1_data)],A1_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A2_data)],A2_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A3_data)],A3_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A4_data)],A4_data,linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(A5_data)],A5_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A6_data)],A6_data,linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(C_data)],C_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V6_data)],V6_data,linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(V5_data)],V5_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V4_data)],V4_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V3_data)],V3_data,linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V2_data)],V2_data,linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(V1_data)],V1_data,linewidth=1,linestyle="-",color='blue')
    else:
        plt.plot(data.t[0:len(data.network)],data.p_diff[0:len(data.network)]/data.p_diff.iloc[0],linewidth=1,linestyle="--",label=None, color='0.8')


        plt.plot(data.t[0:len(A1_data)],A1_data/A1_data.iloc[0],linewidth=1.2,linestyle="-",color=data.colours[0],label='A1')
        # plt.plot(t[0:len(A2_data)],A2_data/A2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A3_data)],A3_data/A3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A4_data)],A4_data/A4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(data.t[0:len(A5_data)],A5_data/A5_data.iloc[0],linewidth=1.2,linestyle="-",color=data.colours[1],label='A5')
        # plt.plot(t[0:len(A6_data)],A6_data/A6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(data.t[0:len(C_data)],C_data/C_data.iloc[0],linewidth=1.2,linestyle="-",color=data.colours[2],label='C')
        # plt.plot(t[0:len(V6_data)],V6_data/V6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(data.t[0:len(V5_data)],V5_data/V5_data.iloc[0],linewidth=1.2,linestyle="-",color=data.colours[3],label='V5')
        # plt.plot(t[0:len(V4_data)],V4_data/V4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V3_data)],V3_data/V3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V2_data)],V2_data/V2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(data.t[0:len(V1_data)],V1_data/V1_data.iloc[0],linewidth=1.2,linestyle="-",color=data.colours[4],label='V1')

    plt.ylabel(column_label)
        
    settings(data,ax)
    plt.ylim(bottom=0)
    ax.legend()

    if not normal:
        save_name = 'pt_vessels'
    else:
        save_name = 'pt_vessels' + '_normal'
    # plt.show()
    fig.tight_layout()
    if not data.path:
        plt.savefig(str(latex_path) + '/' +save_name + '.pgf')
    else:
        plt.savefig(str(data.path) + '/' +save_name + '.pgf')