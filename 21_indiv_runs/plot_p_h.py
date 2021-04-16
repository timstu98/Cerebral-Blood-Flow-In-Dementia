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
import tikzplotlib
import matplotlib
matplotlib.use('pgf')

folder_name = 'just_p_h'

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
# shutil.copyfile(str(text_file_path),str(latex_folder)+'/info.txt')
# shutil.copyfile(str(data_folder)+'/combined_vessels.csv',str(latex_folder)+'/combined_vessels.csv')
# shutil.copyfile(str(data_folder)+'/combined_network.csv',str(latex_folder)+'/combined_network.csv')
print('copied')

network = pd.read_csv(str(latex_folder)+'/combined_network.csv')
vessels = pd.read_csv(str(latex_folder)+'/combined_vessels.csv')
network['Vn'] = 1-network['An']-network['Dn']
network['Vp'] = 1-network['Ap']-network['Dp']
t = network['t']
p_diff = network['pressure_difference']

# 0 - 0.1*t.max()     0 - 0.1*network[column_string].max()


def settings(ax):
    plt.xlabel('Time [s]')
    # plt.xlim(0, t.max()*1.1)
    # plt.ylim(0, network[column_string].max()*1.1)
    plt.grid(which='both')
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)

def plt_single(column_string,column_label):
    fig,ax = plt.subplots(figsize=size)
    plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0],label='pressure difference', color='0.8')
    plt.plot(t,network[column_string],linewidth=1,linestyle="-",color='black')

    plt.ylabel(column_label)

    settings(ax)
    ax.set_xlim(0)
    ax.set_ylim(0)

    # plt.axvline(x=0,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')
    # plt.axhline(y=0,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')

    save_name = column_string
    # plt.show()
    fig.tight_layout()
    plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')

def plt_base(column_string,column_label,baseline):
    fig,ax = plt.subplots(figsize=size)
    plt.plot(t,network[column_string],linewidth=1,linestyle="-",color='black')

    plt.ylabel(column_label)
    settings(ax)

    plt.axhline(y=baseline,linewidth=0.5,linestyle='--',dashes=(5, 10), color='grey')
    # plt.plot(t,np.zeros(len(t))+baseline,linewidth=0.5,linestyle='--',dashes=(5, 10), color='grey')

    ax.set_xlim(0)

    save_name = column_string
    # plt.show()
    fig.tight_layout()
    plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')


def plt_double(column_string_1,column_string_2,column_label):
    fig,ax = plt.subplots(figsize=size)
    plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0],label='pressure difference', color='0.8')
    plt.plot(t,network[column_string_1],linewidth=1,linestyle="-",color='blue')
    plt.plot(t,network[column_string_2],linewidth=1,linestyle="-",color='red')
    plt.ylabel(column_label)
    
    settings(ax)
    ax.set_xlim(0)
    ax.set_ylim(0)

    save_name = column_string_1
    # plt.show()
    fig.tight_layout()
    plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')

    # tikzplotlib.clean_figure()
    # tikzplotlib.save(str(latex_folder)+ '/' +save_name + '.tikz',
    #         axis_height = '\\figureheight',
    #         axis_width = '\\figurewidth')   

def plt_triple(column_string_1,column_string_2,column_string_3,column_label):
    fig,ax = plt.subplots(figsize=size)
    plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0], color='0.8',label=None)

    plt.plot(t,network[column_string_1],linewidth=1,linestyle="-",color='green',label=column_string_1)
    plt.plot(t,network[column_string_2],linewidth=1,linestyle="-",color='blue',label=column_string_2)
    plt.plot(t,network[column_string_3],linewidth=1,linestyle="-",color='red',label=column_string_3)
    plt.ylabel(column_label)
    
    settings(ax)
    ax.set_xlim(0)
    ax.set_ylim(0)

    ax.legend()

    save_name = column_string_1
    # plt.show()
    fig.tight_layout()
    plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')


def plt_vessels(column_string,column_label,normal=False):

    column_string = 'partial pressure blood(mmHg)'
    column_label = 'pressures'



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



    fig,ax = plt.subplots(figsize=size)

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
        plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0],label='pressure difference', color='0.8')


        plt.plot(t[0:len(A1_data)],A1_data/A1_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A2_data)],A2_data/A2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A3_data)],A3_data/A3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A4_data)],A4_data/A4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(A5_data)],A5_data/A5_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(A6_data)],A6_data/A6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(C_data)],C_data/C_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V6_data)],V6_data/V6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(V5_data)],V5_data/V5_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V4_data)],V4_data/V4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V3_data)],V3_data/V3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V2_data)],V2_data/V2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        plt.plot(t[0:len(V1_data)],V1_data/V1_data.iloc[0],linewidth=1,linestyle="-",color='blue')

    plt.ylabel(column_label)
        
    settings(ax)
    ax.set_xlim(0)
    ax.set_ylim(0)

    if not normal:
        save_name = column_label
    else:
        save_name = column_label + '_normal'
    # plt.show()
    fig.tight_layout()
    plt.savefig(str(latex_folder) + '/' +save_name + '.pgf')


size = (3.2,2.625)     


# Index(['Unnamed: 0', 't', 'pressure_difference', 'phi', 'phi_min', 'R_tot',
#        'Q_tot', 'Q_norm', 'pt_volume_averaged', 'kp_p', 'kp_n', 'Ap', 'Dp',
#        'kn_p', 'kn_n', 'An', 'Dn', 'c', 'dphidt', 'dApdt', 'dDpdt', 'dAndt',
#        'dDndt', 'dcdt'],
#       dtype='object')
  
# plt_single('Q_norm', 'Normalised flow')
# plt_single('phi', 'Normalised diamater')
# plt_single('phi_min','Pericyte constriction\nscaling factor')
# plt_single('R_tot','Total resiatance of model [SI]')

# plt_base('pt_volume_averaged','Volume averaged partial\ntissue oxygen pressure',44.113829801035656)
# plt_base('c','AB concentration',0.90909090909)

# plt_double('Ap','Dp','Proportion of pericytes\nin state')
plt_triple('Ap','Vp','Dp','Proportion of pericytes\nin state')

# plt_vessels('partial pressure blood(mmHg)','partial pressure blood')
# plt_vessels('partial pressure blood(mmHg)','partial pressure blood',normal=1)

# for col in vessels.columns:
#     print(col)

print('done')










