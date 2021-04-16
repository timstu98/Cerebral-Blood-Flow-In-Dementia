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

folder_name = '001'

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

t = network['t']

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

    # tikzplotlib.clean_figure()
    # tikzplotlib.save(str(latex_folder)+ '/' +save_name + '.tikz',
    #         axis_height = '\\figureheight',
    #         axis_width = '\\figurewidth')

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

    # tikzplotlib.clean_figure()
    # tikzplotlib.save(str(latex_folder)+ '/' +save_name + '.tikz',
    #         axis_height = '\\figureheight',
    #         axis_width = '\\figurewidth')

def plt_double(column_string_1,column_string_2,column_label):
    fig,ax = plt.subplots(figsize=size)
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
    # 

size = (3.2,2.625)         

plt_single('Q_norm', 'Normalised flow')
plt_single('phi','Normalised Diamater')
print('singles')
plt_double('Ap','Dp','Proportion of pericytes\nin state')
plt_double('An','Dn','Proportion of neurons\nin state')
print('doubles')
plt_base('c','Concentration of AB',0.909090909)
plt_base('pt_volume_averaged','Volume averaged partial\ntissue oxygen pressure',44.113829801035656)
print('done')


