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

def settings(ax):
    plt.xlabel('Time [s]')
    # plt.xlim(0, t.max()*1.1)
    # plt.ylim(0, network[column_string].max()*1.1)
    # plt.grid(which='both')
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.tick_params(left=True,right=True,bottom=True,top=True)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)

all_pts = [65.67883277,62.95718794
,52.97305106
,40.33197468
,24.69455884
,20.38653852
,30.30177666
,30.90815354
,29.68516966
,29.8469117
,29.45720595
,28.89713436
,27.86541401]
all_names = ['A1',
'A2',
'A3',
'A4',
'A5',
'A6',
'C',
'V6',
'V5',
'V4',
'V3',
'V2',
'V1']   

x = np.arange(0,90)
M_constant = 8.2e-4
hill = M_constant*(x/(x+2))
size = (4,3)  
fig,ax = plt.subplots(figsize=size)
settings(ax)
plt.plot(x,hill/M_constant,linewidth=1.2,linestyle="-",color='black',label='Hill style model')
plt.hlines(1,0,x[-1],linewidth=1.2,linestyle="-",label='Constant metabolism', color='0.5')
# for i in range(len(all_pts)):
#     plt.vlines(all_pts[i],0,M_constant,label=all_names[i])

colours = {0:'#e41a1c',1:'#377eb8',2:'#4daf4a', 3:'#984ea3', 4:'#ff7f00'}

plt.vlines(all_pts[0],0,1,linewidth=1,linestyle="--",color=colours[0],label='A1')
plt.vlines(all_pts[4],0,1,linewidth=1,linestyle="--",color=colours[1],label='A5')
plt.vlines(all_pts[6],0,1,linewidth=1,linestyle="--",color=colours[2],label='C')
plt.vlines(all_pts[8],0,1,linewidth=1,linestyle="--",color=colours[3],label='V5')
plt.vlines(all_pts[12],0,1,linewidth=1,linestyle="--",color=colours[4],label='V1')

plt.legend()
plt.xlim(0)
plt.ylim(0)
ax.set_ylabel('Normalised oxygen metabolism')
ax.set_xlabel('$p_[t]$ (mmHg)')
# plt.show()
fig.tight_layout()
plt.savefig('metabolism.pgf')

#  plt.plot(t[0:len(A1_data)],A1_data/A1_data.iloc[0],linewidth=1,linestyle="-",color='red',label='A1')
#         # plt.plot(t[0:len(A2_data)],A2_data/A2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         # plt.plot(t[0:len(A3_data)],A3_data/A3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         # plt.plot(t[0:len(A4_data)],A4_data/A4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         plt.plot(t[0:len(A5_data)],A5_data/A5_data.iloc[0],linewidth=1,linestyle="-",color='green',label='A5')
#         # plt.plot(t[0:len(A6_data)],A6_data/A6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         plt.plot(t[0:len(C_data)],C_data/C_data.iloc[0],linewidth=1,linestyle="-",color='blue',label='C')
#         # plt.plot(t[0:len(V6_data)],V6_data/V6_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         plt.plot(t[0:len(V5_data)],V5_data/V5_data.iloc[0],linewidth=1,linestyle="-",color='magenta',label='V5')
#         # plt.plot(t[0:len(V4_data)],V4_data/V4_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         # plt.plot(t[0:len(V3_data)],V3_data/V3_data.iloc[0],linewidth=1,linestyle="-",color='blue')
#         # plt.plot(t[0:len(V2_data)],V2_data/V2_data.iloc[0],linewidth=1,linestyle="-",color='blue')
        # plt.plot(t[0:len(V1_data)],V1_data/V1_data.iloc[0],linewidth=1,linestyle="-",color='orange',label='V1')