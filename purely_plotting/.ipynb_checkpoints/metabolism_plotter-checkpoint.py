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
# matplotlib.use('pgf')

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