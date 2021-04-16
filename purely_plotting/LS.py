import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
matplotlib.use('pgf')

graph = pd.read_csv('response.csv')
points = pd.read_csv('points.csv')

def settings(ax):
    plt.xlabel('Time [s]')
    ax.set_xlim(0)
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

size = (4,3)
# size = (6,4)
fig,ax = plt.subplots(figsize=size)


# plt.plot(t[0:len(network)],p_diff[0:len(network)]/p_diff.iloc[0],label='pressure difference', color='0.8')
plt.plot(graph['t1'],graph['y1'],linewidth=1.2,linestyle="-",color='black',label='LS best fit')
plt.scatter(points['t'],points['d'], color='r',marker='x',label='Experimental data')

plt.ylabel('$\phi$') #rotation=0
ax.legend()

settings(ax)
# ax.set_xlim(0)
ax.set_ylim(0)

# plt.axvline(x=0,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')
# plt.axhline(y=0,linewidth=2,linestyle='--',dashes=(5, 10), color='grey')

# save_name = column_string
# # plt.show()
fig.tight_layout()
plt.savefig('LS_fit.pgf')

# plt.show()