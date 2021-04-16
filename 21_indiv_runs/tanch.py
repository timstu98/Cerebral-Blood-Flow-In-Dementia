import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use('pgf')

x = np.arange(0,3,0.01)
k_constant = 1

x_50 = 1.5
x_5 = 0.5

k_plus = (k_constant/2)*(1+np.tanh( (x-x_50) / x_5))
k_minus = (k_constant/2)*(1-np.tanh( (x-x_50) / x_5))

def settings(ax):
    # plt.xlabel('Time [s]')
    # plt.xlim(0, t.max()*1.1)
    # plt.ylim(0, network[column_string].max()*1.1)
    # plt.grid(which='both')
    # ax.tick_params(axis="y",direction="in")
    # ax.tick_params(axis="x",direction="in")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.tick_params(bottom=False, labelbottom=False)


    # plt.axis('off')

bottom_value = (k_constant/2)*(1-np.tanh( (x_5) / x_5))
top_value = (k_constant/2)*(1+np.tanh( (x_5) / x_5))

fig,ax = plt.subplots(figsize = (4*1.25,2.5) )
settings(ax)
plt.ylabel('$k_[const]$',rotation=0)
# plt.xlabel('$x$',labelpad=20)
plt.text(2.5,-0.13,'$x$')

plt.plot(x,k_plus,color='green',linewidth=1.5,label='$k_[p]^+$, $k_[n]^-$')
plt.plot(x,k_minus,color='red',linewidth=1.5,label='$k_[p]^-$, $k_[n]^+$')

plt.axvline(x_50,linewidth=1,linestyle='--',dashes=(5, 10), color='grey')
plt.text(x_50-0.02,-0.13,'$x_[50]$',rotation=0)

plt.axvline(x_50+x_5,linewidth=1,linestyle='--',dashes=(5, 10), color='grey')
width=0.001
plt.arrow(x_50+0.02,0,x_5-0.02,0,width=width,length_includes_head=True,ec='black',fc='black',head_width=25*width)
plt.text(x_50+0.38*x_5,0.038,'$x_[5]$')
# plt.text(x_50+x_5,-0.1,rotation=0)

plt.axvline(x_50-x_5,linewidth=1,linestyle='--',dashes=(5, 10), color='grey')
width=0.001
plt.arrow(x_50-0.02,0,-x_5+0.02,0,width=width,length_includes_head=True,ec='black',fc='black',head_width=25*width)
plt.text(x_50-0.55*x_5,0.038,'$x_[5]$')
# plt.text(x_50-x_5,-0.1,rotation=0)

#+'\n*k_[const]'
plt.axhline(top_value,linewidth=1,linestyle='--',dashes=(5, 10), color='grey')
plt.text(-0.5,top_value-0.02,str(round(top_value,3)),rotation=0)
plt.axhline(bottom_value,linewidth=1,linestyle='--',dashes=(5, 10), color='grey')
plt.text(-0.5,bottom_value-0.02,str(round(bottom_value,3)),rotation=0)

ax.legend()


fig.tight_layout()
plt.savefig('tanch.pgf')
# plt.show()