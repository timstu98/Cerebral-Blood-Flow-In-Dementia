import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pylab
# import seaborn
from pylab import *
import tikzplotlib

x = np.arange(0,10,step = 0.001)
y = np.sin(x)
# print(x)
# print(y)
# 
fig = plt.subplots()
plt.plot(x,y)
xlabel('$x$-axis')
ylabel('$y$-axis')
tikzplotlib.clean_figure()
tikzplotlib.save('fig.tikz',
           axis_height = '\\figureheight',
           axis_width = '\\figurewidth')
print('done')
