import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

db = pd.read_csv('db.csv')
print(len(db))

x_axis = ['final Q']
y_axis = ['final c']

x_axis = ['final Q']
y_axis = ['final c']
# point = ['']

# plt.scatter(db[x_axis],db[y_axis])
plt.scatter(db[(db[y_axis]>0).values][x_axis],db[(db[y_axis]>0).values][y_axis])
plt.show()