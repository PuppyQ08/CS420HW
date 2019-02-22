# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:25:06 2019

@author: Aidan
"""

import numpy as np
from matplotlib import pyplot as plt

td = np.array([1, 4, 8, 16])
avg_time = np.array([307.31, 77.45, 39.46, 27.19])


one = np.array([308.097, 304.901, 309.126, 305.047, 309.379])
four = np.array([76.9488, 78.2809, 76.8299, 78.2478, 76.9655])
eight = np.array([39.4189, 39.8951, 39.2331, 39.9419, 38.8113 ])
sixteen = np.array([24.2101, 25.7484, 24.177, 31.6769, 30.1301 ])
to_plot = [one, four, eight, sixteen]
fig=plt.figure(1,figsize=(9,6))
ax=fig.add_subplot(111)
plt.xlabel("Number of threads")
plt.ylabel("Runtime (seconds)")
plt.title("Shared Memory Performance")
bp=ax.boxplot(to_plot)
fig.savefig('boxplot.png',bbox_inches='tight')
plt.show()

plt.xlabel("Number of shared memory threads")
plt.ylabel("Runtime (seconds)")
plt.title("Shared Memory Performance")
plt.scatter(td, avg_time)
plt.show()

efficiency = [0.992, 0.973, 0.706]

plt.xlabel("Number of shared memory threads")
plt.ylabel("Efficiency")
plt.title("Parallel efficiency")
plt.scatter(td[1:], efficiency)
plt.show()


objects = ("Static", "Dynamic")
y_pos = np.arange(len(objects))
time = [32.28, 21.06]
 
plt.bar(y_pos, time, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Runtime (seconds)')
plt.title('Parallel Scheduling')
plt.show()
