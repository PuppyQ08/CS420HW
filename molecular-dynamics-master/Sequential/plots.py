# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 12:51:25 2018

@author: Aidan
"""

import numpy as np
from matplotlib import pyplot

f1 = open("average_speed.txt")
fstr = f1.read()
f1.close()

N = 1000

lst = fstr.split("\n")
step = np.array(range(N))
S = np.zeros(N)
for k in range(N):
    S[k] = float(lst[k])
    
pyplot.plot(step, S)

