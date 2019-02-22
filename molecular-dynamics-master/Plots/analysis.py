import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend('agg')

result = np.loadtxt(open("Book.csv", "rb"), delimiter=",", skiprows=0)

print(result[:,1])
plt.boxplot([result[:,1],result[:,2],result[:,3],result[:,4]])

plt.savefig('1.png')

print(result)

