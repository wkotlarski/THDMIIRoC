import sys
import csv
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (4,3)

data = genfromtxt('scan_over_mAtanb.csv', delimiter=',')

mask = (data[:, 0] >= 0.05)
data = data[mask, :]
#mask = (data[:, 1] <= 1.)
#data = data[mask, :]
mask = (data[:, 2] > 172.57 - 2*0.29) & (data[:, 2] < 172.57 + 2*0.29)
data = data[mask, :]

plt.rcParams['text.usetex'] = True
plt.scatter(data[:,5], data[:, 4], color='limegreen', s=1)
plt.xlabel(r'$m_A$ [GeV]')
plt.ylabel(r'$\tan \beta$')
plt.plot(1000, 1.85,'*',color='red',markersize=10)

plt.savefig('img/mA_vs_tanb.png', bbox_inches='tight', format='png', dpi=2000, pad_inches=0)

