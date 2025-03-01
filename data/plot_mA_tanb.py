import sys
import csv
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axis import Axis

plt.rcParams["figure.figsize"] = (4,3)

data = genfromtxt('scan_over_mAtanb.csv', delimiter=',')

mask = (data[:, 0] >= 0.05)
data = data[mask, :]
#mask = (data[:, 1] <= 1.)
#data = data[mask, :]
mask = (data[:, 2] > 172.57 - 2*0.29) & (data[:, 2] < 172.57 + 2*0.29)
data = data[mask, :]

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots()

m = ax.scatter(data[:,5], data[:, 4], color='limegreen', s=1)
ax.set_xlabel(r'$m_A$ [GeV]')
ax.set_ylabel(r'$\tan \beta$')
ax.plot(1000, 1.85,'*',color='red',markersize=10)
Axis.set_rasterized(m, True)

plt.savefig('img/mA_vs_tanb.pdf', bbox_inches='tight', format='pdf', dpi=500, pad_inches=0)

