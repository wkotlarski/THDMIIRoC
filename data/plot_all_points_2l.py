import sys
import csv
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (4,3)

data = genfromtxt('all_points_2l.csv', delimiter=',')

plt.rcParams['text.usetex'] = True
plt.scatter(data[:,0], data[:, 1], color='limegreen', s=8)
plt.xlabel('p-value (HiggsSignals)')
plt.ylabel(r'$r/r_{95\%~C.L.}$ (HiggsBounds)')
plt.yscale('log')

plt.savefig('img/hs_vs_hb_2l.pdf', bbox_inches='tight', format='pdf', dpi=300, pad_inches = 0)

