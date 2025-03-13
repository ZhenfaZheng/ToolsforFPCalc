#/usr/bin/env python3

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt


figname = 'Absorp_BSE_phsc_1.png'
data = np.loadtxt('absorption_eh.dat')

figname = 'Absorp_BSE_phsc_1_noeh.png'
data = np.loadtxt('absorption_noeh.dat')

figsize_x = 4.8
figsize_y = 3.2 # in inches
fig, ax = plt.subplots()
fig.set_size_inches(figsize_x, figsize_y)

x = data[:,0]
y = data[:,1]

xmax = 5
index = np.where(x<xmax)[0]

ax.plot(x[index], y[index], 'r', lw=1.0)

ax.set_xlim(0, xmax)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Amplitude')

plt.tight_layout()
plt.savefig(figname, dpi=400)

print("\n%s has been saved.\n"%figname)

