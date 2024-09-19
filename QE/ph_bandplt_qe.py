#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt

prefix = 'graphene'
figname = 'Phonon.png'
phbandfile = prefix + '.freq.gp'
phband = np.loadtxt(phbandfile)


fig, ax = plt.subplots()
fig.set_size_inches(3.6,4.8)
mpl.rcParams['axes.unicode_minus'] = False

ax.plot(phband[:,0], phband[:,1:], 'r', lw=0.5)
ax.plot([phband[0,0], phband[-1,0]], [0,0], 'k', ls='--', lw=0.5)
ax.set_xlim(phband[0,0], phband[-1,0])
ax.set_xticks([])

plt.tight_layout()
plt.savefig(figname, dpi=250)

