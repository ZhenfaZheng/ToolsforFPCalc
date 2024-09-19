#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

path = './NATXT'

natxt = np.loadtxt(path)
nac = np.average(np.abs(natxt), axis=0)
n = int(np.sqrt(len(nac)))
nac.resize(n,n)

plt.imshow(nac, cmap='bwr', origin='lower')
plt.colorbar()
plt.title('NA Coupling')
plt.tight_layout()
plt.savefig('COUPLE.png')
