#!/usr/bin/env python

import numpy as np

#-----------------update at 2019/04/05 by zzf-----------------#

iband = 326
max_init_time = 500
nsample = 100

T = np.arange(1, max_init_time)
np.random.shuffle(T)
ini = np.empty([nsample,2])
ini[:,0] = T[:nsample]
ini[:,1] = iband

np.savetxt('INICON', ini, fmt='%5d')
