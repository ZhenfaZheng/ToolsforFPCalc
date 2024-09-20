#!/usr/bin/env python3

import numpy as np

nkx = 9
nky = 9
nkz = 1

kcenter = [0.0, 0.0, 0.0]


kpoints = []
kstepx = 1.0 / nkx
kstepy = 1.0 / nky
kstepz = 1.0 / nkz
for ix in range(nkx):
    offsetx = kstepx * ix
    kx = kcenter[0] + offsetx
    for iy in range(nky):
        offsety = kstepy * iy
        ky = kcenter[1] + offsety
        for iz in range(nkz):
            offsetz = kstepz * iz 
            kz = kcenter[2] + offsetz
            kpoints.append([kx,ky,kz])

nks = len(kpoints)
weight = 1.0/nks

# print(nks)
for kpt in kpoints:
    print('%25.16f%25.16f%25.16f'%(kpt[0], kpt[1], kpt[2]))
