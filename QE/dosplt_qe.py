#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob

def main():
    
    files = glob('*.dos')
    if len(files)==1:
        filename = files[0]
    else:
        filename = raw_input('\nPlease input filename of DOS data file: ')

    Efermi = readEf(filename)
    DOS = np.loadtxt(filename, dtype=float)
    plotDOS(DOS, Efermi)
    print('\nDone!\n')

def readEf(filename):

    with open(filename) as f:
        headline = f.readline()
    Efermi = float(headline.split()[-2])
    return Efermi

def plotDOS(DOS, Efermi):
    
    DOS[:,0] -= Efermi

    fig, ax = plt.subplots()
    fig.set_size_inches(4.8,3.2)
    mpl.rcParams['axes.unicode_minus'] = False

    #ax.plot(DOS[:,0], DOS[:,1], 'r', lw=0.3)
    ax.fill_between(DOS[:,0], DOS[:,1], color='r', lw=0.3)
    ax.set_xlim(-6, 6)
    ax.set_ybound(lower=0)
    ax.set_xlabel('Energy(eV)')
    ax.set_ylabel('DOS')

    plt.tight_layout()
    plt.savefig('DOS.png', dpi=400)

if __name__ == '__main__':
    main()
