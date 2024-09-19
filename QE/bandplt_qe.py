#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt

# Parameters Setting
emin = -6
emax =  6
prefix = 'graphene'
Efermi = -4.4444
figname = 'Band.png'
kpath = 'gmkg'

def main():
    
    bandsfile = prefix + '.bands'
    kpts, eig = readEig(bandsfile)

    kaxis = kptsCoordinate(kpts)
    bandsoutfile = 'bands.out'
    hskpts = HSkpts(bandsoutfile)

    plotBands(kaxis, eig, hskpts, Efermi)


def HSkpts(filename):
    '''
    Get x coordinate of high-symmetry k points from bands.out file.
    '''
    with open(filename) as f:
        lines = f.readlines()

    hskpts = []
    for line in lines:
        if 'high-symmetry point' in line:
            hskpts.append(float(line.split()[-1]))
    return hskpts
    

def readEig(filename):
    '''
    Read eigen value and k points from bands data file.
    '''
    with open(filename) as f:
        lines = f.readlines()

    headline = lines[0].replace(',','')
    nks  = int(headline.split()[-2])
    nbnd = int(headline.split()[2])
    kpts = np.zeros((nks, 3), dtype=float)
    eig  = np.zeros((nks, nbnd), dtype=float)

    # n is number of energy lines for each k point.
    if nbnd%10==0:
        n = nbnd // 10
    else:
        n= nbnd // 10 + 1

    numline = 1
    for i in range(nks):
        lineK = lines[numline]; numline += 1
        kpts[i] = np.array(lineK.split())
        states = []
        for j in range(n):
            lineE = lines[numline]; numline += 1
            states.extend(lineE.split())
        eig[i] = np.array(states)

    return kpts, eig

def kptsCoordinate(kpts):
    '''
    Calculate k axis coordinate of k points.
    '''
    kaxis = [0.0]
    for i in range(1,len(kpts)):
        s = kpts[i] - kpts[i-1]
        length = (s[0]**2 + s[1]**2 + s[2]**2)**0.5
        kaxis.append( kaxis[-1] + length )

    return kaxis


def HSkptLabels(kpath):
    '''
    Instruct k axis tick labels.
    '''
    kticklabels = []
    for s in kpath:
        if s=='g' or s=='G':
            s = r'$\Gamma$'
        else:
            s = s.upper()
        kticklabels.append(s)
    return kticklabels

def plotBands(kaxis, eig, hskpts, Efermi):
    '''
    Plot bands and save as Band.png fig.
    '''
    fig, ax = plt.subplots()
    fig.set_size_inches(3.6,4.8)
    mpl.rcParams['axes.unicode_minus'] = False

    eig -= Efermi
    ax.plot(kaxis, eig, 'r', lw=0.5)
    ax.plot([kaxis[0],kaxis[-1]], [0,0], 'k', ls='--', lw=0.5)
    for kpt in hskpts:
        ax.plot([kpt,kpt], [emin,emax], 'k', ls='--', lw=0.5)
    ax.set_xlim(kaxis[0], kaxis[-1])
    ax.set_ylim(emin, emax)
    ax.set_xticks(hskpts)
    ax.set_xticklabels(HSkptLabels(kpath))

    plt.tight_layout()
    plt.savefig(figname, dpi=250)


if __name__=='__main__':
    main()
