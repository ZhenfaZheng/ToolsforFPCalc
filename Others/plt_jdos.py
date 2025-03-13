#/usr/bin/env python3

import math, os, h5py
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt


def main():

    fileqp = 'eqp.dat'
    kpts, emf, eqp = read_eqp(fileqp)

    nvb = 3
    vbands = eqp[:,:nvb].flatten()
    cbands = eqp[:,nvb:].flatten()

    delist = np.array([np.abs(cc-vv) for cc in cbands for vv in vbands])

    edos = gen_dos(delist, sigma=0.005)
    plt_dos(edos, enlist=delist, figname='JDOS_scph.png')

    print("\nDONE!\n")


def read_eqp(filname='eqp.dat', inp=None):

    data = np.loadtxt(filname)
    n_spin_bands = int(data[0,-1])
    nks = int(data.shape[0] / (n_spin_bands + 1))

    if (nks*(n_spin_bands+1)*4 != data.shape[0]*data.shape[1]):
        error_exit("The number of rows and columns in " \
                   "%s seems to be incorrect!"%filname)

    data = data.reshape(nks, n_spin_bands+1, 4)

    kpts = data[:, 0, :3]
    emf = data[:, 1:, 2] # mean-field energy
    eqp = data[:, 1:, 3] # quasi-particle (GW) energy

    if inp is not None:
        index = data[0, 1:, 1]
        bmin = int(inp['band_index_min'])
        bmax = int(inp['band_index_max'])
        bindex = np.where((index >= bmin) & (index <= bmax))[0]
        emf = emf[:, bindex]
        eqp = eqp[:, bindex]

    return kpts, emf, eqp


def gen_dos(enlist, sigma=0.05):

    enlist = np.array(enlist)

    emin = enlist.min() - 10 * sigma
    emax = enlist.max() + 10 * sigma
    ex = np.linspace(emin, emax, 1000)

    dos = np.zeros_like(ex)
    for en in enlist:
        dos += gaussian(ex, en, sigma)
    dos /= len(enlist)

    edos = np.vstack((ex, dos)).T

    return edos


def gaussian(x, mu, sigma):

    return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((x - mu) / sigma)**2)


def plt_dos(edos, enlist=None, figname='DOS.png'):

    figsize_x = 6.4
    figsize_y = 4.8 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    dos_scale = 1.0
    if (enlist is not None):
        counts, bin_edges = np.histogram(enlist, bins=500)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        ax.vlines(bin_centers, ymin=0, ymax=counts, color='black', linewidth=0.1)
        dos_scale = counts.max() / edos[:,1].max() * 5.0
    
    ax.plot(edos[:,0], edos[:,1] * dos_scale, 'r', lw=1.0)

    ax.set_xlim(2,5)
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('DOS (a.u.)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)

if __name__=='__main__':
    main()
