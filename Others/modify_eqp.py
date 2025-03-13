#/usr/bin/env python3

import numpy as np


def main():

    fileqp = 'eqp_scph.dat'
    kpts, emf, eqp = read_eqp(fileqp)

    eig_uc = np.loadtxt('EIG_uc_ph.txt')

    sorted_indices = np.argsort(np.argsort(eqp.flatten())).reshape(eqp.shape)
    nk, nb = eqp.shape
    for ik in range(nk):
        for ib in range(nb):
            eqp[ik, ib] = eig_uc[sorted_indices[ik, ib]]

    eqp_all_data = np.loadtxt('eqp_scph_all.dat').reshape(nk, -1, 4)
    nb_all = eqp_all_data.shape[1] - 1

    bmin = 1
    bmax = 8
    eqp_all_data[:, bmin:bmax+1, 3] = eqp
    eqp_all_new = eqp_all_data.reshape(-1, 4)
    
    np.savetxt('eqp_sc_new.dat', eqp_all_new, fmt='%15.9f')



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


if __name__=='__main__':
    main()
