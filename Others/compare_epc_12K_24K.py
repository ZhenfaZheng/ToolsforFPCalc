#!/usr/bin/env python

import numpy as np
import h5py


def main():

    bmin = 21
    bmax = 22

    filname = 'epc.h5'
    with h5py.File(filname, 'r') as f:
        dset_name = '/epmat/gmnvkq_real'
        dset = f[dset_name]
        g_real = dset[bmin-1:bmax][:,bmin-1:bmax]
        dset_name = '/epmat/gmnvkq_imag'
        dset = f[dset_name]
        g_imag = dset[bmin-1:bmax][:,bmin-1:bmax]
    g_24K = g_real + 1.0j * g_imag
    kpts_24K= read_h5file(filname, '/header/kpts')
    qpts_24K= read_h5file(filname, '/header/qpts')

    filname = '../4-12K-1B-allQ/epc.h5'
    with h5py.File(filname, 'r') as f:
        dset_name = '/epmat/gmnvkq_real'
        dset = f[dset_name]
        g_real = dset[bmin-1:bmax][:,bmin-1:bmax]
        dset_name = '/epmat/gmnvkq_imag'
        dset = f[dset_name]
        g_imag = dset[bmin-1:bmax][:,bmin-1:bmax]
    g_12K = g_real + 1.0j * g_imag
    kpts_12K= read_h5file(filname, '/header/kpts')
    qpts_12K= read_h5file(filname, '/header/qpts')


    k_ind_12K_in_24K = np.array([kindex(kpts_12K[ik], kpts_24K, 24, 24, 1) for ik in range(144)])

    g_24K_part = g_24K[:, :, :, k_ind_12K_in_24K, :]

    for iq in range(144):
        for ik in range(144):
            a = np.sum(np.abs(g_24K_part[:,:,:,ik,iq])) / 36
            b = np.sum(np.abs(g_12K[:,:,:,ik,iq])) / 36
            print()
            print("iq, ik = ", iq+1, ik+1)
            print('Average |g| of 24 k-grid: ', a)
            print('Average |g| of 12 k-grid: ', b)
            # for im in range(9):
            #     a = np.sum(np.abs(g_24K_part[:, :, im, ik,iq]))
            #     b = np.sum(np.abs(g_12K[:, :, im, ik, iq]))
            #     if (np.abs(a-b)>0.01):
            #         print()
            #         print('XXXXX')
            #         print(iq, ik, im)
            #         print(a)
            #         print(b)
            #     if (b > 0.01 and np.abs(a/b-1)>0.001):
            #         print()
            #         print('XXXXX')
            #         print(iq, ik, im)
            #         print(a)
            #         print(b)
            #     elif (b>0.01):
            #         print()
            #         print(iq, ik, im)
            #         print(a)

    '''
    g_real = read_h5file(filname, '/epmat/gmnvkq_real')
    g_imag = read_h5file(filname, '/epmat/gmnvkq_imag')

    g_real = g_real[12:14][:,12:14]
    g_imag = g_imag[12:14][:,12:14]
    '''


def kindex(kpt, kpts, n1, n2, n3):
    '''
    find index of kpt in kpts.

    Parameters:
        kpt  (ndarray): Array of k-point, shape (3).
        kpts (ndarray): Array of k-points, shape (nk, 3).
        n1 (int): number of k-point grid points in 1st dimension.
        n2 (int): number of k-point grid points in 2ed dimension.
        n3 (int): number of k-point grid points in 3rd dimension.

    Returns:
        index (int): index of kpt in kpts.
    '''

    k1 = int( round(kpt[0] * n1) % n1 )
    k2 = int( round(kpt[1] * n2) % n2 )
    k3 = int( round(kpt[2] * n3) % n3 )

    index = -1

    nk = kpts.shape[0]
    for ik in range(nk):
        i1 = int( round(kpts[ik,0] * n1) % n1 )
        if (i1 != k1): continue
        i2 = int( round(kpts[ik,1] * n2) % n2 )
        if (i2 != k2): continue
        i3 = int( round(kpts[ik,2] * n3) % n3 )
        if (i3 != k3): continue
        index = ik

    return index





def read_h5file(filname, dset_name=""):

    dset = None
    with h5py.File(filname, 'r') as f:
        dset = f[dset_name][()]

    return dset


if __name__=='__main__':
    main()
