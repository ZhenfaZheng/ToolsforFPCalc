#!/usr/bin/env python

import numpy as np
import h5py


def main():

    iq = 5
    im = 8

    ib = 22
    jb = 22
    ik = 1

    gnew = 0.3 + 0.0j

    modify_epc(iq, im, ib, jb, ik, gnew, filname='epc.h5')


def modify_epc(iq, im, ib, jb, ik, gnew, filname='epc.h5'):

    iq -= 1; im -= 1
    ib -= 1; jb -= 1; ik-= 1

    kpts = read_h5file(filname, '/header/kpts')
    qpts = read_h5file(filname, '/header/qpts')

    g_real = read_h5file(filname, '/epmat/gmnvkq_real')
    g_imag = read_h5file(filname, '/epmat/gmnvkq_imag')

    jq = index(-qpts[iq], qpts)
    ikq = index(kpts[ik]+qpts[iq], kpts)

    # gmnv(k,q) = < nk+q | dV_qv | mk >
    g_real[ib, jb, im, ik, iq] = gnew.real
    g_imag[ib, jb, im, ik, iq] = gnew.imag

    # gnmv(k+q,-q) = < mk | dV_-qv | nk+q >
    g_real[jb, ib, im, ikq, jq] = gnew.real
    g_imag[jb, ib, im, ikq, jq] = - gnew.imag

    with h5py.File(filname, 'r+') as ff:
        ff['/epmat/gmnvkq_real'][()] = g_real
        ff['/epmat/gmnvkq_imag'][()] = g_imag


def index(kpt, kpts, norm=1.0e-4):

    nks = kpts.shape[0]
    for ik in range(nks):
        for iax in range(3):
            dk = kpts[ik, iax] - kpt[iax]
            diff = abs(dk - round(dk))
            if (diff > norm):
                break
            if (iax == 2): return ik

    kstr = "(%.4f, %.4f, %.4f)"%(kpt[0], kpt[1], kpt[2])
    error_exit("k/q-point %s not found!"%kstr)
    # return -1


def read_h5file(filname, dset_name=""):

    dset = None
    with h5py.File(filname, 'r') as f:
        dset = f[dset_name][()]

    return dset



if __name__=='__main__':
    main()
