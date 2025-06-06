#/usr/bin/env python

import math, os, h5py
import numpy as np


def main():

    inp = read_inp()
    load_inf(inp, 'bse_1/bsemat.h5')

    qpts = np.array(inp['finite_qpts'])[1:]
    qpts = np.array([inp['finite_qpts'][0]])
    kpts, emf, eqp = read_eqp('bse_1/eqp.dat', inp)

    filepc = 'epc.h5'
    kpts_epc, qpts_epc = load_kq_epc(filepc)
    qlist = [kindex(qpt, qpts_epc, 12, 12, 1) for qpt in qpts]
    gmnvkq = load_epc(filepc, inp, qlist)
    print("Load epc info successfully!")

    lsame = check_kk_same(kpts, kpts_epc, 12, 12, 1)
    if (lsame):
        print("Kpts in eqp.dat and kpts in epc.h5 are the same!")
    else:
        print("Kpts in eqp.dat and kpts in epc.h5 are different!")

    Nkq = [12, 12, 1]
    H = constr_H(inp, eqp, gmnvkq, kpts, qpts, Nkq)

    eig, evec = np.linalg.eig(H)

    idx = np.argsort(eig.real)
    eig = eig[idx]
    evec = evec[:, idx]

    np.savetxt('evec.dat', evec, fmt='%20.5G')

    ndim = H.shape[0]

    fileig = 'EIG_qKK_m1-9_1PH.txt'
    filevec = 'EVEC_qKK_m1-9_1PH.txt'

    fileig = 'EIG_qG_m4-9_1PH.txt'
    filevec = 'EVEC_qG_m4-9_1PH.txt'

    eig_tdgw = np.loadtxt(fileig)
    data = np.loadtxt(filevec)
    evec_tdgw = data[:,0::2] + 1.0j * data[:,1::2]
    reorder = [3, 2, 1, 0, 4, 5, 6, 7]
    evec_tdgw = evec_tdgw.reshape(144, 8, ndim)[:, reorder, :].reshape(ndim, ndim)

    diff = np.sum(np.abs(eig - eig_tdgw)) / ndim
    print("Average difference of EIG: %.4G eV."%diff)
    
    diff = np.sum(np.abs(np.abs(evec) - np.abs(evec_tdgw))) / ndim
    print("Average difference of EVEC: %.4G."%diff)


def constr_H(inp, eqp, gmnvkq, kpts, qpts, Nkq):

    n1, n2, n3 = Nkq
    nk, nb = eqp.shape
    nq = qpts.shape[0]

    ndim = nk * nb
    H = np.zeros((ndim, ndim), dtype=complex)
    # g[in, im, iv, ik, iq] = g_nmv(k, q) = < nk+q | dV_qv | mk >
    # H[ii, jj] = en[ii, jj] + < in k+q | dV_qv | jm k >
    # Read g[in, jm, :, k, q] or g[jm, in, :, k+q, -q]^*

    phA = 1.0
    qread = np.zeros(nq, dtype=int)
    for iq in range(nq):
        iqp = kindex(-qpts[iq], qpts, n1, n2, n3)
        for ik in range(nk):
            kqpt = kpts[ik] + qpts[iq]
            ikq = kindex(kqpt, kpts, n1, n2, n3)
            for ib in range(nb):
                ii = ib + ikq * nb
                for jb in range(nb):
                    jj = jb + ik * nb

                    if (qread[iqp]==0):
                        H[ii, jj] += np.sum(gmnvkq[ib, jb, :, ik, iq] * phA)
                    else:
                        H[ii, jj] += np.sum(np.conj(gmnvkq[jb, ib, :, ikq, iqp]) * phA)

        qread[iq] = 1

    for ik in range(nk):
        for ib in range(nb):
            ii = ib + ik * nb
            H[ii, ii] += eqp[ik, ib]

    # Ry2eV = 13.605703976
    # H = H / Ry2eV

    return H


def check_kk_same(kpts_A, kpts_B, n1, n2, n3):

    lsame = True
    nkA = kpts_A.shape[0]
    for ik in range(nkA):
        jk = kindex(kpts_A[ik], kpts_B, n1, n2, n3)
        if (ik != jk): lsame = False

    return lsame


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


def load_kq_epc(filname):

    with h5py.File(filname, 'r') as f:
        kpts = f['header/kpts'][()]
        qpts = f['header/qpts'][()]

    return kpts, qpts


def load_epc(filname, inp, qlist):

    bmin = int(inp['band_index_min'])
    bmax = int(inp['band_index_max'])

    with h5py.File(filname, 'r') as f:
        dset = f['epmat/gmnvkq_real']
        epc_r = dset[bmin-1:bmax, bmin-1:bmax, :, :, qlist]
        dset = f['epmat/gmnvkq_imag']
        epc_i = dset[bmin-1:bmax, bmin-1:bmax, :, :, qlist]

    epc = epc_r + 1.0j * epc_i

    return epc


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


def read_inp(filname='dynamics.inp'):
    '''
    Read input parameters.

    Parameters:
    filname: string, coupling file.

    Returns: dictionary, the keys and values are all strings.
    '''

 
    inp = { # set default values
                 'calc_dyn' : '1',
               'n_finite_q' : '0',
                     'tmax' : '100.0',
                       'dt' : '0.1',
               'pulse_type' : '1',
                     'eamp' : '0.1',
                    'efreq' : '0.0',
                   'ewidth' : '1.0',
                   'tpulse' : '10.0',
                'iq_efield' : '1 1',
           'band_index_min' : '0',
           'band_index_max' : '0',
                'prop_algo' : '2',
              'finite_qpts' : []
          }

    text = [line for line in open(filname) if line.strip()]

    qpts_tag = 0
    for line in text:

        if (line[0]=='#'):
            continue

        if line.lower().startswith('begin'):
            qpts_tag = 1; continue
        if line.lower().startswith('end'):
            qpts_tag = 0; continue
        if (qpts_tag==1):
            qpt = [float(val) for val in line.split()[:3]]
            inp['finite_qpts'].append(qpt)

        if (qpts_tag==0):
            temp = line.split('#')[0].split(' ', 1)
            key = temp[0].strip().lower()
            value = temp[1].strip().strip('\'').strip('\"')
            inp[key] = value

    return inp


def load_inf(inp, filname='bsemat.h5'):

    inp['nk'] = read_h5file(filname, 'bse_header/kpoints/nk')
    inp['nvb'] = read_h5file(filname, 'bse_header/bands/nvb')
    inp['ncb'] = read_h5file(filname, 'bse_header/bands/ncb')
    inp['nspin'] = read_h5file(filname, 'bse_header/bands/ns')
    inp['efermi'] = read_h5file(filname, 'bse_header/params/efermi')
    inp['celvol'] = read_h5file(filname, 'mf_header/crystal/celvol')
    inp['bvec'] = read_h5file(filname, 'mf_header/crystal/bvec')
    inp['kpts'] = read_h5file(filname, 'bse_header/kpoints/kpts')
    # inp[''] = read_h5file(filname, '')


def read_h5file(filname, dset_name=""):

    dset = None
    with h5py.File(filname, 'r') as f:
        dset = f[dset_name][()]

    return dset



def error_exit(message='Unknown'):
    print('\n', "Error:", message, '\n')
    raise SystemExit


if __name__=='__main__':
    main()
