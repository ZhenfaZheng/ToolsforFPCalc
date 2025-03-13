#/usr/bin/env python3

import struct
import math, os, h5py
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt


def main():

    filvmt = 'bse_1/vmtxel.A'
    vmt, params = read_vmt_bin(filvmt)

    tag = 'qKK_m1-9'
    tag = 'qG_m4-9'

    eig_ph = np.loadtxt('EIG_qKK_m1-9_1PH.txt')
    data = np.loadtxt('EVEC_qKK_m1-9_1PH.txt')
    eig_ph = np.loadtxt('EIG_qG_m4-9_1PH.txt')
    data = np.loadtxt('EVEC_qG_m4-9_1PH.txt')
    evec_ph = data[:,0::2] + 1.0j * data[:,1::2]

    fileqp = 'bse_1/eqp.dat'
    kpts, emf, eqp = read_eqp(fileqp)
    eqp = eqp[:,3:] # select 4vb and 4cb in the 11 bands.
    nk, ncb, nvb, ns, opr = params
    enlist = eqp_to_elist_flat(eqp, nvb, ncb)

    vmt_all = construct_vmt_all(vmt, params)
    vmt_ph = np.matmul(evec_ph.conj().T, np.matmul(vmt_all, evec_ph))

    # en_absp = calc_absp(enlist, vmt_all)
    # plt_spectrum(en_absp, figname='Absp_noph_py.png')
    en_absp = calc_absp(eig_ph, vmt_ph)
    plt_spectrum(en_absp, figname='Absp_%s_ph_py.png'%tag)
    
    print(np.average(np.abs(vmt_ph[:576, :576])))
    print(np.average(np.abs(vmt_ph[:576, 576:])))
    print(np.average(np.abs(vmt_ph[576:, :576])))
    print(np.average(np.abs(vmt_ph[576:, 576:])))
    # for i in range(1000):
    #     for j in range(5):
    #         if (np.abs(vmt_ph[i, j])<0.0001): continue
    #         print(i, j, vmt_ph[i, j], vmt_ph[j, i])

    np.savetxt('absp_%s_ph_py.dat'%tag, en_absp, fmt='%15.9f')

    print("\nDONE!\n")


def plt_spectrum(en_y, figname='Spectrum.png'):

    figsize_x = 6.4
    figsize_y = 3.2 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    ax.plot(en_y[:,0], en_y[:,1], 'r', lw=1.0)

    emin = en_y[:,0].min()
    emax = en_y[:,0].max()
    emax = 8
    ax.set_xlim(emin, emax)
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Amplitude (a.u.)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)



def calc_absp(enlist, vmt, sigma=0.05):

    demin = 0
    demax = enlist.max() - enlist.min()
    dex = np.linspace(demin, demax, 1000)
    absp = np.zeros_like(dex)

    nbas = enlist.shape[0]
    for ibas in range(nbas):
        for jbas in range(nbas):
            Amp = np.abs(vmt[ibas, jbas])**2
            if (Amp==0): continue
            de = enlist[ibas] - enlist[jbas]
            # absp += Amp * gaussian(dex, de, sigma)
            absp += Amp * lorentzian(dex, de, sigma)

    return np.vstack((dex, absp)).T

def eqp_to_elist_flat(eqp, nvb, ncb):

    nk, nb = eqp.shape
    if (nb != (nvb+ncb)): error_exit('The number of bands dost not match!')
    enlist = np.zeros(nk*nb)

    for ik in range(nk):
        for ib in range(nvb):
            ibas = nvb - ib - 1 + ik * nb
            enlist[ibas] = eqp[ik, ib]
        for jb in range(ncb):
            jbas = nvb + jb + ik * nb
            enlist[jbas] = eqp[ik, nvb+jb]

    return enlist


def construct_vmt_all(vmt, params):

    nk, ncb, nvb, ns, opr = params

    nb = nvb + ncb
    nbas = nk * nb

    vmt_all = np.zeros((nbas, nbas), dtype=complex)

    for ik in range(nk):
        for icb in range(ncb):
            ibas = icb + nvb + ik * nb
            for jvb in range(nvb):
                jbas = jvb + ik * nb
                bse_index = jvb + (icb + ik*ncb)*nvb
                vmt_all[ibas, jbas] = vmt[bse_index]

    return vmt_all


def write_vmt_bin(filvmt, params, data, is_complex=True):

    dtype = np.complex128 if is_complex else np.float64

    with open(filvmt, 'wb') as f:

        header_data = struct.pack('5i', *params)
        record_length = len(header_data)
        f.write(struct.pack('i', record_length))  # Write record length
        f.write(header_data)                      # Write data
        f.write(struct.pack('i', record_length))  # Write record length again

        data_formated = data.astype(dtype)
        record_length = data_formated.nbytes  # Total bytes
        f.write(struct.pack('i', record_length))  # Write record length again
        f.write(data_formated)
        f.write(struct.pack('i', record_length))  # Write record length again



def read_vmt_bin(filvmt, is_complex=True):

    dtype = np.complex128 if is_complex else np.float64

    with open(filvmt, 'rb') as f:

        # In the binary file, data is stored as: [record_length] data [record_length]
        # we need to read the record length first

        f.read(4)
        # length = struct.unpack('i', f.read(4))
        nk, nband, mband, ns, opr = struct.unpack('5i', f.read(20))
        params = (nk, nband, mband, ns, opr)
        f.read(4)

        f.read(4)
        nmat = nk * nband * mband * ns
        data = np.fromfile(f, dtype=dtype, count=nmat)
        f.read(4)

    return data, params


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


def gaussian(x, mu, sigma):

    return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((x - mu) / sigma)**2)

def lorentzian(x, x0, gamma):

    return (gamma / np.pi) / ((x - x0)**2 + gamma**2)

def error_exit(message='Unknown'):
    print('\n', "Error:", message, '\n')
    raise SystemExit



if __name__=='__main__':
    main()
