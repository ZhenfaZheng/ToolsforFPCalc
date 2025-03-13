#/usr/bin/env python

import math, os, h5py
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d
from functools import partial


def main():

    inp = read_inp()
    load_inf(inp, 'bse_1/bsemat.h5')
    A = np.array([
        [ 5.979600000, 0.000000000, 0.000000000],
        [-2.989800000, 5.178485504, 0.000000000],
        [ 0.000000000, 0.000000000, 37.794519772]
        ])

    filpol_1 = 'bse_1/polarization_rotate_G0_noeh_qG_m4-9_0.1ph_uc_eig.dat'
    filpol_2 = 'bse_1/polarization_rotate_G0_noeh_qG_m4-9_0.2ph_uc_eig.dat'
    filpol_3 = 'bse_1/polarization_rotate_G0_noeh_qG_m4-9_0.5ph_uc_eig.dat'
    filpol_4 = 'bse_1/polarization_rotate_G0_noeh_qG_m4-9_ph_uc_eig.dat'

    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol_1)
    pol = data[:,1] + 1.0j*data[:,2]
    weps_1 = calc_eps(pol, inp, A)

    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol_2)
    pol = data[:,1] + 1.0j*data[:,2]
    weps_2 = calc_eps(pol, inp, A)

    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol_3)
    pol = data[:,1] + 1.0j*data[:,2]
    weps_3 = calc_eps(pol, inp, A)

    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol_4)
    pol = data[:,1] + 1.0j*data[:,2]
    weps_4 = calc_eps(pol, inp, A)


    filabsp_0 = 'bse_1/absorption_noeh_noph-4b.dat'
    filabsp_1 = 'bse_1/absorption_noeh_qG_m4-9_0.1ph_uc_eig.dat'
    filabsp_2 = 'bse_1/absorption_noeh_qG_m4-9_0.2ph_uc_eig.dat'
    filabsp_3 = 'bse_1/absorption_noeh_qG_m4-9_0.5ph_uc_eig.dat'
    filabsp_4 = 'bse_1/absorption_noeh_qG_m4-9_ph_uc_eig.dat'

    weps_bse_0 = np.loadtxt(filabsp_0)[:, :2]
    weps_bse_1 = np.loadtxt(filabsp_1)[:, :2]
    weps_bse_2 = np.loadtxt(filabsp_2)[:, :2]
    weps_bse_3 = np.loadtxt(filabsp_3)[:, :2]
    weps_bse_4 = np.loadtxt(filabsp_4)[:, :2]

    wepss = [weps_1, weps_2, weps_3, weps_4]
    weps_bses = [weps_bse_1, weps_bse_2, weps_bse_3, weps_bse_4]
    tags = ['0.1-Phmod', '0.2-Phmod', '0.5-Phmod', '1-Phmod']

    plt_absp_cmp(weps_bse_0, wepss, weps_bses, tags, figname='Absp_comp.png')

    figname = "PolAbs_qG_m49_compare.png"
    en_absp_py = np.loadtxt('absp_ph_py.dat')
    plt_pol_abs(data[:,0], pol, weps_4[:,0], weps_4[:,1], weps_bse_4, en_absp_py, figname=figname)


    figname = "PolAbs_qKK_m19_compare.png"
    filpol = 'bse_1/polarization_rotate_G0_noeh_qKK_m1-9_ph_uc_eig.dat'
    filabsp = 'bse_1/absorption_noeh_qKK_m1-9_ph_uc_eig.dat'
    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol)
    pol = data[:,1] + 1.0j*data[:,2]
    weps = calc_eps(pol, inp, A)
    weps_bse = np.loadtxt(filabsp)[:, :2]
    en_absp_py = np.loadtxt('absp_qKK_m1-9_ph_py.dat')
    plt_pol_abs(data[:,0], pol, weps[:,0], weps[:,1], weps_bse, en_absp_py, figname=figname)

    figname = "PolAbs_qG_m49_eh_compare.png"
    filpol = 'bse_1/polarization_rotate_G0_eh_qG_m4-9_ph_uc_eig_vmt_gbar.dat'
    filabsp = 'bse_1/absorption_eh_qG_m4-9_ph_uc_eig.dat'
    inp['eamp'] = '0.05'
    data = np.loadtxt(filpol)
    pol = data[:,1] + 1.0j*data[:,2]
    weps = calc_eps(pol, inp, A)
    weps_bse = np.loadtxt(filabsp)[:, :2]
    filabsp = 'bse_1/absorption_eh_qG_m4-9_ph_uc_eig_vmt.dat'
    en_absp_py = np.loadtxt(filabsp)[:, :2]
    plt_pol_abs(data[:,0], pol, weps[:,0], weps[:,1], weps_bse, en_absp_py, figname=figname)



def plt_pol_abs(time, pol, wlist, epsw, weps_bse=None, en_absp_py=None, figname='PolAbs.png'):

    figsize_x = 4.8
    figsize_y = 4.0 # in inches
    fig, axes = plt.subplots(2,1)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    ax = axes[0]
    cls = ['r', 'b']
    Ry2fs = 0.0484
    x = time * Ry2fs
    ax.plot(x, pol.real, 'r', lw=1.0, label='real part')
    ax.plot(x, pol.imag, 'b', lw=1.0, label='imaginary part')
    ax.legend(loc=1, ncol=2)

    xmin = 0 if (x.min() < 0.1) else x.min()
    xmax = x.max()

    ax.set_xlim(xmin, xmax)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Polarization')

    wmin = 1 # eV
    wmax = 8 # eV
    index = np.where(wlist<wmax)
    index = np.where((wlist<wmax) & (wlist>wmin))
    wlist = wlist[index]; epsw = epsw[index]

    sort = np.argsort(wlist)
    wlist = wlist[sort]; epsw = epsw[sort]

    ax = axes[1]
    # num = int(len(wlist)/10)
    # epsw[:num] = 0.0
    # ax.plot([wmin, wmax], [0, 0], 'k', ls='--', lw=0.5)
    ax.plot(wlist, epsw, 'gray', lw=1.5, label='TD-aGW-ph')

    if(weps_bse is not None):
        index = np.where(weps_bse[:,0] < wmax)[0]
        # ax.plot(weps_bse[index,0], weps_bse[index,1], 'g', ls=':', lw=2.2, label='GW-BSE')
        ax.plot(weps_bse[index,0], weps_bse[index,1], 'g', ls='-', lw=1.0, label='GW-BSE')
        # ax.legend(loc=9)

    if(en_absp_py is not None):
        index = np.where(en_absp_py[:,0] < wmax)[0]
        scale = epsw.max()/en_absp_py[:,1].max()
        ax.plot(en_absp_py[index,0], en_absp_py[index,1]*scale, 'orange', ls='--', lw=0.7, label='Absp_1st_order')
        # ax.plot(en_absp_py[index,0], en_absp_py[index,1]*scale, '', ls=':', lw=2.2, label='Absp_1st_order')
        ax.legend(loc=1)

    # ax.set_xlim(0, wmax)
    ax.set_xlim(wmin, wmax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Absorption')
    # ax.set_ylabel('Amplitude')

    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    print("\n%s has been saved."%figname)



def plt_absp_cmp(weps_bse_0, wepss, weps_bses, tags, figname='Absp_comp.png'):

    npair = len(wepss)

    figsize_x = 4.8 * 2
    figsize_y = 2.4 * npair # in inches
    fig, axes = plt.subplots(npair,2)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    wmin = 2
    wmax = 8

    for ipair in range(npair):

        ax = axes[ipair, 0]

        ax.plot([wmin, wmax], [0, 0], 'k', ls='--', lw=0.5)

        weps = weps_bse_0
        index = np.where((weps[:,0]<wmax) & (weps[:,0]>wmin))[0]
        ax.plot(weps[index,0], weps[index,1], 'gray', lw=1.0, label='GW-BSE-NOPH')

        weps = weps_bses[ipair]
        index = np.where((weps[:,0]<wmax) & (weps[:,0]>wmin))[0]
        ax.plot(weps[index,0], weps[index,1], 'r', lw=1.0, label='GW-BSE')

        weps = wepss[ipair]
        index = np.where((weps[:,0]<wmax) & (weps[:,0]>wmin))[0]
        ax.plot(weps[index,0], weps[index,1], 'b', lw=1.0, label='TDaGW-PH')

        ax.set_title(tags[ipair])
        ax.set_xlim(wmin, wmax)
        ax.set_xlabel('Energy (eV)')


        ax = axes[ipair, 1]

        wlist = np.arange(wmin, wmax, 0.01)
        ax.plot([wmin, wmax], [0, 0], 'k', ls='--', lw=0.5)

        eps_diff = calc_diff_absp(weps_bse_0, weps_bses[ipair], wlist)
        ax.plot(wlist, eps_diff, 'r', lw=1.0, label='GW-BSE')

        eps_diff = calc_diff_absp(weps_bse_0, wepss[ipair], wlist)
        ax.plot(wlist, eps_diff, 'b', lw=1.0, label='TDaGW-PH')

        ax.set_title(tags[ipair])
        ax.set_xlim(wmin, wmax)
        ax.set_xlabel('Energy (eV)')


    plt.tight_layout()
    plt.savefig(figname, dpi=300)
    print("\n%s has been saved."%figname)


def calc_diff_absp(weps_1, weps_2, x):

    interp_y1 = interp1d(weps_1[:,0], weps_1[:,1], kind='linear', fill_value="extrapolate")
    interp_y2 = interp1d(weps_2[:,0], weps_2[:,1], kind='linear', fill_value="extrapolate")

    y1_interp = interp_y1(x)
    y2_interp = interp_y2(x)

    y_diff = y2_interp - y1_interp

    return y_diff


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



def calc_eps(pol, inp, A, fileps='eps.dat'):

    nk1 = 12
    nk2 = 12
    nk3 = 1
    nk = nk1 * nk2 * nk3

    cell_vol = np.dot(A[0], np.cross(A[1], A[2]))

    eamp = float(inp['eamp'])
    efreq = float(inp['efreq'])
    ewidth = float(inp['ewidth'])
    tpulse = float(inp['tpulse'])
    dt = float(inp['dt'])

    it_start = int( tpulse / dt )

    dlen = 100000
    dlen = max(dlen, len(pol))

    Ry2eV = 13.605703976
    gamma = 0.05/Ry2eV

    # processing polarization data
    avg = np.average(pol)
    Pfull = np.zeros(dlen, dtype=complex)
    Pfull[:len(pol)] = (pol[:] - avg) * 4.0 * np.pi * np.sqrt(2.0) / nk / cell_vol
    P_decay = np.zeros(dlen, dtype=complex)
    for it in range(dlen):
        P_decay[it] = Pfull[it] * decay_exp(it*dt, tpulse, gamma)

    P_decay = P_decay[it_start:]
    # fft of polarization
    Pw = np.fft.ifft(P_decay, axis=0)
    # Pw = np.fft.irfft(P_decay.real, axis=0)

    # fft of external field (light)
    tlist = np.arange(dlen, dtype=float) * dt + 2.0 * dt
    Et = eamp * gaussian(tlist, tpulse, ewidth) * np.cos(efreq * tlist)
    for it in range(dlen):
        Et[it] = Et[it] * decay_exp(it*dt, tpulse, gamma)
    Et = Et[it_start:]
    Ew = np.fft.ifft(Et, axis=0)
    # Ew = np.fft.irfft(Et.real, axis=0)
    if (eamp==0.0): Ew[:] = 1.0

    # calculate frequency list
    wlist = np.fft.fftfreq(dlen-it_start,dt)*np.pi*2.0*Ry2eV

    # ind = np.argwhere(Ew==Ew.max())[0]
    # print("Photo energy: %.2f eV\n"%(-wlist[ind]))

    epsw = np.sqrt(2.0) * Pw / Ew
    epsw *= 2.0

    weps = np.vstack((wlist, epsw.imag)).T

    return weps


def gaussian(x, mu, sigma):

    y = np.exp(-0.5 * ((x - mu) / sigma) **2)
    y /= sigma * np.sqrt(2.0 * np.pi)

    return y

def decay_exp(t,tpulse,gamma):

    # exponential decay -> Lorentzian broadening
    factor = 1
    if t >= tpulse:
       factor = np.exp(-gamma*(t-tpulse))
   
    return factor


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

