#!/usr/bin/env python

import numpy as np

def main():

    kpts = np.loadtxt('eph.kpt', skiprows=1)[:,:3]
    qpts_tot = kgenerate(450)
    # qpts = kqmap(kpts, qpts_tot)
    qpts = kq_match(kpts, kpts, qpts_tot, 450, 450, 1)
    savekpts('eph.qpt.new', qpts)


def kgenerate(numstep):

    kstep = 1.0 / numstep
    kcenter = np.array([0.0, 0.0, 0.0])
    kpoints = []

    for i in range(numstep):
        offsetx = kstep * i
        kx = kcenter[0] + offsetx
        for j in range(numstep):
            offsety = kstep * j
            ky = kcenter[1] + offsety
            kz = kcenter[2]
            kpoints.append([kx,ky,kz])

    return np.array(kpoints)


def kq_match(kpts_A, kpts_B, qpts_tot, n1, n2, n3):
    """
    Find all q-points in `qpts_tot` that satisfy qpt = kpt_A - kpt_B
    for any kpt_A in `kpts_A` and any kpt_B in `kpts_B`.

    Parameters:
        kpts_A (ndarray): Array of k-points from set A, shape (nk_A, 3).
        kpts_B (ndarray): Array of k-points from set B, shape (nk_B, 3).
        qpts_tot (ndarray): Array of q-points to search in, shape (nq, 3).

    Returns:
        qpts: Array of q-points in `qpts_tot` that match the condition, shape(nq, 3).
    """


    nk_A = kpts_A.shape[0]
    nk_B = kpts_B.shape[0]
    nq = qpts_tot.shape[0]

    k_index_A = np.zeros((nk_A, 3))
    k_index_B = np.zeros((nk_B, 3))
    q_tag = np.zeros(nq, dtype=int)
    q_loc = np.zeros((n1, n2, n3), dtype=int)
    q_loc += -1

    for ik in range(nk_A):
        i1 = round(kpts_A[ik,0] * n1)
        i2 = round(kpts_A[ik,1] * n2)
        i3 = round(kpts_A[ik,2] * n3)
        k_index_A[ik] = np.array([i1, i2, i3])

    for ik in range(nk_B):
        i1 = round(kpts_B[ik,0] * n1)
        i2 = round(kpts_B[ik,1] * n2)
        i3 = round(kpts_B[ik,2] * n3)
        k_index_B[ik] = np.array([i1, i2, i3])

    for iq in range(nq):
        i1 = round(qpts_tot[iq,0] * n1)
        i2 = round(qpts_tot[iq,1] * n2)
        i3 = round(qpts_tot[iq,2] * n3)
        q_loc[i1, i2, i3] = iq

    for ik in range(nk_A):
        for jk in range(nk_B):
            dk_1 = int((k_index_A[ik,0] - k_index_B[jk,0]) % n1)
            dk_2 = int((k_index_A[ik,1] - k_index_B[jk,1]) % n2)
            dk_3 = int((k_index_A[ik,2] - k_index_B[jk,2]) % n3)
            iq = q_loc[dk_1, dk_2, dk_3]
            if (iq>=0): q_tag[iq] = 1

    qpts= qpts_tot[q_tag==1]

    return qpts


def kqmap(kpts, qpts_tot, mapax='xy'):

    std = 1.0e-5

    nks = kpts.shape[0]
    nqs_tot = qpts_tot.shape[0]

    axtag = np.zeros(3)
    if 'x' in mapax: axtag[0] = 1
    if 'y' in mapax: axtag[1] = 1
    if 'z' in mapax: axtag[2] = 1

    qtag = np.zeros((nqs_tot, 3))

    for ik in range(nks):
        if (ik%5==0): print(ik)
        for jk in range(nks):

            dk = kpts[ik] - kpts[jk]
            qtag[np.sum(qtag, axis=1)<3, :] = 0
            qtag[:,axtag==0] = 1.0

            for ia in range(3):
                if (axtag[ia]==0): continue
                if (dk[ia]<0): dk[ia] += 1.0
                dkqa = qpts_tot[:,ia] - dk[ia]
                qtag[np.abs(dkqa)<std, ia] = 1.0
            '''
            for iq in range(nqs_tot):
                dkq = np.abs(qpts_tot[iq] - dk)
                if (dkq[0]<std and dkq[1]<std and dkq[1]<std):
                    qtag[iq] = 1
                    break
            '''

    qpts = qpts_tot[np.sum(qtag, axis=1)==3, :]
    return qpts


def savekpts(filname, kpoints, weight=None):
    
    nks = kpoints.shape[0]

    if weight is None:
        weight = np.ones((nks,1))

    header = '%d'%nks
    data = np.append(kpoints, weight, axis=1)

    np.savetxt(filname, data,
               header ='%d'%nks, comments='',
               fmt='%15.8f%15.8f%15.8f%6d')


if __name__=='__main__':
    main()
