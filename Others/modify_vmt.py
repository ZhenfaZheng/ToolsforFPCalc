#/usr/bin/env python3

import struct
import numpy as np


def main():

    filvmt = 'vmtxel_origin'
    vmt_origin, params = read_vmt_bin(filvmt)

    filvmt = 'vmtxel_uc'
    vmt, params = read_vmt_bin(filvmt)
    vmt_all = construct_vmt_all(vmt, params)

    data = np.loadtxt('EVEC_uc_ph.txt')
    evec_ph = data[:,0::2] + 1.0j * data[:,1::2]
    vmt_ph = np.matmul(evec_ph.conj().T, np.matmul(vmt_all, evec_ph))
    # the oder of states in evec and vmt_ph is in energy order

    print(np.average(np.abs(vmt_all[:576, :576])))
    print(np.average(np.abs(vmt_all[:576, 576:])))
    print(np.average(np.abs(vmt_all[576:, :576])))
    print(np.average(np.abs(vmt_all[576:, 576:])))

    print(np.average(np.abs(vmt_ph[:576, :576])))
    print(np.average(np.abs(vmt_ph[:576, 576:])))
    print(np.average(np.abs(vmt_ph[576:, :576])))
    print(np.average(np.abs(vmt_ph[576:, 576:])))

    # for i in range(1000):
    #     for j in range(5):
    #         if (np.abs(vmt_ph[i, j])<0.0001): continue
    #         print(i, j, vmt_ph[i, j], vmt_ph[j, i])

    kpts, emf, eqp = read_eqp('eqp.dat')
    kpts, emf, eqp = read_eqp('eqp_scph.dat')
    eig_ph = np.loadtxt('EIG_uc_ph.txt')
    egstates_ph_ik = get_ik_of_eig_states_with_ph(evec_ph, params)
    # np.savetxt('xx', egstates_ph_ik, fmt='%d')
    vmt_new = construct_vmt_flatten(vmt_ph, params, egstates_ph_ik, eig_ph, eqp, vmt_origin)

    print(np.average(np.abs(vmt_all)))
    print(np.average(np.abs(vmt_ph)))
    print(np.average(np.abs(vmt)))
    print(np.average(np.abs(vmt_new)))

    filvmt = 'vmtxel'
    write_vmt_bin(filvmt, params, vmt_new)


def get_ik_of_eig_states_with_ph(evec_ph, params):
    
    nstates = evec_ph.shape[0]
    nk, ncb, nvb, ns, opr = params
    nb = nvb + ncb

    egstates_ph_ik = []
    evec_ph_abs = np.abs(evec_ph)
    for i in range(nstates):
        index = np.argmax(evec_ph_abs[:,i])
        ik = index // nb # floor of index/nb
        egstates_ph_ik.append(ik)

    return np.array(egstates_ph_ik)


def construct_vmt_flatten(vmt_ph, params, egstates_ph_ik, eig_ph, eqp , vmt_origin):
    # construct flatten vmt in the same order of eqp

    nk, ncb, nvb, ns, opr = params

    nb = nvb + ncb
    nbas = nk * nb

    vmt = np.zeros(nk*ncb*nvb, dtype=complex)

    vb_r_index = range(nvb-1, -1, -1)
    eqp[:, :nvb] = eqp[:, vb_r_index]
    # sorted_indices = np.argsort(np.argsort(eqp.flatten())).reshape(eqp.shape)
    # sorted_indices[ik, ib] = number of states in energy order

    sorted_indices = np.zeros((nk, nb), dtype=int)
    for ik in range(nk):
        indices_ik = np.where(egstates_ph_ik==ik)[0]
        eig_ph_k = eig_ph[indices_ik]
        sorted_indices[ik,:] = indices_ik[np.argsort(eig_ph_k)]
    sorted_indices[:, :nvb] = sorted_indices[:, vb_r_index]
    print(sorted_indices[52])

    for ik in range(nk):
        for icb in range(ncb):
            ib = icb + nvb
            ibas = sorted_indices[ik, ib]
            for jvb in range(nvb):
                jb = jvb
                jbas = sorted_indices[ik, jb]
                bse_index = jvb + (icb + ik*ncb)*nvb
                vmt[bse_index] = vmt_ph[ibas, jbas]
                scale = np.abs(vmt[bse_index]) / np.abs(vmt_origin[bse_index])
                vmt[bse_index] = vmt_origin[bse_index] * scale

    return vmt



def construct_vmt_flatten_old(vmt_ph, params, eqp):
    # construct flatten vmt in the same order of eqp

    nk, ncb, nvb, ns, opr = params

    nb = nvb + ncb
    nbas = nk * nb

    vmt = np.zeros(nk*ncb*nvb, dtype=complex)

    vb_r_index = range(nvb-1, -1, -1)
    eqp[:, :nvb] = eqp[:, vb_r_index]
    sorted_indices = np.argsort(np.argsort(eqp.flatten())).reshape(eqp.shape)
    # sorted_indices[ik, ib] = number of states in energy order

    # eig_ph = np.loadtxt('EIG_uc_ph.txt')

    for ik in range(nk):
        for icb in range(ncb):
            ib = icb + nvb
            ibas = sorted_indices[ik, ib]
            # print(ibas, ik, ib, eqp[ik, ib], eig_ph[ibas])
            for jvb in range(nvb):
                jb = jvb
                jbas = sorted_indices[ik, jb]
                bse_index = jvb + (icb + ik*ncb)*nvb
                vmt[bse_index] = vmt_ph[ibas, jbas]
                # if (ik%10==0): print(ibas, jbas, np.abs(vmt_ph[ibas, jbas]), np.abs(vmt_ph[ibas-2:ibas+2, jbas-2:jbas+2]))

    return vmt


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


if __name__=='__main__':
    main()

