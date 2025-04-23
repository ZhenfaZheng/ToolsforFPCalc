import h5py
import numpy as np
import os

def main():

    data = np.loadtxt('EVEC_tdagw-ph.txt')
    evec_ph = data[:,0::2] + 1.0j * data[:,1::2]

    nb = 8
    kb_index = get_ik_of_eig_states_with_ph(evec_ph, nb)

    rotate_wfn('wfn.h5', evec_ph, kb_index, fil_out='wfn_ph.h5')


def rotate_wfn(fil_in, evec_ph, kb_index, fil_out='wfn_ph.h5'):

    # Copy file1 to output using os.system
    os.system(f'cp "{fil_in}" "{fil_out}"')  # Linux/Mac
    # For Windows use: os.system(f'copy "{fil_in}" "{fil_out}"')

    with h5py.File(fil_in, 'r') as f:

        gvecs = f['/wfns/gvecs'][:]
        nGk = f['/mf_header/kpoints/ngk'][:]
        coeffs = f['/wfns/coeffs'][:]

    kG_index = np.hstack((np.array([0]), np.cumsum(nGk)))

    coeffs_cplx = coeffs[..., 0] + 1j * coeffs[..., 1]
    coeffs_new = coeffs_cplx * 1.0

    # from tdagw dynamics.inp
    bmin = 18
    bmax = 25
    nb = bmax - bmin + 1
    nk = nGk.shape[0]

    for ik in range(nk):

        ikG = kG_index[ik]
        jkG = kG_index[ik+1]

        for ib in range(nb):

            coeffs_new[ib+bmin-1, :, ikG:jkG] = 0.0

            idim = kb_index[ik, ib]
            # evec_ph[ik*nb:ik*nb+nb, idim] = 0.0
            # evec_ph[ik*nb+ib, idim] = 1.0
            # print(ik, ib, evec_ph[ik*nb:ik*nb+nb, idim])
            reorder = [3, 2, 1, 0, 4, 5, 6, 7]
            for jb in range(nb):
                jbp = reorder[jb]
                coeffs_new[ib+bmin-1, :, ikG:jkG] += coeffs_cplx[jb+bmin-1, :, ikG:jkG] * evec_ph[ik*nb+jbp, idim]

            norm_factor = np.sum(np.abs(coeffs_new[ib+bmin-1, :, ikG:jkG])**2)
            print(ik, ib, norm_factor)
            if (norm_factor==0): norm_factor = 1.0
            coeffs_new[ib+bmin-1, :, ikG:jkG] /= np.sqrt(norm_factor)

    coeffs[..., 0] = coeffs_new.real
    coeffs[..., 1] = coeffs_new.imag

    with h5py.File(fil_out, 'r+') as f_out:
        del f_out['/wfns/coeffs']
        f_out.create_dataset('/wfns/coeffs', data=coeffs, dtype=coeffs.dtype)


def get_ik_of_eig_states_with_ph(evec_ph, nb):
    
    nstates = evec_ph.shape[0]
    nk = int(nstates / nb)
    kb_index = np.zeros((nk, nb), dtype=int)

    egstates_ph_ik = []
    evec_ph_abs = np.abs(evec_ph)
    for i in range(nstates):
        index = np.argmax(evec_ph_abs[:,i])
        ik = index // nb # floor of index/nb
        egstates_ph_ik.append(ik)

    egstates_ph_ik = np.array(egstates_ph_ik)

    for ik in range(nk):
        s_index = np.argwhere(egstates_ph_ik==ik)[:,0]
        kb_index[ik,:] = s_index

    # reorder = [3, 2, 1, 0, 4, 5, 6, 7]
    # kb_index = kb_index[:, reorder]
    # print(kb_index[0])
    return kb_index


if __name__=='__main__':
    main()
