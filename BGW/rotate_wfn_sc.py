import h5py
import numpy as np
import os

def main():

    data = np.loadtxt('EVEC_tdagw-ph.txt')
    evec_ph = data[:,0::2] + 1.0j * data[:,1::2]

    nb_uc = 8
    nb_sc = 24
    k_uc_sc_map = get_kk_map()
    kb_index = get_ikb_of_eig_states_with_ph(evec_ph, nb_uc, nb_sc, k_uc_sc_map)

    # rotate_wfn('wfn.h5', evec_ph, kb_index, fil_out='wfn_ph.h5')
    rotate_wfn_sc('wfn_uc.h5', 'wfn_sc.h5', evec_ph, kb_index, fil_out='wfn_ph.h5')


def rotate_wfn_sc(f_ucw, f_scw, evec_ph, kb_index, fil_out='wfn_ph.h5'):

    # from tdagw dynamics.inp
    bmin = 18
    bmax = 25
    nb_uc = bmax - bmin + 1
    nb_sc = nb_uc * 3

    with h5py.File(f_ucw, 'r') as f:

        kpts_uc = f['/mf_header/kpoints/rk'][:]
        nGk_uc = f['/mf_header/kpoints/ngk'][:]
        gvecs_uc = f['/wfns/gvecs'][:]
        coeffs_dset = f['/wfns/coeffs']
        coeffs_uc = coeffs_dset[bmin-1:bmax, :, :, :]

    nk_uc = kpts_uc.shape[0]
    kG_index_uc = np.hstack((np.array([0]), np.cumsum(nGk_uc)))
    kgvecs_uc = 1.0 * gvecs_uc
    for ik in range(nk_uc):
        ikG = kG_index_uc[ik]
        jkG = kG_index_uc[ik+1]
        kgvecs_uc[ikG:jkG, 0]  += kpts_uc[ik, 0]
        kgvecs_uc[ikG:jkG, 1]  += kpts_uc[ik, 1]
        kgvecs_uc[ikG:jkG, 2]  += kpts_uc[ik, 2]


    with h5py.File(f_scw, 'r') as f:

        kpts_sc = f['/mf_header/kpoints/rk'][:]
        nGk_sc = f['/mf_header/kpoints/ngk'][:]
        gvecs_sc = f['/wfns/gvecs'][:]
        coeffs_dset = f['/wfns/coeffs']
        coeffs_sc = coeffs_dset[(bmin-1)*3:bmax*3, :, :, :]

    nk_sc = kpts_sc.shape[0]
    kG_index_sc = np.hstack((np.array([0]), np.cumsum(nGk_sc)))
    kgvecs_sc = 1.0 * gvecs_sc
    for ik in range(nk_sc):
        ikG = kG_index_sc[ik]
        jkG = kG_index_sc[ik+1]
        kgvecs_sc[ikG:jkG, 0]  += kpts_sc[ik, 0]
        kgvecs_sc[ikG:jkG, 1]  += kpts_sc[ik, 1]
        kgvecs_sc[ikG:jkG, 2]  += kpts_sc[ik, 2]

    # gmap[nkg_sc] --> loc in kgvecs_uc
    gmap = kG_match(kgvecs_uc, kgvecs_sc, n1=4, n2=12, n3=1)

    # print(gmap[nGk_sc[0]-3:nGk_sc[0]+300])
    # print(gmap[-100:])
    # print(kgvecs_sc[nGk_sc[0]-2:nGk_sc[0]+100,:])
    # print('gvecs_sc')
    # print(gvecs_sc[nGk_sc[0]-100:nGk_sc[0]+300,:])
    # print('coeffs_sc')
    # print(coeffs_sc[0, 0, :20, :])

    c_new_sc = coeffs_sc[..., 0] + 1j * coeffs_sc[..., 1]
    c_cplx_uc = coeffs_uc[..., 0] + 1j * coeffs_uc[..., 1]
    c_temp_uc = c_cplx_uc * 0.0

    ndim = nb_uc * nk_uc
    reorder = np.array([3, 2, 1, 0, 4, 5, 6, 7], dtype=int)

    for ik in range(nk_sc):
        ikG = kG_index_sc[ik]
        jkG = kG_index_sc[ik+1]
        for ib in range(nb_sc):

            idim = kb_index[ik, ib]
            c_temp_uc[:] = 0.0

            for jdim in range(ndim):

                jk_uc = jdim // nb_uc # floor of jdim/nb
                jb_uc = int(jdim%nb_uc)
                jbp_uc = reorder[jb_uc]

                ikG_uc = kG_index_uc[jk_uc]
                jkG_uc = kG_index_uc[jk_uc+1]

                if (np.abs(evec_ph[jdim, idim]) < 1.0e-10):
                    continue
                # else:
                #     print(ik, ib)
                #     print(gmap[ikG:ikG+10])
                #     print(ikG_uc)

                c_temp_uc[jbp_uc, :, ikG_uc:jkG_uc] = c_cplx_uc[jbp_uc, :, ikG_uc:jkG_uc] * evec_ph[jdim, idim]

            c_new_sc[ib, :, ikG:jkG] = np.sum(c_temp_uc[:, :, gmap[ikG:jkG]], axis=0)
            norm_factor = np.sum(np.abs(c_new_sc[ib, :, ikG:jkG])**2)
            print(ik, ib, norm_factor)
            if (norm_factor==0): norm_factor = 1.0
            c_new_sc[ib, :, ikG:jkG] /= np.sqrt(norm_factor)

    coeffs_sc[..., 0] = c_new_sc.real
    coeffs_sc[..., 1] = c_new_sc.imag

    # Copy file1 to output using os.system
    os.system(f'cp "{f_scw}" "{fil_out}"')  # Linux/Mac
    # For Windows use: os.system(f'copy "{f_scw}" "{fil_out}"')

    with h5py.File(fil_out, 'r+') as f:
        coeffs_dset = f['/wfns/coeffs']
        coeffs_dset[(bmin-1)*3:bmax*3, :, :, :] = coeffs_sc


def kG_match(kgvecs_uc, kgvecs_sc, n1=4, n2=12, n3=1):

    '''
    Match kgvecs_uc * [3,1,1] = kgvecs_sc
    return gmap[nkg_sc] --> loc in kgvecs_uc
    '''

    nkg_uc = kgvecs_uc.shape[0]
    nkg_sc = kgvecs_sc.shape[0]

    kgvecs_uc[:,0] *= 3

    grange_1 = int( (kgvecs_uc[:,0].max() - kgvecs_uc[:,0].min()) * n1 + 2 )
    grange_2 = int( (kgvecs_uc[:,1].max() - kgvecs_uc[:,1].min()) * n2 + 2 )
    grange_3 = int( (kgvecs_uc[:,2].max() - kgvecs_uc[:,2].min()) * n3 + 2 )

    g_loc_uc = np.zeros((grange_1, grange_2, grange_3), dtype=int)
    gmap = np.zeros(nkg_sc, dtype=int)

    for ikg in range(nkg_uc):
        i1 = int( round(kgvecs_uc[ikg,0] * n1) )
        i2 = int( round(kgvecs_uc[ikg,1] * n2) )
        i3 = int( round(kgvecs_uc[ikg,2] * n3) )
        g_loc_uc[i1, i2, i3] = ikg

    for ikg in range(nkg_sc):
        i1 = int( round(kgvecs_sc[ikg,0] * n1) )
        i2 = int( round(kgvecs_sc[ikg,1] * n2) )
        i3 = int( round(kgvecs_sc[ikg,2] * n3) )
        gmap[ikg] = g_loc_uc[i1, i2, i3]

    return gmap


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


def get_kk_map():
    
    nk_uc = 144
    k_uc_sc_map = np.zeros(nk_uc, dtype=int)

    for ix in range(12):
        for iy in range(12):
            ii = ix + 12*iy
            k_uc_sc_map[ii] = ix%4 + 4*iy

    return k_uc_sc_map


def get_ikb_of_eig_states_with_ph(evec_ph, nb_uc, nb_sc, k_uc_sc_map):
    
    nstates = evec_ph.shape[0]
    nk_sc = int(nstates / nb_sc)
    kb_index = np.zeros((nk_sc, nb_sc), dtype=int)

    egstates_ph_ik = []
    evec_ph_abs = np.abs(evec_ph)
    for i in range(nstates):
        index = np.argmax(evec_ph_abs[:,i])
        ik_uc = index // nb_uc # floor of index/nb
        ik = k_uc_sc_map[ik_uc]
        egstates_ph_ik.append(ik)

    egstates_ph_ik = np.array(egstates_ph_ik)

    for ik in range(nk_sc):
        s_index = np.argwhere(egstates_ph_ik==ik)[:,0]
        kb_index[ik,:] = s_index
        # print(ik, s_index)

    # reorder = [3, 2, 1, 0, 4, 5, 6, 7]
    # kb_index = kb_index[:, reorder]
    # print(kb_index[0])
    return kb_index


if __name__=='__main__':
    main()
