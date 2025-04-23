import h5py
import numpy as np
import os

def main():

    ndim = 9
    fil_dwfns = ['dwfn%d.h5'%(idim+1) for idim in range(ndim)]
    fil_wfn = 'wfn.h5'
    displ = np.zeros(ndim)
    # displ = displ + 0.005
    displ[7] = 0.01
    fil_out = 'wfn_dR8_01_qG_norm.h5'

    # displ = np.loadtxt('displ_qG_m4-9.dat').flatten()
    # fil_out = 'wfn_qG_m4-9.h5'
    print(displ)

    constr_wfn_ph(fil_dwfns, fil_wfn, displ, fil_out)
    # check_wfns(fil_dwfns, fil_wfn, fil_out)


def constr_wfn_ph(fil_dwfns, fil_wfn, displ, fil_out='wfn_ph.h5'):

    # Copy file1 to output using os.system
    os.system(f'cp "{fil_wfn}" "{fil_out}"')  # Linux/Mac
    # os.system(f'cp "wfn_sc.h5" "{fil_out}"')  # Linux/Mac
    # For Windows use: os.system(f'copy "{fil_wfn}" "{fil_out}"')

    ndim = len(fil_dwfns)

    with h5py.File(fil_out, 'r+') as f_out:
        coeffs = f_out['/wfns/coeffs'][:]
        gvecs1 = f_out['/wfns/gvecs'][:]
        nkG = f_out['/mf_header/kpoints/ngk'][:]


        # kG_index = np.cumsum(nkG)
        # print(kG_index)
        kG_index = find_kG_index(gvecs1)
        # print(kG_index)
        
        for idim in range(ndim):

            fil_dwfn = fil_dwfns[idim]
            with h5py.File(fil_dwfn, 'r') as f:
                dcoeffs = f['/wfns/coeffs'][:]
                gvecs2 = f['/wfns/gvecs'][:]

            if coeffs.shape != dcoeffs.shape:
                raise ValueError("Error: /wfns/coeffs datasets have different shapes!")

            if (not np.array_equal(gvecs1, gvecs2)):
                raise ValueError("Error: /wfns/gvecs mismatch!")

            coeffs = coeffs + dcoeffs * displ[idim]

        coeffs = normalize_coeffs(coeffs, kG_index)
        # print(coeffs.shape)

        del f_out['/wfns/coeffs']
        f_out.create_dataset('/wfns/coeffs', data=coeffs, dtype=coeffs.dtype)


def normalize_coeffs(coeffs, kG_index):

    coeffs_complex = coeffs[..., 0] + 1j * coeffs[..., 1]

    nk = kG_index.shape[0] - 1

    for ik in range(nk):

        ikG = kG_index[ik]
        jkG = kG_index[ik+1]

        norm_factor = np.sum(np.abs(coeffs_complex[..., ikG:jkG])**2, axis=2, keepdims=True)
        norm_factor[norm_factor == 0] = 1.0
        print(np.sum(norm_factor))
        coeffs_complex[..., ikG:jkG] = coeffs_complex[..., ikG:jkG] / np.sqrt(norm_factor)
    
    coeffs_normalized = coeffs
    coeffs_normalized[..., 0] = coeffs_complex.real
    coeffs_normalized[..., 1] = coeffs_complex.imag

    return coeffs_normalized


def find_kG_index(gvecs):

    nkG = gvecs.shape[0]

    kG_index = []
    for ikG in range(nkG):
        if (np.sum(np.abs(gvecs[ikG,:]))==0):
            kG_index.append(ikG)
    kG_index.append(nkG)

    return np.array(kG_index)


def check_wfns(fil_dwfns, fil_wfn, fil_out='wfn_ph.h5'):

    # ib, ispin, iG, real/imag
    index = tuple((32, 0, 0))
    
    ndim = len(fil_dwfns)

    for idim in range(ndim):
        fil_dwfn = fil_dwfns[idim]
        with h5py.File(fil_dwfn, 'r') as f:
            dcoeffs = f['/wfns/coeffs'][:]
            print(idim+1, dcoeffs[index])

    with h5py.File(fil_wfn, 'r') as f:
        coeffs = f['/wfns/coeffs'][:]
        print(coeffs[index])

    with h5py.File(fil_out, 'r') as f:
        coeffs = f['/wfns/coeffs'][:]
        print(coeffs[index])


if __name__=='__main__':
    main()

