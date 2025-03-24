import h5py
import numpy as np
import os

def main():

    ndim = 9
    fil_dwfns = ['dwfn%d.h5'%(idim+1) for idim in range(ndim)]
    fil_wfn = 'wfn.h5'
    displ = np.zeros(ndim)
    displ = 1.0
    fil_out = 'wfn_dR1_qG.h5'

    bohr = 0.529177210903 # Angstrom
    displ = np.loadtxt('displ_qG_m4-9.dat').flatten() / bohr
    fil_out = 'wfn_qG_m4-9.h5'

    constr_wfn_ph(fil_dwfns, fil_wfn, displ, fil_out)
    # check_wfns(fil_dwfns, fil_wfn, fil_out)


def constr_wfn_ph(fil_dwfns, fil_wfn, displ, fil_out='wfn_ph.h5'):

    # Copy file1 to output using os.system
    os.system(f'cp "{fil_wfn}" "{fil_out}"')  # Linux/Mac
    # For Windows use: os.system(f'copy "{fil_wfn}" "{fil_out}"')

    ndim = len(fil_dwfns)

    with h5py.File(fil_out, 'r+') as f_out:
        coeffs = f_out['/wfns/coeffs'][:]

        for idim in range(ndim):

            fil_dwfn = fil_dwfns[idim]
            with h5py.File(fil_dwfn, 'r') as f:
                dcoeffs = f['/wfns/coeffs'][:]

            if coeffs.shape != dcoeffs.shape:
                raise ValueError("Error: /wfns/coeffs datasets have different shapes!")

            coeffs = coeffs + dcoeffs

        del f_out['/wfns/coeffs']
        f_out.create_dataset('/wfns/coeffs', data=coeffs, dtype=coeffs.dtype)


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

