import numpy as np


def main():

    eqp_12K = np.loadtxt('eqp1.dat')

    nk = 144
    nb = 18
    eqp_12K = eqp_12K.reshape(nk, nb+1, 4)

    eqp_6K  = []
    for ik in range(nk):
        kpt = eqp_12K[ik, 0]
        if (int( round(kpt[0]*12) )%2 == 1 ): continue
        if (int( round(kpt[1]*12) )%2 == 1 ): continue
        eqp_6K.append(eqp_12K[ik])
    
    eqp_6K = np.array(eqp_6K).reshape(-1, 4)

    write_eqp(eqp_6K, filname='eqp1_6K.dat')


def write_eqp(eqp, filname='eqp1.dat'):
    
    nb = int(eqp[0, -1])
    nk = int(eqp.shape[0] / (nb + 1))
    eqp = eqp.reshape(nk, nb+1, 4)

    fm1 = "{:13.9f}{:13.9f}{:13.9f}{:8.0f}"
    fm2 = "{:8.0f}{:8.0f}{:15.9f}{:15.9f}"
    with open(filname, 'w') as f:

        for ik in range(nk):
            line = fm1.format(*eqp[ik, 0])
            f.write(line)
            f.write('\n')
            for ib in range(nb):
                line = fm2.format(*eqp[ik, ib+1])
                f.write(line)
                f.write('\n')




if __name__=='__main__':
    main()

