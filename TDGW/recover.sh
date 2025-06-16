#!/bin/bash

NQ=3

# recover the data files.
unlink epc.h5
mv epc_backup.h5 epc.h5
for iq in `seq 1 $NQ`
do
    cd bse_$iq
    unlink kernel.h5
    mv kernel_${iq}.h5 kernel.h5
    cd ..
done

