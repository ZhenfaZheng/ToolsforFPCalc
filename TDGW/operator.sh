#!/bin/bash

prefix='test'
prefix='3Q_eh_noph_pulse_light'

folder="results/1b_24K_${prefix}"
sub_folder="1b_24K_${prefix}"

nq=3
lrestart=0

if [ "$lrestart" -eq 1 ]; then
    echo 'Copy data from backup files.'
    cp $folder/dynamics.inp .
    cd bse_1
    cp occupation_${prefix}.h5 occupation.h5
    cp tdinfo_${prefix}.h5 tdinfo.h5
    cd ..
    for iq in `seq 1 $nq`
    do
        cd bse_$iq
        cp occupation_${prefix}.h5 occupation.h5
        cp polarization_${prefix}.dat polarization.dat
        cd ..
    done
fi

# python3 tdagw.py

cd bse_1
cp occupation.h5 occupation_${prefix}.h5
cp tdinfo.h5 tdinfo_${prefix}.h5
cd ..
for iq in `seq 1 $nq`
do
    cd bse_$iq
    cp polarization.dat polarization_${prefix}.dat
    cd ..
done

mkdir -p $folder
mv bse_*/*png $folder/
mv PE.png PolAbs.png TDBAND.png TDEN.png TDPOP.png TDPOP_FFT.png $folder
mv TRPES*png $folder
# cp *gif $folder
cp tdoccsum.dat $folder/
cp dynamics.inp $folder/

mkdir -p $folder/codes
cp ~/software/BerkeleyGW/TDGW/*f90 $folder/codes


source ~/.bashrc
cd results
send $sub_folder/*
cd ..
