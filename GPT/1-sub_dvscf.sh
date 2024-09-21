#!/bin/bash
# #SBATCH --nodes=1
# #SBATCH --time=02:00:00
# #SBATCH --constraint=cpu
# #SBATCH --qos=interactive
# #SBATCH --account=m2651

module load PrgEnv-gnu cray-libsci cray-fftw cray-hdf5-parallel cray-parallel-netcdf
ABINIT_PATH=/global/homes/z/zfzheng/src/abinit/abinit-modified/build/src/98_main
# NPROC=$SLURM_NTASKS
NPROC=512

qfile='./q.list'

qstart=1
nq=$(wc -l < $qfile)

for iq in `seq $qstart $nq`
do
 
    PATHQ="phonon-q-$iq/Dvscf"

    # Skip this directory if job is running or has ended
    if [[ -f $PATHQ/RUNNING || -f $PATHQ/ENDED ]]; then
        # Do nothing, go back to the parent directory
        continue
    fi

    phonq1=$(awk -v iq=$iq 'NR == iq {print $1}' $qfile)
    phonq2=$(awk -v iq=$iq 'NR == iq {print $2}' $qfile)
    phonq3=$(awk -v iq=$iq 'NR == iq {print $3}' $qfile)

    echo
    echo "================================================================================"
    echo "Running in the directory $PATHQ."
    echo "Q-point: ($phonq1, $phonq2, $phonq3)"
    echo "================================================================================"
    echo
 
    mkdir -p $PATHQ
    cp -a prep-phonon-q/Dvscf/* $PATHQ/
 
    cd $PATHQ

    sed -i "s/_phonq1_/$phonq1/g" calc.in
    sed -i "s/_phonq2_/$phonq2/g" calc.in
    sed -i "s/_phonq3_/$phonq3/g" calc.in
 
    touch RUNNING
 
    srun -n $NPROC $ABINIT_PATH/abinit < calc.files
 
    sleep 1
 
    # Check if jobs ended correctly
    if grep 'Calculation completed' calc.out >& /dev/null; then
        # If so, mark the job as ENDED
        touch ENDED
    else
        rm ENDED 2> /dev/null
    fi

    rm -f RUNNING
 
    cd ../..
done

