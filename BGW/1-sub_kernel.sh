#!/usr/bin/env bash
# #SBATCH --nodes=1
# #SBATCH --time=02:00:00
# #SBATCH --constraint=cpu
# #SBATCH --qos=interactive
# #SBATCH --account=m2651

export CRAYPE_LINK_TYPE=static
module load cray-fftw cray-hdf5-parallel cray-libsci python
BGW_PATH=/global/homes/z/zfzheng/src/BerkeleyGW/bin

EXE=$BGW_PATH/kernel.cplx.x
INFILE='kernel.inp'
OUTFILE='kernel.out'

qfile='./q.list'
# NPROC=$SLURM_NTASKS
NPROC=512

qstart=1
nq=$(wc -l < $qfile)
nq=2

for iq in `seq $qstart $nq`
do
 
    PATHQ="bse-q-$iq/kernel"

    # Skip this directory if job is running or has ended
    if [[ -f $PATHQ/RUNNING || -f $PATHQ/ENDED ]]; then
        # Do nothing, go back to the parent directory
        continue
    fi

    Q1=$(awk -v iq=$iq 'NR == iq {print $1}' $qfile)
    Q2=$(awk -v iq=$iq 'NR == iq {print $2}' $qfile)
    Q3=$(awk -v iq=$iq 'NR == iq {print $3}' $qfile)

    echo
    echo "================================================================================"
    echo "Running in the directory $PATHQ."
    echo "Q-point: ($Q1, $Q2, $Q3)"
    echo "================================================================================"
    echo
 
    mkdir -p $PATHQ
    cp -a prep-bse-q/kernel/* $PATHQ/
 
    cd $PATHQ

    sed -i "s/_Q1_/$Q1/g" kernel.inp
    sed -i "s/_Q2_/$Q2/g" kernel.inp
    sed -i "s/_Q3_/$Q3/g" kernel.inp
 
    if [ $iq -eq 1 ]; then
        sed -i "s/_comm_tag_/# /g" kernel.inp
    else
        sed -i "s/_comm_tag_//g" kernel.inp
    fi
 
    touch RUNNING
 
    srun -n $NPROC $EXE $INFILE > $OUTFILE
 
    sleep 1
 
    # Check if jobs ended correctly
    if grep 'TOTAL' kernel.out >& /dev/null; then
        # If so, mark the job as ENDED
        touch ENDED
    else
        rm ENDED 2> /dev/null
    fi

    rm -f RUNNING
 
    cd ../..
done

