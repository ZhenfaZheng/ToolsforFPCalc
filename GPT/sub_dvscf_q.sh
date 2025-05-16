#!/bin/bash
#SBATCH -J Dvscf
#SBATCH -N 36
#SBATCH --ntasks-per-node 56
#SBATCH -p flex
#SBATCH -t 24:00:00
##SBATCH -p development
##SBATCH -t 02:00:00
#SBATCH -A DMR22042

module load parallel-netcdf/4.6.2
export PATH=/home1/09917/zfzheng/software/abinit-modified/src/98_main:$PATH
EXE=abinit

QSTART=49
NQ=$(wc -l < q.list)

NNODES_PER_Q=$(($SLURM_JOB_NUM_NODES / $NQ))
if (($NNODES_PER_Q == 0)); then
  NNODES_PER_Q=1
fi
NNODES_PER_Q=3

NTASKS_PER_Q=$(($NNODES_PER_Q * $SLURM_NTASKS_PER_NODE))
NQ_PER_LOOP=$(($SLURM_JOB_NUM_NODES / $NNODES_PER_Q))


for iq in `seq $QSTART $NQ`
do

  DIRQ=phonon-q-$iq
  cd $DIRQ/Dvscf

  # Copy input files to node local /tmp directory, 144 GB per node
  # ibrun -n 1 -o $OFFSET cp -L $PATH_WFK/wfn.cplx /tmp/WFN
  # ln -sf /tmp/WFN WFN_co

  Q1=$(awk -v iq=$iq 'NR == iq {print $1}' ../../q.list)
  Q2=$(awk -v iq=$iq 'NR == iq {print $2}' ../../q.list)
  Q3=$(awk -v iq=$iq 'NR == iq {print $3}' ../../q.list)
  echo "Nq = " $iq
  echo "Calculating at q-point: " $Q1 $Q2 $Q3
  echo

  OFFSET=$((($iq - 1) % $NQ_PER_LOOP * $NTASKS_PER_Q))

  ibrun -n $NTASKS_PER_Q -o $OFFSET task_affinity $EXE < calc.files &> calc.log 2> calc.err &

  cd ../..

  sleep 5

  if (( $iq % $NQ_PER_LOOP == 0 )); then
    wait
  fi

done

wait

echo "All tasks have completed."

