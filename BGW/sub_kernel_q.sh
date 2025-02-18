#!/bin/bash
#SBATCH -J kernel
#SBATCH -N 3
#SBATCH --ntasks-per-node 56
#SBATCH -p development
#SBATCH -t 01:45:00
#SBATCH -A DMR22042


module load arpack impi intel/19.1.1 phdf5
export PATH=/home1/09917/zfzheng/software/BerkeleyGW/bin:$PATH

JOB=kernel
JOB_FOLDER=kernel
EXE=kernel.cplx.x

QSTART=1
NQ=$(wc -l < q.list)

PATH_WFK='../Wfk/'
PATH_EPSILON='../gw/1-epsilon/'
PATH_PREPARE='prep_bse-q/'

PATH_WFK=$(cd $PATH_WFK && pwd)
PATH_EPSILON=$(cd $PATH_EPSILON && pwd)
PATH_PREPARE=$(cd $PATH_PREPARE && pwd)

NNODES_PER_Q=$(($SLURM_JOB_NUM_NODES / $NQ))
if (($NNODES_PER_Q == 0)); then
  NNODES_PER_Q=1
fi

NTASKS_PER_Q=$(($NNODES_PER_Q * $SLURM_NTASKS_PER_NODE))
NQ_PER_LOOP=$(($SLURM_JOB_NUM_NODES / $NNODES_PER_Q))

echo -e "\nThe number of q-points: $NQ"
echo

for iq in `seq $QSTART $NQ`
do

  DIRQ=bse-q-$iq
  if [ ! -d "$DIRQ" ]; then
    mkdir $DIRQ
  fi
  cd $DIRQ

  if [ ! -d "$JOB_FOLDER" ]; then
    mkdir $JOB_FOLDER
  fi
  cd $JOB_FOLDER

  OFFSET=$((($iq - 1) * $NTASKS_PER_Q))

  # Copy input files to node local /tmp directory, 144 GB per node
  # ibrun -n 1 -o $OFFSET cp -L $PATH_WFK/wfn.cplx /tmp/WFN
  # ln -sf /tmp/WFN WFN_co

  cp -L $PATH_WFK/wfn.cplx WFN_co
  cp -L $PATH_WFK/wfn.cplx WFNq_co
  cp -L $PATH_EPSILON/eps0mat.h5 .
  cp -L $PATH_EPSILON/epsmat.h5 .

  cp $PATH_PREPARE/$JOB_FOLDER/${JOB}.inp .

  Q1=$(awk -v iq=$iq 'NR == iq {print $1}' ../../q.list)
  Q2=$(awk -v iq=$iq 'NR == iq {print $2}' ../../q.list)
  Q3=$(awk -v iq=$iq 'NR == iq {print $3}' ../../q.list)
  sed -i "s/_Q1_/$Q1/g" ${JOB}.inp
  sed -i "s/_Q2_/$Q2/g" ${JOB}.inp
  sed -i "s/_Q3_/$Q3/g" ${JOB}.inp
  echo "Calculating at q-point: " $Q1 $Q2 $Q3
  echo

  ibrun -n $NTASKS_PER_Q -o $OFFSET task_affinity $EXE < $JOB.inp > $JOB.out &

  cd ../..

  if (( $iq % $NQ_PER_LOOP == 0 )); then
    wait
  fi

done

wait

echo "All tasks have completed."
