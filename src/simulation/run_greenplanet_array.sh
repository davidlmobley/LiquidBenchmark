#!/bin/bash
#
# Job name:
#----------------
#SBATCH -J "dens_array"
#----------------

# Array info
#SBATCH --array=1-247%5

#----------------
#SBATCH -p mf_titanx
#----------------

#----------------
#SBATCH -o out.txt
#----------------

#----------------
#SBATCH -e err.txt
#----------------

# Specifying resources needed for run:
#
#--------------
#SBATCH --nodes=1
# #SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=72:00:00
#SBATCH --distribution=block:cyclic
#SBATCH --partition=mf_titanx
#SBATCH --gres=gpu:titanxa:1
# If you want to specify a spectific gpu partiion use #SBATCH --gres=gpu:titanxb:1
# or #SBATCH --gres=gpu:titanxa:1

#--------------

# Informational output
echo "=================================== SLURM JOB ==================================="
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             $SLURM_JOB_ID"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST"
echo "Submit directory:   $SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo

cd $SLURM_SUBMIT_DIR

echo 'Working Directory:'
pwd
echo 'Array ID number:' $SLURM_ARRAY_TASK_ID
 
date
#insert submitted commands here
python simulate_thermoml.py $SLURM_ARRAY_TASK_ID

###
date
