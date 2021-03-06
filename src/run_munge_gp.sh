#!/bin/bash
#
# Author: Gaetano Calabro, UCI gacabro@uci.edu
# Job name:
#----------------
#SBATCH -J "munge3"
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
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=12:00:00
#SBATCH --distribution=block:cyclic
#SBATCH --partition=mf_m-c1.9

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

 
date
# or #SBATCH --gres=gpu:titanxa:1
#insert submitted commands here
python munge_output_amber.py

###
date
