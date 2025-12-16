#!/bin/bash
#SBATCH -A che240225         # Allocation name
#SBATCH -p shared
#SBATCH -t 12:00:00              # Time limit
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 40                    # Total tasks
#SBATCH -c 1                     # Cores per task
#SBATCH --job-name=alo2.4
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out
#SBATCH --hint=nomultithread

# Load environment
module purge
module load gcc/11.2.0 openmpi/4.1.6
module load quantum-espresso/7.3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running pw.x (SCF)..."
srun pw.x -inp topw.al5.scf.in -npool 10 -ndiag 4 > topw.al5.30.scf.out

echo "job complete"
