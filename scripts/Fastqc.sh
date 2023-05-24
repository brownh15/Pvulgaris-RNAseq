#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########
 
#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=2G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name fastqc      # you can give your job a name for easier identification (same as -J)
#SBATCH --ntasks-per-node=1
 
########## Command Lines to Run ##########
 
module load FastQC/0.11.7-Java-1.8.0_162
 
cd ${PBS_O_WORKDIR}

fastqc -o FastQC/ E11B_S28_L004_R1_001_trimmed.fastq.gz