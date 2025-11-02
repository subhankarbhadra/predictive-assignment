#!/bin/bash

#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=20       # Number of CPU cores per task
#SBATCH --mem-per-cpu=16G        # Memory per CPU
#SBATCH --time=6:00:00           # Time limit
#SBATCH --partition=open         # Partition name
#SBATCH --output=out.%J          # Standard output (J is job ID)
#SBATCH --error=err.%J           # Standard error (J is job ID)

module load matlab               # Load the MATLAB module

# Run your .m file in MATLAB
matlab -nodisplay -r "run('newsim3.m'); exit;"

