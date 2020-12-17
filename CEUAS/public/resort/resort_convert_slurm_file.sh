#!/bin/bash
#SBATCH --job-name=CEUAS	     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ulrich.voggenberger@univie.ac.at    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # Number of processes
#SBATCH --mem=100gb                    # Total memory limit
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=CEUAS_%j.log # Standard output and error log
module load anaconda3
python convert.py
