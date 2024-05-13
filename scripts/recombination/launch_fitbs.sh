#!/bin/bash
#$ -S /bin/bash
#$ -M sylvain.glemin@univ-rennes.fr
#$ -m bea
#SBATCH --job-name=fitbs
#SBATCH --chdir=workingdirectory
#SBATCH --nodes=4 
#SBATCH --time=48:00

/local/env/envr-4.1.3.sh

Rscript 6_linkedselection_fitmodel.R

exit 0