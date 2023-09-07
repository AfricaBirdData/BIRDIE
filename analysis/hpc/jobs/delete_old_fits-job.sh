#!/bin/sh

# This example submission script contains several important directives, please examine it thoroughly

# Do not put spaces between the start of the line and #SBATCH, the line must start exactly with #SBATCH, no spaces.
# Do not put spaces between the # and SBATCH

# The line below indicates which accounting group to log your job against
#SBATCH --account=biosci

# The line below selects the group of nodes you require
#SBATCH --partition=ada

# The line below reserves 1 worker node and 40 cores
#SBATCH --nodes=1 --ntasks=1

# The line below indicates the wall time your job will need, 10 hours for example. NB, this is a mandatory directive!
#SBATCH --time=00:05:00

# A sensible name for your job, try to keep it short
#SBATCH --job-name="brd_remove_fits"

#Modify the lines below for email alerts. Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-user=crvfra001@myuct.ac.za
#SBATCH --mail-type=ALL

# The cluster is configured primarily for OpenMPI and PMI. Use srun to launch parallel jobs if your code is parallel aware.
# To protect the cluster from code that uses shared memory and grabs all available cores the cluster has the following
# environment variable set by default: OMP_NUM_THREADS=1
# If you feel compelled to use OMP then uncomment the following line:
# export OMP_NUM_THREADS=$SLURM_NTASKS

# NB, for more information read https://computing.llnl.gov/linux/slurm/sbatch.html

# Use module to gain easy access to software, typing module avail lists all packages.

module load software/R-4.2.0 software/geos-3.10.2 software/gdal-2.4.2 software/jags-4.3.0

R CMD BATCH ~/birdie/scripts/delete_old_fits.R

