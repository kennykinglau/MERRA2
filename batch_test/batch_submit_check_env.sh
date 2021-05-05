#! /usr/bin/env bash
#
#SBATCH --job-name=env_test
#
#SBATCH --time=00:03:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

python3 -c "import matplotlib; import numpy; import scipy; import netCDF4; import pandas"

