#! /usr/bin/env bash
#
#SBATCH --job-name=internal_winter
#
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G

python3 internal_pager_master_winter.py