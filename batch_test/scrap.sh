#!/bin/env bash
#
#SBATCH --job-name=SPconsistency_job
#
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

srun python SPconsistency_check.py -fb MERRA2_400.inst1_2d_int_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.inst3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg1_2d_slv_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31


#!/bin/env bash
#
#SBATCH --job-name=SPconsistency_job
#
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

module load python/3.6
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r requirements.txt

srun python SPconsistency_check.py -fb MERRA2_400.inst1_2d_int_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.inst3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg1_2d_slv_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31

# Note that requirements.txt contains:
#py-h5py
#py-pandas
#py-jupyter
#py-numpy
#py-netcdf4
#py-matplotlib
#py-scipy

# Then I tried adding the "ml load" statements from my .bashrc:
#!/bin/env bash
#
#SBATCH --job-name=SPconsistency_job
#
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

ml load python/3.6.1
ml load viz
ml load py-h5py/2.10.0_py36
ml load py-pandas/1.0.3_py36
ml load py-jupyter/1.0.0_py36
ml load py-numpy/1.18.1_py36
ml load py-netcdf4/1.3.1_py36
ml load py-matplotlib/3.2.1_py36
ml load py-scipy/1.4.1_py36

srun python SPconsistency_check.py -fb MERRA2_400.inst1_2d_int_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.inst3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg1_2d_slv_Nx. -fe .SUB.nc -fn SPconsistencyAll2019 -v TQI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31
srun python SPconsistency_check.py -fb MERRA2_400.tavg3_3d_asm_Nv. -fe .SUB.nc -fn SPconsistencyAll2019 -v QI -mi -sy 2019 -sm 1 -sd 1 -ey 2019 -em 12 -ed 31

