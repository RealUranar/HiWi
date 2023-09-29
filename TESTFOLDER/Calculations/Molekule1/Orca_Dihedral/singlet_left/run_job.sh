#!/usr/bin/zsh

########## start of batch directives #######
#SBATCH --job-name=Molekule1
#SBATCH --output=output.%J.txt
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --constraint=Rocky8
#SBATCH --mem-per-cpu=3900
#SBATCH --partition=c18m
#SBATCH --account=p0020506
###### start of shell commands ######

module load GCC/11.3.0 OpenMPI/4.1.4  ORCA/5.0.4

/cvmfs/software.hpc.rwth.de/Linux/RH8/x86_64/intel/skylake_avx512/software/ORCA/5.0.4-gompi-2022a/bin/orca  orca.inp