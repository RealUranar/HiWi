#!/usr/bin/zsh

########## start of batch directives #######
#SBATCH --job-name=Molekule1
#SBATCH --output=GAUSSIANJOB.%J
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M

###### start of shell commands ######

module load Gaussian/16.C.01-AVX2

g16 < orca.com > gauss.log