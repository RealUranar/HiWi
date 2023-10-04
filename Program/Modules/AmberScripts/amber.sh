#!/usr/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
source ~/src/HiWi/.venv/bin/activate
conda activate AmberTools23

antechamber -fi gout -fo prepi -c resp -i orca.log -o orca.prep -rn F1 -at gaff2
parmchk2 -i orca.prep -f prepi -o orca.frcmod

python3.10 unitcell.py

PropPDB -p SHIFTED.PDB -o NEWPDB4x4.PDB -ix 1 -iy 4 -iz 4

tleap -f tleap.in

python3.10 convert2Gromax.py
rm SHIFTED.PDB