import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import subprocess
from Sbatch import JobScripts
from task import Task
from excel import Excel


class Orca_opt(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Amber"

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        gaussOutFile = glob.glob(f"{self.job.location}Gaussian/*.log")[0]  #Get Path to the log file
        shutil.copy(gaussOutFile, f"{self.newPath}/") #Copy log File to new directory

    def things(self):
        console = subprocess.Popen(["source ~/miniconda3/etc/profile.d/conda.sh"],
                        shell=True,
                        stdin =subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True,
                        bufsize=0)
        console.stdin.write("source ~/miniconda3/etc/profile.d/conda.sh")
        console.stdin.write("source ~/src/.venv/bin/activate")
        console.stdin.write("conda activate AmberTools23")
        console.stdin.close()
        for line in console.stdout:
            print(line.strip())


if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Orca_opt(jobs[0])
    task.things()


#!/usr/bin/bash

# source ~/miniconda3/etc/profile.d/conda.sh
# source ~/src/.venv/bin/activate
# conda activate AmberTools23

# antechamber -fi gout -fo prepi -c resp -i orca.log -o orca.prep -rn F1 -at gaff2
# parmchk2 -i orca.prep -f prepi -o orca.frcmod

# python3.10 unitcell.py

# PropPDB -p SHIFTED.PDB -o NEWPDB4x4.PDB -ix 1 -iy 4 -iz 4

# tleap -f tleap.in

# python3.10 convert2Gromax.py
# rm SHIFTED.PDB