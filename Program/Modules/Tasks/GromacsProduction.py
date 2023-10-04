import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import numpy as np
from Sbatch import JobScripts
from task import Task
from excel import Excel
import subprocess

class GromacsProd(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"


    def writeInputFile(self):
        shutil.copy("Modules/GromacsScripts/prod.mdp", f"{self.newPath}")
        shutil.copy("Modules/GromacsScripts/plumed.dat", f"{self.newPath}")


    def generateJobScript(self):
        with open(f"{self.newPath}/prod.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f prod.mdp -c nvt.gro -r nvt.gro -p System.top -o prod.tpr\n",
            ])
        os.chmod(f"{self.newPath}/prod.sh", 0o755)
        JobScripts().writeGausianJob(name = self.job.name, location=self.newPath)



    def submit(self):
        ret = subprocess.run(f"./prod.sh",
                    capture_output = True, 
                    text = True,
                    cwd=self.newPath)
        print(ret)
        return super().submit(self.newPath)
        

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = GromacsProd(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
