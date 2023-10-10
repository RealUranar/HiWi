import os, shutil, sys
sys.path.append("Modules/Misc")
from Sbatch import JobScripts
from task import Task

import subprocess

class GromacsProd(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def writeInputFile(self):
        shutil.copy("Modules/GromacsScripts/prod.mdp", f"{self.newPath}")
        shutil.copy("Modules/GromacsScripts/plumed.dat", f"{self.newPath}")
        dihedral = self._findNNDihedral(f"Gromacs/System.gro")

        with open(f"{self.newPath}/plumed.dat", "r") as file:
            lines = file.readlines()
        
        for line in lines:
            if "ATOMS=" in line:
                print(line)

    def generateJobScript(self):
        with open(f"{self.newPath}/prod.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f prod.mdp -c nvt.gro -r nvt.gro -p System.top -o prod.tpr\n",
            ])
        os.chmod(f"{self.newPath}/prod.sh", 0o755)
        JobScripts().writeGromacsJob(name = self.job.id, location=self.newPath)


    def submit(self):
        ret = subprocess.run(f"./prod.sh",
                    capture_output = True,
                    text = True,
                    cwd=self.newPath)
        print(f"Setup for Gromacs production job {self.job.name} finished with code {ret.returncode}")
        self.job.updateJob(GromacsProduction = 2)
        print(f"Submitted Gromacs job {self.job.name}")
        return super().submit(self.newPath)
        
        
    def isFinished(self):
        hasFinished, succesfull = False, False

        tail = self._readTail(self.newPath)
        if "Segmentation fault" in str(tail) or "DUE TO TIME LIMIT" in str(tail):
            hasFinished = True
        elif "Writing final coordinates." in str(tail): 
            hasFinished = True
            succesfull = True

        if hasFinished:
            if succesfull:
                self.job.updateJob(GromacsProduction = 1)
                print(f"Gromacs Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(GromacsProduction = -1)
                print(f"Gromacs Job {self.job.name} run into a problem")
        else:
            print(f"Gromacs Job {self.job.name} is still running")


if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})

    task = GromacsProd(job)
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
