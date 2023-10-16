import os, shutil, sys
from task import Task

sys.path.append("Modules/Misc")
from Sbatch import JobScripts
import subprocess

class GromacsEquill(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]
        
    def writeInputFile(self, temp =310):
        shutil.copy("Modules/GromacsScripts/nvt.mdp", f"{self.newPath}")

    def generateJobScript(self):
        with open(f"{self.newPath}/nvt.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f nvt.mdp -c em.gro -p System.top -r em.gro -o nvt.tpr",
            ])
        os.chmod(f"{self.newPath}/nvt.sh", 0o755)
        JobScripts().writeGromacsJob(name = self.job.id, location=self.newPath,  inputFile= "nvt.tpr",plumed=False, jobtype="nvt", time= "0-01:00:00")


    def submit(self):
        ret = subprocess.run(f"./nvt.sh",
                    capture_output = True, 
                    text = True,
                    cwd=self.newPath)
        print(f"Gromacs equillibration returned code: {ret.returncode}")
        if ret.returncode != 0:
             self.job.updateJob(GromacsEquil = -1)
        else:
            self.job.updateJob(GromacsEquil = 2)
            print(f"Submitted Gromacs Equillibration job {self.job.name}")
            return super().submit(self.newPath)
        
    def isFinished(self):
        tail = self._readTail(self.newPath, file= "nvt.log")
        hasFinished = "Finished mdrun" in str(tail)
        succesfull = "Finished mdrun" in str(tail)

        if hasFinished:
            if succesfull:
                self.job.updateJob(GromacsEquil = 1, GromacsProduction= 3)
                print(f"Gromacs Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(GromacsEquil = -1)
                print(f"Gromacs Job {self.job.name} run into a problem")
        else:
            print(f"Gromacs Job {self.job.name} is still running")

if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Modules/TESTING/", tasks={"Amber":1})

    task = GromacsEquill(job)
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
